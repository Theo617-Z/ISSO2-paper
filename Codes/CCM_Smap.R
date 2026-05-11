required_packages <- c("dplyr", "lubridate", "rEDM")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(rEDM)
})

script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script <- sub(file_arg, "", cmd[grepl(file_arg, cmd)])
  if (length(script) > 0) {
    return(dirname(normalizePath(script[1])))
  }
  getwd()
}

base_dir <- dirname(script_dir())
data_file <- file.path(base_dir, "Data", "analysis_data.csv")

df_analysis <- read.csv(data_file)
dis_name <- "isc"
plt_name <- c("fsp", "o3", "so2", "temp", "rh")
lag_days <- 0:3
E_value <- 9
num_surr <- as.integer(Sys.getenv("ISSO2_CCM_SURROGATES", "500"))
num_sample <- 100
seed_value <- 2019
theta_grid <- c(0, 0.01, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9)

run_ccm_task <- function(plt, lag_day) {
  pair_df <- df_analysis %>% select(date, all_of(dis_name), all_of(plt))
  lib_size <- nrow(pair_df) - E_value - 5

  set.seed(seed_value)
  surr_data <- SurrogateData(
    unlist(pair_df[, plt]),
    method = "ebisuzaki",
    num_surr = num_surr,
    alpha = 1
  )

  all_data <- as.data.frame(cbind(pair_df, surr_data))
  names(all_data) <- c("date", dis_name, "T1", paste0("T", 2:(num_surr + 1)))

  out <- data.frame(
    i = seq_len(num_surr + 1),
    rho = NA_real_,
    dis = dis_name,
    plt = plt,
    E = E_value,
    lag_day = lag_day
  )

  for (i in seq_len(num_surr + 1)) {
    target_col <- paste0("T", i)
    ccm_result <- CCM(
      dataFrame = all_data,
      E = E_value,
      Tp = -lag_day,
      columns = dis_name,
      target = target_col,
      libSizes = lib_size,
      random = FALSE,
      sample = num_sample,
      seed = seed_value
    )
    out$rho[i] <- ccm_result[1, paste0(dis_name, ":", target_col)]
  }

  out
}

summarise_ccm <- function(ccm_df) {
  ccm_df %>%
    group_by(dis, plt, E, lag_day) %>%
    summarise(
      rho_obs = rho[i == 1][1],
      n_surrogates = sum(i != 1),
      n_lower = sum(rho[i != 1] <= rho[i == 1][1], na.rm = TRUE),
      n_upper = sum(rho[i != 1] >= rho[i == 1][1], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      p_lower = (n_lower + 1) / (n_surrogates + 1),
      p_upper = (n_upper + 1) / (n_surrogates + 1),
      p_value = pmin(1, 2 * pmin(p_lower, p_upper)),
      p_tail = ifelse(p_lower <= p_upper, "lower", "upper"),
      significant = p_value <= 0.05 + sqrt(.Machine$double.eps)
    ) %>%
    select(-n_lower, -n_upper) %>%
    arrange(match(plt, plt_name), lag_day)
}

theta_out <- bind_rows(lapply(theta_grid, function(theta) {
  y_df <- df_analysis %>% select(date, all_of(dis_name))
  PredictNonlinear(
    dataFrame = y_df,
    embedded = FALSE,
    columns = dis_name,
    target = dis_name,
    Tp = 1,
    theta = theta,
    lib = c(1, nrow(y_df)),
    pred = c(1, nrow(y_df)),
    showPlot = FALSE,
    E = E_value
  ) %>%
    mutate(dis = dis_name, theta = theta, E = E_value)
}))

best_theta <- theta_out %>%
  slice_max(order_by = rho, n = 1, with_ties = FALSE) %>%
  pull(theta)

as_safe_date <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }
  if (is.numeric(x)) {
    return(lubridate::as_date(x, origin = lubridate::origin))
  }
  as.Date(x)
}

run_smap_task <- function(plt, lag_day) {
  lagged_df <- df_analysis %>%
    mutate(plt_lag = dplyr::lag(.data[[plt]], lag_day)) %>%
    select(date, all_of(dis_name), plt_lag)

  embedded_y <- Embed(dataFrame = lagged_df, E = E_value, tau = -1, columns = dis_name)
  past_y <- embedded_y[, 2:E_value, drop = FALSE]
  names(past_y) <- paste0(dis_name, "_lag", seq_len(E_value - 1))

  exposure_df <- data.frame(lagged_df$plt_lag)
  names(exposure_df) <- plt

  smap_df <- cbind(
    lagged_df[, c("date", dis_name)],
    as.data.frame(past_y),
    exposure_df
  ) %>%
    filter(complete.cases(.)) %>%
    as.data.frame()

  smap_columns <- paste(c(names(past_y), plt), collapse = " ")
  smap_result <- SMap(
    dataFrame = smap_df,
    embedded = TRUE,
    columns = smap_columns,
    target = dis_name,
    lib = c(1, nrow(smap_df)),
    pred = c(1, nrow(smap_df)),
    theta = best_theta,
    Tp = 0,
    E = E_value
  )

  smap_result$coefficients[c(1, 2 + E_value)] %>%
    setNames(c("date", "effect")) %>%
    mutate(
      date = as_safe_date(date),
      dis = dis_name,
      plt = plt,
      E = E_value,
      lag_day = lag_day,
      theta = best_theta
    )
}

ccm_raw <- bind_rows(lapply(plt_name, function(plt) {
  bind_rows(lapply(lag_days, function(lag_day) run_ccm_task(plt, lag_day)))
}))
ccm_res <- summarise_ccm(ccm_raw)

smap_raw <- bind_rows(lapply(plt_name, function(plt) {
  bind_rows(lapply(lag_days, function(lag_day) run_smap_task(plt, lag_day)))
}))
smap_res <- smap_raw %>%
  group_by(dis, plt, E, lag_day, theta) %>%
  summarise(effect = mean(effect, na.rm = TRUE), n_effect_dates = sum(!is.na(effect)), .groups = "drop") %>%
  arrange(match(plt, plt_name), lag_day)

write.csv(ccm_res, file.path(base_dir, "ccm_res.csv"), row.names = FALSE)
write.csv(smap_res, file.path(base_dir, "smap_res.csv"), row.names = FALSE)
message("Wrote ccm_res.csv and smap_res.csv")
