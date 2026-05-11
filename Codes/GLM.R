required_packages <- c("dplyr", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages(library(dplyr))

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
out_file <- file.path(base_dir, "glm_res.csv")

data <- read.csv(data_file) %>%
  select(date, dow, holiday, isc, temp, rh, o3, fsp, so2) %>%
  arrange(as.Date(date)) %>%
  mutate(
    across(
      c(o3, fsp, temp, rh, so2),
      list(
        lag0 = ~ dplyr::lag(.x, 0),
        lag1 = ~ dplyr::lag(.x, 1),
        lag2 = ~ dplyr::lag(.x, 2),
        lag3 = ~ dplyr::lag(.x, 3)
      ),
      .names = "{.col}_{.fn}"
    )
  )

fit_lag_glm <- function(data, plt, lag_day) {
  lag_terms <- paste0(plt, "_lag", c(lag_day, setdiff(0:3, lag_day)))
  extra_terms <- character(0)

  if (plt == "fsp" && lag_day == 2) {
    extra_terms <- c(extra_terms, "so2_lag2")
  }
  if (plt == "so2" && lag_day == 0) {
    extra_terms <- c(extra_terms, "temp")
  }

  formula_text <- paste(
    "isc ~",
    paste(c(lag_terms, extra_terms, "factor(dow)", "factor(holiday)"), collapse = " + ")
  )

  fit <- glm(
    formula = as.formula(formula_text),
    data = data,
    family = gaussian(),
    na.action = na.omit
  )

  term <- paste0(plt, "_lag", lag_day)
  coef_tab <- summary(fit)$coefficients

  data.frame(
    beta = coef_tab[term, 1],
    se = coef_tab[term, 2],
    t = coef_tab[term, 3],
    p = coef_tab[term, 4],
    plt = plt,
    lag_day = lag_day
  )
}

glm_res <- bind_rows(lapply(
  c("temp", "rh", "o3", "fsp", "so2"),
  function(plt) bind_rows(lapply(0:3, function(lag_day) fit_lag_glm(data, plt, lag_day)))
)) %>%
  mutate(
    effectlow = beta - 1.96 * se,
    effecthigh = beta + 1.96 * se,
    sig = ifelse(p < 0.05, "sig", "not_sig")
  )

iqr_tbl <- data %>%
  select(-c(date, dow, holiday, isc), -matches("_lag[0-3]$")) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "plt", values_to = "value") %>%
  group_by(plt) %>%
  summarise(iqr = IQR(value, na.rm = TRUE), .groups = "drop")

glm_res <- glm_res %>%
  left_join(iqr_tbl, by = "plt") %>%
  transmute(
    plt,
    lag_day,
    effect = beta * iqr,
    effectlow = effectlow * iqr,
    effecthigh = effecthigh * iqr,
    p,
    sig,
    method = "GLM"
  ) %>%
  arrange(match(plt, c("fsp", "o3", "so2", "temp", "rh")), lag_day)

write.csv(glm_res, out_file, row.names = FALSE)
message("Wrote ", out_file)
