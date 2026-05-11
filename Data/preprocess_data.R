required_packages <- c("dplyr", "lubridate")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
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

data_dir <- script_dir()
raw_file <- file.path(data_dir, "hko_23ap_ha_0516.rda")
out_file <- file.path(data_dir, "analysis_data.csv")

load(raw_file)

normalize_series <- function(x, normalization = TRUE, dtrend = TRUE,
                             dseason = TRUE, season_sd = FALSE, sea = 365,
                             dTtype = "linear") {
  x <- as.numeric(x)
  xt <- x

  if (dtrend && dTtype == "first") {
    xt <- diff(xt)
  } else if (dtrend && dTtype == "linear") {
    fit <- lm(xt ~ seq_along(xt))
    xt <- xt - (fit$coefficients[1] + fit$coefficients[2] * seq_along(xt))
  }

  if (dseason) {
    full_year_length <- sea * length(xt) %/% sea
    seasonal_matrix <- matrix(xt[seq_len(full_year_length)], ncol = sea, byrow = TRUE)
    seasonal_mean <- as.numeric(apply(seasonal_matrix, 2, mean, na.rm = TRUE))
    seasonal_sd <- as.numeric(apply(seasonal_matrix, 2, sd, na.rm = TRUE))
    xt <- xt - rep(seasonal_mean, 1 + length(xt) %/% sea)[seq_along(xt)]
    if (season_sd) {
      xt <- xt / rep(seasonal_sd, 1 + length(xt) %/% sea)[seq_along(xt)]
    }
  }

  if (normalization) {
    xt <- (xt - mean(xt, na.rm = TRUE)) / sd(xt, na.rm = TRUE)
  }

  xt
}

analysis_data <- hko_ap_ha_0516 %>%
  select(
    date, DOW, Holiday, Tmean, Humidity, co, fsp, no, no2, nox, o3,
    o3h8max, rsp, so2, isc, hem
  ) %>%
  filter(year(date) >= 2008, year(date) <= 2016)

names(analysis_data) <- c(
  "date", "dow", "holiday", "temp", "rh", "co", "fsp", "no", "no2",
  "nox", "o3", "o3h8max", "rsp", "so2", "isc", "hem"
)

analysis_data$isc <- log(analysis_data$isc)
analysis_data <- as.data.frame(analysis_data)
for (nm in names(analysis_data)[4:ncol(analysis_data)]) {
  analysis_data[[nm]] <- normalize_series(analysis_data[[nm]])
}

write.csv(analysis_data, out_file, row.names = FALSE)
message("Wrote ", out_file)
