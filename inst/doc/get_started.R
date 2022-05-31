## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
#################### required #################### 
library(portvine)
#################### optional #################### 
## the following packages are only necessary for the visualizations
# for visualization:
library(ggplot2)
library(patchwork)
theme_set(
  theme_minimal() +
  theme(plot.title = ggtext::element_markdown(size = 11),
        plot.subtitle = ggtext::element_markdown(size = 9))
)
# some data wrangling for visualizations
library(dplyr)
library(tidyr)

## -----------------------------------------------------------------------------
set.seed(2)

## -----------------------------------------------------------------------------
data("sample_returns_small")
head(sample_returns_small)
dim(sample_returns_small)

## -----------------------------------------------------------------------------
default_garch_spec()

## -----------------------------------------------------------------------------
uncond_marg_settings <- marginal_settings(
  train_size = 960,
  refit_size = 20
)

## -----------------------------------------------------------------------------
cond_marg_settings <- marginal_settings(
  train_size = 960,
  refit_size = 20,
  individual_spec = list("AMZN" = default_garch_spec()),
  default_spec = default_garch_spec(ma = 2)
)

## -----------------------------------------------------------------------------
uncond_vine_settings <- vine_settings(
  train_size = 100,
  refit_size = 10,
  family_set = "onepar", # valid bivariate building blocks
  vine_type = "rvine"
)

## -----------------------------------------------------------------------------
cond_vine_settings <- vine_settings(
  train_size = 100,
  refit_size = 20,
  family_set = c("gumbel", "joe", "t"),
  vine_type = "dvine"
)

## -----------------------------------------------------------------------------
## specify parallel strategy
# future::plan("multisession", workers = 4)
uncond_risk_roll <- estimate_risk_roll(
  data = sample_returns_small,
  weights = NULL,
  marginal_settings = uncond_marg_settings,
  vine_settings = uncond_vine_settings,
  alpha = c(0.01, 0.025),
  risk_measures = c("VaR", "ES_mean"),
  n_samples = 50,
  trace = TRUE
)
## return to default sequential settings (cut off any background processes)
# future::plan("sequential")

## -----------------------------------------------------------------------------
uncond_risk_roll

## -----------------------------------------------------------------------------
summary(uncond_risk_roll)

## -----------------------------------------------------------------------------
cond_risk_roll <- estimate_risk_roll(
  data = sample_returns_small,
  weights = c("AAPL" = 1, "GOOG" = 2, "AMZN" = 0),
  marginal_settings = cond_marg_settings,
  vine_settings = cond_vine_settings,
  alpha = c(0.01, 0.025),
  risk_measures = c("VaR", "ES_mean"),
  n_samples = 50,
  cond_vars = "AMZN",
  cond_u = 0.1
)

## -----------------------------------------------------------------------------
cond_risk_roll

## -----------------------------------------------------------------------------
summary(cond_risk_roll)

## -----------------------------------------------------------------------------
names(fitted_marginals(uncond_risk_roll))

## -----------------------------------------------------------------------------
# An opinionated function for residual analysis. Note this function has no 
# input checks and is not tested in any ways.
# Output: named list with a composition of residual plots for each asset
marg_resid_viz_list <- function(
  roll, # portvine_roll or cond_portvine_roll object
  asset_names = NULL, # filter for certain  assets
  marg_window = 1, # specify a marginal window
  squared = FALSE # if set to true the results for the squared stand. residuals
  # are displayed
  ) {
  fitted_marginals <- fitted_marginals(roll)
  if (is.null(asset_names)) asset_names <- fitted_vines(roll)[[1]]$names
  sapply(
    asset_names,
    function(asset_name) {
      # use again a utility function from the portvine package to extract the 
      # fitted standardized residuals
      model_resid <- roll_residuals(
        fitted_marginals[[asset_name]], marg_window
      )
      if (squared) model_resid <- model_resid^2

      simple_exploratory <- data.frame(resid = model_resid) %>%
        mutate(id = seq(length(model_resid))) %>%
        ggplot(aes(x = id, y = resid)) +
        geom_line(size = 0.2) +
        labs(x = "t", y = ifelse(squared, expression(z[t]^2),expression(z[t])),
             title = asset_name) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

      acf_plot <- data.frame(
        acf = as.numeric(acf(model_resid, type = "cor", lag.max = 20,
                             plot = FALSE)$acf),
        lag = 0:20
      ) %>%
        filter(lag != 0 & lag <= 10) %>%
        ggplot() +
        geom_hline(yintercept = 0, col = "black", size = 0.3) +
        geom_hline(yintercept = qnorm(c(0.025, 0.975)) /
                     sqrt(length(model_resid)),
                   linetype = "longdash", col = "#92B8DE", size = 0.5) +
        geom_segment(aes(x = lag, xend = lag, y = 0, yend = acf)) +
        geom_point(aes(x = lag, y = acf)) +
        scale_x_continuous(breaks = seq(1, 10, 1)) +
        ylim(-1, 1) +
        labs(x = "h", y = "ACF(h)")

      ljungplot <- data.frame(
        pval = sapply(
          1:10,
          function(i) Box.test(model_resid, lag = i,
                               type = "Lju")$p.value),
        lag = 1:10) %>%
        ggplot() +
        geom_hline(yintercept = 0, col = "black", size = 0.3) +
        geom_hline(yintercept = 0.05,
                   linetype = "longdash", col = "#92B8DE", size = 0.5) +
        geom_line(aes(x = lag, y = pval)) +
        geom_point(aes(x = lag, y = pval)) +
        scale_x_continuous(breaks = seq(1, 10, 1)) +
        labs(x = "h", y = "p-value of Ljung-Box test at lag h")

      (simple_exploratory / (ljungplot + acf_plot)) +
        plot_layout(nrow = 2)
    }, USE.NAMES = TRUE, simplify = FALSE)
}

## ---- results='hide'----------------------------------------------------------
uncond_residual_viz <- marg_resid_viz_list(uncond_risk_roll)
uncond_residual_viz_squared <- marg_resid_viz_list(
  uncond_risk_roll, squared = TRUE
)

## -----------------------------------------------------------------------------
uncond_residual_viz$AAPL
uncond_residual_viz_squared$GOOG

## -----------------------------------------------------------------------------
# The function creates a heatmap of Ljung Box test p-values for a 
# (cond-)portvine_roll object. Note this function again has no 
# input checks and is not tested in any ways.
# Output: The heatmap (ggplot2 object)
ljung_heatmap <- function(roll, roll_num = 1, squared = FALSE) {
  asset_names <- fitted_vines(roll)[[1]]$names
  roll_marginals <- fitted_marginals(roll)
  ljung_data <- sapply(asset_names, function(asset_name) {
    model_resid <- roll_residuals(
      roll_marginals[[asset_name]], roll_num = roll_num
    )
    if (squared) model_resid <- model_resid^2
    sapply(
      1:10,
      function(i) Box.test(model_resid, lag = i, type = "Lju")$p.value
    )
  }, USE.NAMES = TRUE, simplify = TRUE)
  ljung_data <- ljung_data %>%
    as.data.frame() %>%
    mutate(lag = seq(nrow(ljung_data))) %>%
    pivot_longer(-lag, names_to = "asset", values_to = "pval")
  if (all(ljung_data$pval >= 0.05)) {
    legend_scale <- scale_fill_gradient(
      low = "#92B8DE", high = "#2a82db"
    )
  } else {
    legend_scale <- scale_fill_gradientn(
      colours = c("#db4f59","#C37285" ,
                  "#92B8DE", "#2a82db"),
      values = scales::rescale(c(0, 0.05 - 0.01, 0.05, 1)),
      breaks = c(0.05),
      labels = c(0.05),
      guide = guide_colourbar(nbin = 1000))
  }
  ljung_data %>%
    ggplot(aes(x = lag, y = asset, fill = pval)) +
    geom_tile() +
    scale_x_continuous(breaks = 1:10) +
    labs(y = "", x = "h", fill = "p-value",
         title = "Results of the Ljung-Box tests",
         caption = paste("Rolling window:", roll_num)) +
    legend_scale +
    theme(legend.position = "right",
          panel.grid.minor.x = element_blank())
}

## -----------------------------------------------------------------------------
ljung_heatmap(uncond_risk_roll)

## -----------------------------------------------------------------------------
uncond_fitted_vines <- fitted_vines(uncond_risk_roll)
uncond_fitted_vines[[3]]

## -----------------------------------------------------------------------------
## a glimpse at the tidy data.frame with all estimated risk measures
head(risk_estimates(uncond_risk_roll), 10)

## a glimpse at the estimated VaR at confidence level 1% with exceeded column
head(
  risk_estimates(
    uncond_risk_roll,
    risk_measures = "VaR",
    alpha = 0.01,
    exceeded = TRUE
  ),
  6
)

## of course also applies for the conditional case
head(
  risk_estimates(
    cond_risk_roll,
    risk_measures = "ES_mean",
    alpha = 0.01,
    cond_u = c("resid", 0.1)
  ),
  6
)

## -----------------------------------------------------------------------------
risk_estimates(uncond_risk_roll, alpha = 0.025) %>%
  ggplot() +
  geom_line(aes(x = row_num, y = realized), col = "grey") +
  geom_line(aes(x = row_num, y = risk_est,
                col = factor(risk_measure))) +
  scale_color_manual(values = c("#92B8DE", "#477042")) +
  labs(x = "trading day",
       y = "portfolio log returns",
       col = "Risk measure",
       title = "Comparison of unconditional risk measure estimates at 
       confidence level 2.5%",
       subtitle = paste0("Realized log returns in ",
                     "<span style='color:",
                     "grey",
                     "'>**grey**</span>",
                     "."))

## -----------------------------------------------------------------------------
risk_estimates(
  uncond_risk_roll,
  risk_measures = c("VaR"),
  alpha = 0.01,
  exceeded = TRUE) %>%
  ggplot() +
  geom_line(aes(x = row_num, y = realized), col = "lightgrey") +
  geom_line(aes(x = row_num, y = risk_est), col = "#92B8DE") +
  geom_point(aes(x = row_num, y = realized), col = "#db4f59", 
             inherit.aes = FALSE, data = . %>% filter(exceeded)) +
  labs(x = "trading days",
       y = "portfolio log returns",
       col = "Exceeded",
       subtitle = paste0("Exceedances are highlighted in ",
                     "<span style='color:",
                     "#db4f59",
                     "'>**red**</span>",
                     "."),
       title = "Risk measure: VaR with confidence level 1%"
       ) +
  theme(legend.position = "none")

## -----------------------------------------------------------------------------
risk_estimates(
  cond_risk_roll,
  risk_measures = "ES_mean",
  alpha = 0.01) %>%
  ggplot() +
  geom_line(
    data = sample_returns_small[961:1000, "AMZN"],
    aes(x = 961:1000,
        y = AMZN), col = "grey", size = .3) +
  geom_line(aes(x = row_num, y = AMZN, col = factor(cond_u)),
            size = 0.5) +
  scale_color_manual(values = c("#92B8DE", "#477042"),
                     name = "Conditional strategy") +
  labs(title = "Conditional variable: AMZN",
       x = "trading days", y = "log returns of the asset AMZN",
       subtitle = paste0("Realized log returns in ",
                     "<span style='color:",
                     "grey",
                     "'>**grey**</span>",
                     "."))

## -----------------------------------------------------------------------------
risk_estimates(cond_risk_roll, alpha = 0.025, risk_measures = "ES_mean") %>%
  ggplot() +
  geom_line(aes(x = row_num, y = realized), col = "grey") +
  geom_line(aes(x = row_num, y = risk_est,
                col = factor(cond_u))) +
  scale_color_manual(values = c("#477042", "#92B8DE")) +
  labs(x = "trading day",
       y = "portfolio log returns",
       col = "Conditional strategy",
       title = "Comparison of conditional risk measure estimates at 
       confidence level 2.5%",
       subtitle = paste0("Realized log returns in ",
                     "<span style='color:",
                     "grey",
                     "'>**grey**</span>",
                     "."))

