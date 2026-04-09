library(tidyverse)
library(readxl)
library(car)
library(patchwork)

workbook_path <- "/Users/saumyajindal17/DISSERTAION/blanks_lod_analysis_lod2_filled.xlsx"
sheet_name <- "metal_scatter_data"
output_dir <- file.path("analysis_outputs", "hypo1_gamma_glm")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_plot_dual <- function(plot_obj, file_stub, width = 12, height = 9) {
  ggsave(
    filename = file.path(output_dir, paste0(file_stub, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )

  ggsave(
    filename = file.path(output_dir, paste0(file_stub, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height
  )
}

build_glm_diagnostic_panel <- function(fit, model_label) {
  diag_data <- tibble(
    fitted = fitted(fit),
    deviance_resid = residuals(fit, type = "deviance"),
    pearson_resid = residuals(fit, type = "pearson"),
    std_resid = rstandard(fit),
    leverage = hatvalues(fit),
    cooks_d = cooks.distance(fit)
  ) %>%
    mutate(
      sqrt_abs_std_resid = sqrt(abs(std_resid))
    )

  base_theme <- theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      plot.margin = margin(8, 8, 8, 8)
    )

  resid_plot <- ggplot(diag_data, aes(x = fitted, y = deviance_resid)) +
    geom_point(size = 2, alpha = 0.8, color = "#2F2F2F") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7, color = "#7A7A7A") +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.9, color = "#A23E48") +
    labs(
      title = "Residuals vs Fitted",
      x = "Fitted values",
      y = "Deviance residuals"
    ) +
    base_theme

  qq_plot <- ggplot(diag_data, aes(sample = std_resid)) +
    stat_qq(size = 1.8, alpha = 0.8, color = "#2F2F2F") +
    stat_qq_line(linewidth = 0.9, color = "#A23E48") +
    labs(
      title = "Normal Q-Q",
      x = "Theoretical quantiles",
      y = "Standardized residuals"
    ) +
    base_theme

  scale_location_plot <- ggplot(diag_data, aes(x = fitted, y = sqrt_abs_std_resid)) +
    geom_point(size = 2, alpha = 0.8, color = "#2F2F2F") +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.9, color = "#A23E48") +
    labs(
      title = "Scale-Location",
      x = "Fitted values",
      y = expression(sqrt("|Standardized residuals|"))
    ) +
    base_theme

  leverage_max <- max(diag_data$leverage, na.rm = TRUE) * 1.1
  leverage_grid <- seq(0.001, max(leverage_max, 0.01), length.out = 200)
  num_params <- length(coef(fit))
  cooks_levels <- c(0.5, 1)

  cooks_df <- purrr::map_dfr(cooks_levels, function(level) {
    tibble(
      leverage = leverage_grid,
      std_resid = sqrt(level * num_params * (1 - leverage_grid) / leverage_grid),
      cook_label = paste0("Cook's D = ", level)
    )
  })

  leverage_plot <- ggplot(diag_data, aes(x = leverage, y = std_resid)) +
    geom_point(size = 2, alpha = 0.8, color = "#2F2F2F") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7, color = "#7A7A7A") +
    geom_line(
      data = cooks_df,
      aes(x = leverage, y = std_resid, linetype = cook_label),
      inherit.aes = FALSE,
      linewidth = 0.8,
      color = "#A23E48"
    ) +
    geom_line(
      data = cooks_df,
      aes(x = leverage, y = -std_resid, linetype = cook_label),
      inherit.aes = FALSE,
      linewidth = 0.8,
      color = "#A23E48"
    ) +
    scale_linetype_manual(values = c("Cook's D = 0.5" = "dashed", "Cook's D = 1" = "solid")) +
    labs(
      title = "Residuals vs Leverage",
      x = "Leverage",
      y = "Standardized residuals",
      linetype = NULL
    ) +
    base_theme +
    theme(legend.position = "bottom")

  (resid_plot + qq_plot) / (scale_location_plot + leverage_plot) +
    plot_annotation(
      title = paste(model_label, "diagnostics"),
      theme = theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
      )
    )
}

plot_data <- read_excel(workbook_path, sheet = sheet_name) %>%
  transmute(
    WaterType = case_when(
      str_to_lower(as.character(Water_Type)) %in% c("lake", "lake water") ~ "Lake",
      str_to_lower(as.character(Water_Type)) %in% c("di", "deionised", "deionized", "di water") ~ "Deionised",
      TRUE ~ as.character(Water_Type)
    ),
    Gas = case_when(
      str_to_lower(as.character(Gas_Condition)) %in% c("air", "oxygen", "o2", "oxic") ~ "Oxygen",
      str_to_lower(as.character(Gas_Condition)) %in% c("n2", "nitrogen", "anoxic") ~ "Nitrogen",
      TRUE ~ as.character(Gas_Condition)
    ),
    Day = factor(str_replace_all(as.character(Day), "DAY\\s*", "Day "), levels = c("Day 0", "Day 1", "Day 2")),
    Filter = as.character(Filter),
    Mn = as.numeric(Mn_mg_l),
    Fe = as.numeric(Fe_mg_l)
  ) %>%
  filter(Filter == "F") %>%
  filter(!is.na(WaterType), !is.na(Gas), !is.na(Day), !is.na(Mn), !is.na(Fe)) %>%
  mutate(
    WaterType = factor(WaterType, levels = c("Lake", "Deionised")),
    Gas = factor(Gas, levels = c("Oxygen", "Nitrogen"))
  )

if (any(plot_data$Mn <= 0, na.rm = TRUE) || any(plot_data$Fe <= 0, na.rm = TRUE)) {
  stop("Gamma GLM requires strictly positive Mn and Fe values.")
}

# Gamma(log) is usually more appropriate than ANOVA for strictly positive,
# right-skewed concentration data because it models the mean-variance pattern
# directly rather than trying to force normally distributed residuals.
mn_fit <- glm(Mn ~ WaterType * Day * Gas, data = plot_data, family = Gamma(link = "log"))
fe_fit <- glm(Fe ~ WaterType * Day * Gas, data = plot_data, family = Gamma(link = "log"))

mn_diagnostic_panel <- build_glm_diagnostic_panel(mn_fit, "Gamma(log) Mn")
fe_diagnostic_panel <- build_glm_diagnostic_panel(fe_fit, "Gamma(log) Fe")

mn_anova <- car::Anova(mn_fit, type = 2, test.statistic = "LR") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  mutate(model = "Gamma(log) Mn ~ WaterType * Day * Gas", .before = 1)

fe_anova <- car::Anova(fe_fit, type = 2, test.statistic = "LR") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term") %>%
  mutate(model = "Gamma(log) Fe ~ WaterType * Day * Gas", .before = 1)

coeff_tbl <- bind_rows(
  broom::tidy(mn_fit) %>% mutate(model = "Gamma(log) Mn", .before = 1),
  broom::tidy(fe_fit) %>% mutate(model = "Gamma(log) Fe", .before = 1)
)

fit_tbl <- tibble(
  model = c("Gamma(log) Mn", "Gamma(log) Fe"),
  aic = c(AIC(mn_fit), AIC(fe_fit)),
  residual_deviance = c(deviance(mn_fit), deviance(fe_fit)),
  df_residual = c(df.residual(mn_fit), df.residual(fe_fit)),
  dispersion = c(summary(mn_fit)$dispersion, summary(fe_fit)$dispersion)
)

summary_lines <- c(
  "Hypothesis 1 Gamma GLM check",
  "",
  "Why this model:",
  "Gamma(log) is suitable for positive, right-skewed concentration data and does not require normal residuals in the same way as ANOVA.",
  "",
  "Model fit summary",
  capture.output(print(fit_tbl)),
  "",
  "Mn likelihood-ratio ANOVA",
  capture.output(print(mn_anova)),
  "",
  "Fe likelihood-ratio ANOVA",
  capture.output(print(fe_anova))
)

write_csv(bind_rows(mn_anova, fe_anova), file.path(output_dir, "H1_gamma_glm_anova_tables.csv"), na = "")
write_csv(coeff_tbl, file.path(output_dir, "H1_gamma_glm_coefficients.csv"), na = "")
write_csv(fit_tbl, file.path(output_dir, "H1_gamma_glm_fit_summary.csv"), na = "")
write_lines(summary_lines, file.path(output_dir, "H1_gamma_glm_summary.txt"))
save_plot_dual(mn_diagnostic_panel, "H1_gamma_glm_Mn_diagnostic_panel")
save_plot_dual(fe_diagnostic_panel, "H1_gamma_glm_Fe_diagnostic_panel")

cat("Gamma GLM fit summary:\n")
print(fit_tbl)
cat("\nMn Gamma GLM likelihood-ratio ANOVA:\n")
print(mn_anova)
cat("\nFe Gamma GLM likelihood-ratio ANOVA:\n")
print(fe_anova)
