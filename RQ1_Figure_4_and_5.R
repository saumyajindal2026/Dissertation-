library(tidyverse)
library(readxl)
library(patchwork)

# Step 1: Define the source workbook and the output folder.
workbook_path <- "ICP MS DATA - working copy.xlsx"
sheet_name <- "metals_left_censored"
output_dir <- file.path("analysis_outputs", "icp_ms_point_ranges")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

treatment_order <- c("LWO2", "LWN2", "DIO2", "DIN2")
day_levels <- c("DAY 0", "DAY 1", "DAY 2")
filter_levels <- c("Filtered", "Unfiltered")

filter_colors <- c("Filtered" = "#1B9E77", "Unfiltered" = "#D95F02")
filter_shapes <- c("Filtered" = 16, "Unfiltered" = 17)
filter_linetypes <- c("Filtered" = "solid", "Unfiltered" = "dashed")
mean_shape <- 18

metal_specs <- list(
  list(
    metal = "Mn",
    source_col = "55 Mn Conc. [mg/l]",
    y_label = expression("Mn concentration (mg L"^{-1}*")"),
    file_stub = "Mn_point_range_4panel"
  ),
  list(
    metal = "Fe",
    source_col = "56 Fe Conc. [mg/l]",
    y_label = expression("Fe concentration (mg L"^{-1}*")"),
    file_stub = "Fe_point_range_4panel"
  )
)

save_plot_dual <- function(plot_obj, file_stub, width = 11, height = 8.5) {
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

# Step 2: Read the ICP-MS table and standardize the labels used in every figure.
icp_data <- read_excel(workbook_path, sheet = sheet_name) %>%
  filter(Replicate %in% c("3", "4", "5", "6")) %>%
  mutate(
    Treatment = substr(Sample_ID, 1, 4),
    Day = factor(trimws(Day), levels = day_levels),
    Filter = case_when(
      Filter == "F" ~ "Filtered",
      Filter == "UF" ~ "Unfiltered",
      TRUE ~ as.character(Filter)
    ),
    Treatment = factor(Treatment, levels = treatment_order),
    Filter = factor(Filter, levels = filter_levels)
  ) %>%
  filter(!is.na(Treatment), !is.na(Day), !is.na(Filter))

# Step 3: Build a reusable plotting function so the same code can produce Mn and Fe.
make_metal_point_range_plot <- function(df, metal_col, y_label, file_stub) {
  plot_data <- df %>%
    transmute(
      Treatment,
      Day,
      Filter,
      Concentration = as.numeric(.data[[metal_col]])
    ) %>%
    filter(!is.na(Concentration))

  summary_data <- plot_data %>%
    group_by(Treatment, Day, Filter) %>%
    summarise(
      n = n(),
      mean_concentration = mean(Concentration, na.rm = TRUE),
      se_concentration = sd(Concentration, na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    ) %>%
    mutate(se_concentration = if_else(is.finite(se_concentration), se_concentration, NA_real_))

  y_limits <- range(
    c(
      plot_data$Concentration,
      summary_data$mean_concentration - summary_data$se_concentration,
      summary_data$mean_concentration + summary_data$se_concentration
    ),
    na.rm = TRUE
  )

  dodge_width <- 0.35

  treatment_plots <- lapply(treatment_order, function(treatment_name) {
    current_data <- plot_data %>%
      filter(Treatment == treatment_name)

    ggplot(
      current_data,
      aes(
        x = Day,
        y = Concentration,
        color = Filter,
        shape = Filter,
        group = Filter
      )
    ) +
      geom_point(
        position = position_jitterdodge(
          jitter.width = 0.12,
          jitter.height = 0,
          dodge.width = dodge_width
        ),
        alpha = 0.6,
        size = 2
      ) +
      stat_summary(
        fun = mean,
        geom = "line",
        aes(linetype = Filter),
        position = position_dodge(width = dodge_width),
        linewidth = 0.8
      ) +
      stat_summary(
        fun = mean,
        geom = "point",
        position = position_dodge(width = dodge_width),
        shape = mean_shape,
        size = 3.4,
        stroke = 0.8,
        fill = "white"
      ) +
      stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        position = position_dodge(width = dodge_width),
        width = 0.12,
        linewidth = 0.55
      ) +
      scale_color_manual(values = filter_colors, drop = FALSE) +
      scale_shape_manual(values = filter_shapes, drop = FALSE) +
      scale_linetype_manual(values = filter_linetypes, drop = FALSE) +
      scale_x_discrete(drop = FALSE) +
      scale_y_continuous(limits = y_limits) +
      labs(
        title = treatment_name,
        x = "Day",
        y = y_label,
        color = "Fraction",
        shape = "Fraction",
        linetype = "Fraction"
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        legend.position = "top",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.margin = margin(8, 8, 8, 8)
      )
  })

  combined_plot <- wrap_plots(treatment_plots, ncol = 2, byrow = TRUE) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "top",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )

  save_plot_dual(combined_plot, file_stub)
}

# Step 4: Create and save the Mn and Fe figures.
for (spec in metal_specs) {
  make_metal_point_range_plot(
    df = icp_data,
    metal_col = spec$source_col,
    y_label = spec$y_label,
    file_stub = spec$file_stub
  )
}

cat("Saved ICP-MS 4-panel point-range figures to", output_dir, "\n")
