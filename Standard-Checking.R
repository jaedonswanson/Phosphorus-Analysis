#### Packages used in this code ####
library(tidyverse)
library(ggpmisc)

#### User inputs ####
# File path
data <- read.csv("06182026_Standards_only.csv")


# Potential problem standards
ROWS_TO_DROP <- c()   # only change this after checking the curve with all points

#### Concentration Vector creation ####
# Concentrations for the standards (ppm)
concentrations <- c(0, 3.1, 6.1, 11.9, 28.2, 51.7, 88.6)
data <- data %>% rename(Absorbance = X880nm..Abs.)
#### Building the data frame of just standards to be graphed ####
# New CSV format is already long: one row per measurement.
# Columns: #, Sample ID, User Name, Date and Time, 880nm (Abs)
# Each standard (std1–std7) appears 3 times — one row per replicate.

# Rename the absorbance column to something R-friendly
names(data)[names(data) == "880nm (Abs)"] <- "Absorbance"
names(data)[names(data) == "Sample ID"]   <- "Sample.ID"

standards_clean <- data %>%
  group_by(Sample.ID) %>%
  mutate(Replication = paste0("Rep_", row_number())) %>%  # label reps within each standard
  ungroup() %>%
  mutate(
    Standard_Level = as.integer(factor(Sample.ID, levels = unique(Sample.ID))), # 1–7 in order of appearance
    Concentration  = concentrations[Standard_Level],
    Flagged        = Standard_Level %in% ROWS_TO_DROP
  ) %>%
  select(Standard_Level, Concentration, Replication, Absorbance, Flagged)

#### Plot 1: All standards (Showing flagged if applicable) ####
standards_good    <- standards_clean %>% filter(!Flagged)
standards_flagged <- standards_clean %>% filter(Flagged)

fit_all       <- lm(Absorbance ~ Concentration, data = standards_clean)
r2_all        <- summary(fit_all)$r.squared
slope_all     <- coef(fit_all)[["Concentration"]]
intercept_all <- coef(fit_all)[["(Intercept)"]]

label_all <- paste0(
  "R² = ", round(r2_all, 4), " (all standards) \n ",
  "Slope = ", round(slope_all, 4), "\n",
  "Intercept = ", round(intercept_all, 4)
)

plot_all <- ggplot() +
  geom_line(
    data = standards_clean,
    aes(x = Concentration, y = Absorbance, group = Standard_Level),
    color = "gray80", linetype = "dashed"
  ) +
  geom_point(
    data = standards_good,
    aes(x = Concentration, y = Absorbance, color = Replication),
    size = 3, alpha = 0.85
  ) +
  geom_point(
    data = standards_flagged,
    aes(x = Concentration, y = Absorbance),
    color = "firebrick", size = 4, shape = 4,
    stroke = 1.5
  ) +
  geom_label(
    data = standards_flagged %>%
      group_by(Standard_Level, Concentration) %>%
      summarise(Absorbance = max(Absorbance), .groups = "drop"),
    aes(x = Concentration, y = Absorbance, label = "FLAGGED"),
    color = "firebrick", fill = "white", size = 3,
    vjust = -0.6,
    label.size = 0.3
  ) +
  geom_smooth(
    data = standards_clean,
    aes(x = Concentration, y = Absorbance),
    method = "lm", se = TRUE,
    color = "black", fill = "steelblue", alpha = 0.15,
    linewidth = 0.9
  ) +
  annotate(
    "label",
    x     = min(standards_good$Concentration, na.rm = TRUE),
    y     = max(standards_clean$Absorbance,   na.rm = TRUE),
    label = label_all,
    hjust = 0, vjust = 1,
    size  = 3.8, fontface = "bold",
    fill  = "white", label.size = 0.3
  ) +
  labs(
    title    = "Full Calibration Curve (Flagged Standards Shown)",
    subtitle = paste(
      "Red ✕ marks = flagged rows:", paste(ROWS_TO_DROP, collapse = ", "),
      "| Regression fitted to non-flagged points only"
    ),
    x     = "Concentration (µM)",
    y     = "Absorbance",
    color = "Replicate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40", size = 10)
  )
print(plot_all)

#### Plot 2 — Calibration curve with flagged standards removed ####
standards_filtered <- standards_clean %>% filter(!Flagged)

fit_filtered       <- lm(Absorbance ~ Concentration, data = standards_filtered)
r2_filtered        <- summary(fit_filtered)$r.squared
slope_filtered     <- coef(fit_filtered)[["Concentration"]]
intercept_filtered <- coef(fit_filtered)[["(Intercept)"]]

r2_direction <- if (r2_filtered > r2_all) "IMPROVED" else if (r2_filtered < r2_all) "WORSENED" else "unchanged"

label_filtered <- paste0(
  "R² = ", round(r2_filtered, 4), "\n",
  "Slope = ", round(slope_filtered, 4), "\n",
  "Intercept = ", round(intercept_filtered, 4)
)

cat("=== R² Comparison ===\n")
cat("All standards (excl. flagged):", round(r2_all, 4), "\n")
cat("Filtered standards only:       ", round(r2_filtered, 4), "\n\n")

plot_filtered <- ggplot(data = standards_filtered, aes(x = Concentration, y = Absorbance)) +
  geom_line(
    aes(group = Standard_Level),
    color = "gray80", linetype = "dashed"
  ) +
  geom_point(
    aes(color = Replication),
    size = 3, alpha = 0.85
  ) +
  geom_smooth(
    method = "lm", se = TRUE,
    color = "black", fill = "steelblue", alpha = 0.15,
    linewidth = 0.9
  ) +
  annotate(
    "label",
    x     = min(standards_filtered$Concentration, na.rm = TRUE),
    y     = max(standards_filtered$Absorbance,    na.rm = TRUE),
    label = label_filtered,
    hjust = 0, vjust = 1,
    size  = 3.8, fontface = "bold",
    fill  = "white", label.size = 0.3
  ) +
  labs(
    title    = "Filtered Calibration Curve",
    subtitle = paste0(
      "Rows removed: ", paste(ROWS_TO_DROP, collapse = ", "),
      " | R² ", r2_direction, " from ", round(r2_all, 4),
      " → ", round(r2_filtered, 4)
    ),
    x     = "Concentration (µM)",
    y     = "Absorbance (AU)",
    color = "Replicate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40", size = 10)
  )
print(plot_filtered)