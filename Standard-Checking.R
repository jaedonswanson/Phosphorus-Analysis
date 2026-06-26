#### Packages ####
library(tidyverse)

#### User inputs ####
# Path to your CSV file (change as needed)
data <- read.csv("Data/06262026.csv")

# Which standard levels (1–7) should be flagged as problematic?
# Example: ROWS_TO_DROP <- c(3, 5)  # drops std3 and std5
ROWS_TO_DROP <- c()   # modify after inspecting the full curve

#### Concentration vector (µM) – index matches standard level ####
# std1 → 0, std2 → 3.1, …, std7 → 88.6
concentrations <- c(0, 3.1, 6.1, 11.9, 28.2, 51.7, 88.6)
JITTER_WIDTH <- 0.15
#### Clean and rename columns ####
data <- data %>%
  rename(
    Absorbance = `X880nm..Abs.`,
    Sample.ID  = `Sample.ID`
  )

#### Extract only standards and parse level/replicate ####
standards_raw <- data %>%
  filter(str_detect(Sample.ID, "^std")) %>%          # keep only standards
  mutate(
    # Extract the numeric level (e.g., "std1_1" → 1)
    Standard_Level = as.numeric(str_extract(Sample.ID, "(?<=std)\\d+")),
    # Extract replicate number (e.g., "std1_1" → 1)
    Replication    = as.numeric(str_extract(Sample.ID, "\\d+$"))
  ) %>%
  filter(!is.na(Standard_Level)) %>%                 # safety check
  mutate(
    Concentration = concentrations[Standard_Level],  # map to concentration
    Flagged       = Standard_Level %in% ROWS_TO_DROP
  ) %>%
  select(Standard_Level, Concentration, Replication, Absorbance, Flagged)

#### Plot 1: All standards (flagged shown) ####
standards_good    <- standards_raw %>% filter(!Flagged)
standards_flagged <- standards_raw %>% filter(Flagged)

# Fit on all data (but flagged points are still in the regression? 
# The code uses standards_clean for the fit, so it includes flagged points.
# We'll fit using all points, but the plot will show flagged separately.
fit_all       <- lm(Absorbance ~ Concentration, data = standards_raw)
r2_all        <- summary(fit_all)$r.squared
slope_all     <- coef(fit_all)[["Concentration"]]
intercept_all <- coef(fit_all)[["(Intercept)"]]

label_all <- paste0(
  "R² = ", round(r2_all, 4), " (all standards)\n",
  "Slope = ", round(slope_all, 4), "\n",
  "Intercept = ", round(intercept_all, 4)
)

plot_all <- ggplot() +
  # Dashed lines connecting replicates of the same standard level
  geom_line(
    data = standards_raw,
    aes(x = Concentration, y = Absorbance, group = Standard_Level),
    color = "gray80", linetype = "dashed"
  ) +
  # Good points (colored by replicate)
  geom_point(
    data = standards_good,
    aes(x = Concentration, y = Absorbance, color = as.factor(Replication)),
    size = 3, alpha = 0.85,
    position = position_jitter(width = JITTER_WIDTH, height = 0)
  ) +
  # Flagged points (red X)
  geom_point(
    data = standards_flagged,
    aes(x = Concentration, y = Absorbance),
    color = "firebrick", size = 4, shape = 4, stroke = 1.5,
    position = position_jitter(width = JITTER_WIDTH, height = 0)
  ) +
  # Label for flagged points
  geom_label(
    data = standards_flagged %>%
      group_by(Standard_Level, Concentration) %>%
      summarise(Absorbance = max(Absorbance), .groups = "drop"),
    aes(x = Concentration, y = Absorbance, label = "FLAGGED"),
    color = "firebrick", fill = "white", size = 3,
    vjust = -0.6, label.size = 0.3
  ) +
  # Regression line (fitted on all points)
  geom_smooth(
    data = standards_raw,
    aes(x = Concentration, y = Absorbance),
    method = "lm", se = TRUE,
    color = "black", fill = "steelblue", alpha = 0.15,
    linewidth = 0.9
  ) +
  annotate(
    "label",
    x = min(standards_raw$Concentration, na.rm = TRUE),
    y = max(standards_raw$Absorbance,   na.rm = TRUE),
    label = label_all,
    hjust = 0, vjust = 1,
    size = 3.8, fontface = "bold",
    fill = "white", label.size = 0.3
  ) +
  labs(
    title    = "Full Calibration Curve (Flagged Standards Shown)",
    subtitle = paste(
      "Red ✕ marks = flagged levels:",
      paste(ROWS_TO_DROP, collapse = ", "),
      "| Regression uses all points"
    ),
    x = "Concentration (µM)",
    y = "Absorbance",
    color = "Replicate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40", size = 10)
  )
print(plot_all)

#### Plot 2: Filtered curve (flagged levels removed) ####
standards_filtered <- standards_raw %>% filter(!Flagged)

fit_filtered <- lm(Absorbance ~ Concentration, data = standards_filtered)
r2_filtered  <- summary(fit_filtered)$r.squared
slope_filtered     <- coef(fit_filtered)[["Concentration"]]
intercept_filtered <- coef(fit_filtered)[["(Intercept)"]]

r2_direction <- if (r2_filtered > r2_all) "IMPROVED" else if (r2_filtered < r2_all) "WORSENED" else "unchanged"

label_filtered <- paste0(
  "R² = ", round(r2_filtered, 4), "\n",
  "Slope = ", round(slope_filtered, 4), "\n",
  "Intercept = ", round(intercept_filtered, 4)
)

cat("=== R² Comparison ===\n")
cat("All standards (incl. flagged):", round(r2_all, 4), "\n")
cat("Filtered standards only:      ", round(r2_filtered, 4), "\n\n")

plot_filtered <- ggplot(standards_filtered, aes(x = Concentration, y = Absorbance)) +
  geom_line(
    aes(group = Standard_Level),
    color = "gray80", linetype = "dashed"
  ) +
  geom_point(
    aes(color = as.factor(Replication)),
    size = 3, alpha = 0.85,
    position = position_jitter(width = JITTER_WIDTH, height = 0)
  ) +
  geom_smooth(
    method = "lm", se = TRUE,
    color = "black", fill = "steelblue", alpha = 0.15,
    linewidth = 0.9
  ) +
  annotate(
    "label",
    x = min(standards_filtered$Concentration, na.rm = TRUE),
    y = max(standards_filtered$Absorbance,    na.rm = TRUE),
    label = label_filtered,
    hjust = 0, vjust = 1,
    size = 3.8, fontface = "bold",
    fill = "white", label.size = 0.3
  ) +
  labs(
    title    = "Filtered Calibration Curve",
    subtitle = paste0(
      "Removed levels: ", paste(ROWS_TO_DROP, collapse = ", "),
      " | R² ", r2_direction, " from ", round(r2_all, 4),
      " → ", round(r2_filtered, 4)
    ),
    x = "Concentration (µM)",
    y = "Absorbance (AU)",
    color = "Replicate"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "gray40", size = 10)
  )
print(plot_filtered)