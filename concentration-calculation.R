#### Setting Up the Environment ####
library(tidyverse)
library(readxl)

#### User Inputs (Edit these as needed) ####
# Path to your absorbance CSV file (e.g., "06262026.csv")
plate_csv <- "Data/06262026.csv"

# Path to the phosphorus key Excel file
key_excel <- "phosphorus_key.xlsx"

# Minimum acceptable R² for the calibration curve
min_r2 <- 0.98

# Stop execution if R² is below min_r2? (TRUE = stop, FALSE = continue with warning)
stop_if_low_r2 <- TRUE

# Standard levels (1–7) to drop from the calibration (e.g., c(3,5) drops std3 and std5)
bad_standards <- c()   # change after inspecting the curve

#### Concentration vector for standards (µM) ####
# std1 = 0, std2 = 3.1, ..., std7 = 88.6
conc_vector <- c(0, 3.1, 6.1, 11.9, 28.2, 51.7, 88.6)

#### Read the absorbance data ####
data <- read.csv(plate_csv, check.names = FALSE)
# Rename columns for clarity
colnames(data) <- c("Index", "Sample_ID", "User", "DateTime", "Absorbance")

#### Determine the correct sheet in the key Excel ####
# The key sheet name matches the date part of the CSV file name (e.g., "06262026")
sheet_name <- str_extract(basename(plate_csv), "^[0-9]+")
if (is.na(sheet_name)) {
  stop("Could not extract date from file name. Please ensure CSV name starts with a date (e.g., 06262026.csv).")
}

# Read the key (vial number → sample name)
key <- read_excel(key_excel, sheet = sheet_name)
colnames(key) <- c("Vial", "Sample", "Total")   # 'Total' is not used

#### Separate standards and unknowns from the absorbance data ####
# Standards: Sample_ID starts with "std"
standards_raw <- data %>%
  filter(str_detect(Sample_ID, "^std")) %>%
  mutate(
    Standard_Level = as.numeric(str_extract(Sample_ID, "(?<=std)\\d+")),
    Replication    = as.numeric(str_extract(Sample_ID, "\\d+$")),
    Concentration  = conc_vector[Standard_Level],
    Flagged        = Standard_Level %in% bad_standards
  ) %>%
  select(Standard_Level, Replication, Concentration, Absorbance, Flagged)

# Unknowns: Sample_ID is a numeric vial number
unknowns_raw <- data %>%
  filter(!str_detect(Sample_ID, "^std")) %>%
  mutate(Vial = as.numeric(Sample_ID)) %>%
  filter(!is.na(Vial))   # drop any non‑numeric entries

# Join unknowns with the key to get sample names
unknowns_with_key <- unknowns_raw %>%
  left_join(key, by = c("Vial" = "Vial")) %>%
  filter(!is.na(Sample))   # remove any vial not found in the key

# Assign a replication number within each sample (based on vial order)
unknowns_final <- unknowns_with_key %>%
  group_by(Sample) %>%
  mutate(Replication = row_number()) %>%
  ungroup() %>%
  select(Sample, Replication, Vial, Absorbance)

#### Fit the calibration curve (using non‑flagged standards) ####
standards_for_model <- standards_raw %>%
  filter(!Flagged) %>%
  drop_na(Absorbance)

fit_model <- lm(Absorbance ~ Concentration, data = standards_for_model)
slope     <- coef(fit_model)[["Concentration"]]
intercept <- coef(fit_model)[["(Intercept)"]]
current_r2 <- summary(fit_model)$r.squared

#### Quality Control and diagnostic output ####
cat("\n--- Calibration Performance ---\n")
cat("Concentration range: 0 – 88.6 µM\n")
cat("Omitted Standards:   ", if (length(bad_standards) == 0) "None" else paste(bad_standards, collapse = ", "), "\n")
cat("Achieved R²:         ", round(current_r2, 5), "\n")

if (current_r2 >= min_r2) {
  cat("✅ QC PASSED: R² meets the", min_r2, "threshold.\n\n")
} else {
  cat("\n=============================================================================\n")
  cat("❌ QC WARNING: R² (", round(current_r2, 4), ") is BELOW the", min_r2, "target!\n")
  cat("=============================================================================\n")
  cat("Running diagnostic: try dropping one standard level at a time...\n")
  
  for (i in 1:7) {
    sim_standards <- standards_raw %>%
      filter(Standard_Level != i) %>%
      drop_na(Absorbance)
    if (nrow(sim_standards) > 3) {   # need at least 3 points for a reasonable fit
      sim_fit <- lm(Absorbance ~ Concentration, data = sim_standards)
      sim_r2  <- summary(sim_fit)$r.squared
      status  <- if (sim_r2 >= min_r2) " [🎯 TARGET MET!]" else ""
      cat("  Omit level", i, ": R² =", round(sim_r2, 4), status, "\n")
    }
  }
  cat("\n💡 Recommendation: Add the problematic level(s) to 'bad_standards' and re‑run.\n")
  if (stop_if_low_r2) {
    stop("Execution halted due to low R² quality.")
  }
}

#### Calculate concentrations for unknowns ####
unknowns_with_conc <- unknowns_final %>%
  mutate(Concentration = (Absorbance - intercept) / slope) %>%
  select(Sample, Replication, Vial, Absorbance, Concentration)

#### Export the results ####
# Create output folder "Concentrations" if it doesn't exist
out_dir <- "Concentrations"
if (!dir.exists(out_dir)) dir.create(out_dir)

# Build output file name: same base as input CSV + "_concentrations.csv"
base_name <- str_remove(basename(plate_csv), "\\.csv$")
out_file <- file.path(out_dir, paste0(base_name, "_concentrations.csv"))

write.csv(unknowns_with_conc, out_file, row.names = FALSE)

cat("=========================================================================\n")
cat("✅ Concentration calculation complete.\n")
cat("Saved results to:", out_file, "\n")
cat("=========================================================================\n")