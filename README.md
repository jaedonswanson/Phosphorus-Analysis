# Analysis of Phosphate via Spec
This repository contains everything needed to prepate, analyze, and caluclate PO$_{4}^{3-}$ concentrations for freshwater samples

## A note on sample collection
If only collecting water samples, they should be filtered through a 0.45 micron filter. If the samples cannot be processed the same day, freeze them.

# Part 1: Setting up the experiment
1. Read though the [procudure](SRP-Procedure.md) 

2. Set up an Excel file to act as a key for the sample locations.
   1. If not formatted like [this key](phosphorus_key.xlsx), the code following this will break.
3. After the incubation time is complete, begin to run all samples at 880 nm on the spec 
4. Copy the data to a CSV and name it `MMDDYYYY.csv`
   1. This is essential to make sure the plate ID matches the corresponding key.

## Part 2: Assessing the Standards
1. Open the [standard checking](standard-checking.qmd) and update the relevant information.
2. Run the script and examine the $r^2$ value.
   1. Ideally, all standards produce an $r^2$ ≥ 0.98.
3. If $r^2$ is < 0.98, examine the graph and remove problem standards.
   1. **Notes:**
      1. If possible, do not remove a standard that best approximates sample concentration.
      2. If you have to remove more than 3 standards (only 5 points remaining), something likely went wrong and you should start over.

# Part 3: Calculating concentrations
1. Open the [concentration calculation](concentration-calculation.qmd) script.
2. Update the information to exactly match the settings from above (e.g., plate ID, removed standards).
3. Run the script. 
4. It should output a file with the concentrations of your samples.

# Part 4: Tidying
1. Whenever you want a clean file of your data, run the [tidying script](tidying.qmd).
   1. **Note:** Update site names and output paths to reflect your desired output.