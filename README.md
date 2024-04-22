# 202402-MRI-behavioral-analysis

This R script is designed to analyze behavioral data collected from subjects undergoing MRI scans. It compares the subjects' actions and statements to their brain activity. Below is an overview of the functionality and steps performed by the script:

## Input Data
- Scanner Files: The script reads scanner files containing data from subjects' MRI sessions.
- Behavioral Questions: Behavioral questions are loaded from a CSV file.
- Meants Files: Meants files containing data related to brain activity are used.

## Data Processing
- The script preprocesses the scanner files and behavioral questions to prepare them for analysis.
- Dissimilarity matrices are calculated for various factors including targets, valence, words, and response time.

## Linear Modeling
- Linear models are fitted to predict neural dissimilarity from various behavioral dissimilarities.
- Separate linear models are run for different predictors such as targets, valence, words, and response time.

## Mantel Tests
- Mantel tests are conducted to assess the correlation between neural dissimilarity and behavioral dissimilarities.
- Mantel tests are run for each predictor separately.

## Export
- The results of linear models and Mantel tests are exported as CSV files for further analysis.

## Libraries Used
- dplyr: Data manipulation and transformation.
- rio: Importing and exporting data files.
- ade4: Dissimilarity matrix calculation.
- broom: Tidying model outputs.
- reshape2: Reshaping data for analysis.
