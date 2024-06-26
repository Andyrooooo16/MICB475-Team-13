# Meeting Notes 2024/02/15

## Data Preparation and Initial Analysis

- [ ] **Filter Dataset for Time Point 1**
  - Filter the dataset to include only time point 1.
  - Remove any rows with NA values in the `resection` column.
  - Assess the amount of available data in the `resection` column to ensure there are enough values for meaningful analysis.

- [ ] **Basic Statistics on Calprotectin Values**
  - Calculate basic statistics for calprotectin values, including:
    - Mean
    - Median
    - Standard Deviation (std)
    - Minimum (min) value
    - Maximum (max) value

## Literature Review

- [ ] **Review Calprotectin in the Context of IBD**
  - Investigate existing literature on calprotectin, focusing on its role and significance in Inflammatory Bowel Disease (IBD).

## Further Analysis and Threshold Determination

- [ ] **Determine Value Thresholds for Calprotectin Levels**
  - Based on the initial analysis and literature review, decide on the value thresholds that categorize different levels of inflammation.

## Data Filtering and Manifest File Preparation

- [ ] **Filter Manifest File Post-Analysis**
  - After completing the filtering steps and analyses, ensure to filter the manifest file accordingly to reflect the refined dataset.

## Correlation Analysis

- [ ] **Correlate Calprotectin with CD Resection**
  - Perform a Mann-Whitney U test or Wilcoxon rank sum test to correlate calprotectin levels with CD (Crohn's Disease) resection.
  - This analysis aims to provide directional insight for our proposal.
     
## Goals for next week 

- [ ] **Proposal --> aim to get QIIME processing done (complete denoising step)**
