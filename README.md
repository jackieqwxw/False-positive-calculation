# False Positive Calculation

**False-positive-calculation** is a tool designed for calculating false positives based on a confusion matrix table. It is particularly useful for determining thresholds to distinguish between active and inactive compounds using multiple free energy methods, including MM/GBSA, MM/PBSA, and free energy perturbation.

## Installation
Install Jupyter Notebook using Anaconda in a Windows environment.

## Explanation
1. Input files were saved as `.csv`.
   
2. Specify the target names using `target_titles = ['ACE', 'ADRB1', 'FAK1', 'GRIK1', 'HMDH', 'MCR', 'PGH2', 'PRGR', 'TRYB1']`. Adjust the target names to match your specific targets.

3. The tool initially calculates false positives for all targets based on a range of threshold values.

4. It then selects 7 out of 9 targets through combinations to determine the minimum thresholds. Simultaneously, the false positive values for each group are saved as a text file with a prefix name `false_positive_`.
