# EVSI Trial Design Analysis

R functions for calculating Expected Value of Sample Information (EVSI), trial costs, and Expected Net Benefit of Sampling (ENBS) for the early health economic evaluation of a biomarker test.

## Overview

- **`evsi_diagnostic.R`** - EVSI and ENBS calculation using importance sampling
- **`trial_costs.R`** - Trial cost estimation with uncertainty quantification based on van Asselt et al. (2018)

## Dependencies

- `mgcv` - GAM fitting
- `ggplot2` - Visualization
- `dplyr` - Data manipulation
- `voi` - EVSI methods

## References

- Fernández Coves A., Ramaekers B., Joore M., Grimm S. What We Don’t Model Matters: Quantifying Bias in Research Prioritization for Biomarkers Based on Value of Information Analysis. (Manuscript submitted to Society of Medical Decision Making)
- van Asselt, T., B. Ramaekers, I. Corro Ramos, M. Joore, M. Al, I. Lesman-Leegte, M. Postma, P. Vemer, and T. Feenstra, Research Costs Investigated: A Study Into the Budgets of Dutch Publicly Funded Drug-Related Research. PharmacoEconomics, 2018. 36(1): p. 105-113. https://doi.org/10.1007/s40273-017-0572-7
