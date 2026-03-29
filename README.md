**This repository contains the updated files to replicate the results in "The rising cost of climate change: evidence from the bond market" by M. Bauer and G. Rudebusch (Review of Economics and Statistics, 2023).**

The codes in this folder are an improved version of the original replication package provided by the authors.

# DESCRIPTION OF THE DATASETS
- GS1.csv: one-year constant-maturity Treasury yield from FRED
- GS10.csv: ten-year constant-maturity Treasury yield from FRED
- HISTDATA_FRBUS.TXT: data package for FRB/US model of the Federal Reserve Board
    - See https://www.federalreserve.gov/econres/us-models-package.htm
- MedianGrowthRate_CPI.csv: median CPI inflation expectation from Livingston survey
    - See https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/livingston-historical-data
- r* estimates from literature, for references see Online Appendix 
    - rstar_delnegro.csv: Del Negro et al. (2017)
    - rstar_jm.csv: Johannsen and Mertense (2016)
    - rstar_kiley.csv: Kiley (2020)
    - rstar_lw (2016): Laubach and Williams (2016)
    - Note: These are all *smoothed estimates* based on the state space models in these published paper.
- Damage profiles for DICE models from the literature, for details see the paper and the Online Appendix
    - dice_2016_damages.csv: DICE-2016 model (Nordhaus, 2017)
    - dice_newell_pizer_damages.csv: DICE-94 model (Nordhaus, 1994) used in Newell and Pizer (2003)
    - dice_dfg_damages.csv: DICE-Fair-Geoffroy model of Dietz et al. (2020)
    - dice_haensel_damages.csv: DICE model of Haensel et al. (2020)

# SOFTWARE DESCRIPTION

The following codes have been created with R version: 4.5.1

R packages (also load dependencies)
- dplyr (1.2.0)
- dynlm (0.3-6)
- sandwich (3.0-1)
- KFAS (1.6.0)
- ggplot2 (4.0.2)
- sandwich (3.1-1)
- zoo (1.8-14)

How to set up your folder
- The R code assumes that the working directory is the ROOT folder, which must be organized as follows

ROOT (your main working directory)\
&nbsp;&nbsp;|- data (stores data)\
&nbsp;&nbsp;|- R (stores scripts)\
&nbsp;&nbsp;|- results (stores .RData files produced by scripts)\
&nbsp;&nbsp;&nbsp;&nbsp;       |- plots (stores plots produced by scripts)\
&nbsp;&nbsp;&nbsp;&nbsp;       |- tabs (stores tables produced by scripts)\

# Replication scripts
All the scripts can be called from 'main.R', in the following order:

0. Utility code files
These files contain functions that are used by the other scripts.
- **setup.R**: makes sure the relevant R packages are installed
- **data_fns.R**: contains functions for loading data
- **uc_fns.R**: contains functions for estimating the UC model
- **sdr_fns.R**: contains functions for generating term structures of SDRs using simulation

1. **estimate_uc_real_rate.R**
  - Estimates the UC model. Allows to toggle between 1y and 10y rate (main runs both)
  - Saves estimation results to 'results' folder.

2. **figure_1.R**
  - produces and saves Figure 1 (estimates of r*)
  
3. **sdr_uc.R**
  - Calculates term structures of SDRs for UC model. Allows to toggle between 1y and 10y rate (main runs both)
  - Saves term structures to 'results' folder

4. **ar_meanshift.R**
  - Estimates AR model with break in mean and calculate term structures of SDRs. Allows to toggle between 1y and 10y rate (main runs both)
  - Saves term structures to 'results' folder.

5. **ar_learning.R**
  - Estimates AR model with exponential smoothing/learning and calculate term structures of SDRs. Allows to toggle between 1y and 10y rate (main runs both)
  - Saves term structures to 'results' folder.

6. **table_1.R**
  - produces and saves in .tex format Table 1

7. **figure_2.R** 
  - produces and saves Figure 2 (estimates of SDR term structure for 1y and 10y rate)

8. **scc.R**
  - Calculate the estimates of the social cost of carbon.
  - Requires as input term structures of SDRs in 'results' folder (see above scripts) and DICE model damages in 'data' folder.
  - Allows to toggle between DICE-H and DICE-D model (main runs both)
  - Produces Table 2


