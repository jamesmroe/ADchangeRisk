![LCBC logo](docs/images/LCBC_logo.png)


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13861219.svg)](https://doi.org/10.5281/zenodo.13861219)


Repository associated with the paper

**Brain change trajectories in healthy adults correlate with Alzheimer's related genetic variation and memory decline across life**
================

[Preprint](https://www.biorxiv.org/content/10.1101/2023.10.09.559446v1)\
Status: Accepted (Nature Communications)

<br>
<hr>

**The following scripts will run with simulated input:**

* **_01a-ADchangeRisk_GAMMprocPRS_simulate.r_**  

`R` script to run simulated GAMM trajectory modelling example and PRS-AD association test procedure with simulated data.\
_Instructions:_ Open Rproject file `ADchangeRisk.Rproj`. Script will run with the provided simulated data.

<br>
<hr>

**The following scripts will run with provided summary-level data as input:**

* **_02-GAMMprocPRS_merge.r_**  

_Instructions:_  
Set `simulated = 0` to load the provided summary-level source data underlying Fig. 1D-E and Fig. 2 (PRS-AD tests), run multiple testing correction, and reproduce plots.\
Set `simulated = 1` to visualize the simulated PRS-AD assocation tests from `01a-ADchangeRisk_GAMMprocPRS_simulate.r`.

* **_plotUnivariate.R_**  

_Instructions:_  
set to `analysis_sample` to "main" or "replication" to plot the univariate PRS-AD results for that sample. This script is similar to the one above with added functionality for plotting the replication results.

* **_plotMultivariate.R_**  

_Instructions:_  
set to `analysis_sample` to "main" or "replication" to plot the multivariate PRS-AD results for that sample.

<br>
<hr>

**The following scripts reproduce paper results but require access to restricted individual-level data and are _not executable_:**

* **_01-GAMMprocPRS.r_**

`R` script to reproduce GAMM trajectory modelling and PRS-AD association tests in paper.
(02-ADchangeRisk_GAMMprocPRS_merge.r visualizes output)


* **_03-prepSlopesADNI.r_**

`R` script to prepare individual-specific slopes in longitudinal AD-control data (ADNI) for machine learning with xgboost.

* **_04-XGBoostCV.r_**

`R` script to run cross-validation procedure with xgboost.

* **_05-runMLmultivariate.r_**

`R` script to run ML analysis in AD-control data (ADNI), apply AD-control weights to healthy adult lifespan brain change estimates, and test PRS-AD-change associations.
Performs multivariate analysis using brain features with accelerated change in AD to test PRS-AD-brain change assocations in healthy adults.

* **_06-memorySlopes.r_**

`R` script to run memory change analysis in LCBC healthy adult lifespan data

