# Simulation studies in supporting information

Several simulation experiments are detailed in the manuscript's supporting
information (see Section S3), with the results presented establishing desirable
performance of our estimators in RCTs, observational studies (with positivity
violations or near-violations), against doubly robust estimators, and in
settings with many covariates. R scripts to reproduce these results are located
in the `R/` directory, in particular

* Performance against doubly robust estimators and in settings with many
  covariates is evaluated in `Rcode_Biometrics1.R`, `Rcode_Biometrics2.R`, and
  `Rcode_Biometrics3.R`.
  * `Rcode_Biometrics1.R`:
  * `Rcode_Biometrics2.R`:
  * `Rcode_Biometrics3.R`:
* Performance in RCTs and the two observational study settings is evaluated in
   all other scripts in the `R/` directory; the numbered scripts must be run
   sequentially, and all other scripts contain either helper functions or code
   for creating the plots in the manuscript.
