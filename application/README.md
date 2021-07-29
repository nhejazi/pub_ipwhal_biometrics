# NHEFS study application

The main manuscript demonstrates the use of our proposed IPW estimators by
application to the NHEFS data in Section 6. The R scripts required to reproduce
the results of this empirical illustration are located in the `R/` directory,
each numbered and with extensive code comments. The analysis is performed
entirely within the script `01_ipw_analysis.R`, and the results summarized
entirely by `02_summarize.R`. When replicating these data analysis results, note
the [`renv` virtual
environment](https://rstudio.github.io/renv/articles/renv.html) set up in this
directory to ensure that package versions used are kept consistent.
