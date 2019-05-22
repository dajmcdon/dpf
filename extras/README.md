# Reproducing parameter estimates.

This short file describes how to reproduce parameter estimates from the data included with our `dpf` package.

## Preliminaries

* Insure that `dpf` has been installed properly.
* This also requires the `R` packages `dplyr` and `optimr` (or `optimrx`, the development version on Github)
* Finally, you must install and configure the `batchtools` package. See https://mllg.github.io/batchtools/ for setup instructions to run on your own machine or an appropriate managed HPC system.
* Note that estimation for each performance required around 5 hours on a 2018 MacBook Pro.

## Workflow

1. Open an `R` session in this directory. The file "my_model.R" should be in the same location as the `R` session.
2. Source the file "cluster_analysis_my_model.R". This will not actually submit the jobs to the cluster (or laptop). The final commented line of that file submits for my configuration based on walltime and processor limitations. See the `batchtools` documentation for appropriate changes.
3. Once the analysis has completed, source "collect_cluster_results.R". This will generate the dataframe "mazurkaResults.Rdata" included in this directory.