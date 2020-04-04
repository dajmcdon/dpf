# 1. This file estimates parameters for all performances in parallel using
# a cluster. See https://mllg.github.io/batchtools/ for setup instructions
# to run on your own machine or an appropriate managed HPC system.
# Note that, as specified, each performance takes around 5 hours on my laptop.

library(dpf)
library(batchtools)
library(dplyr)
library(optimr) # cran version
# library(optimrx) # dev version (installed on my laptop)
# Note: optimrx correctly handles inability to find an optimum within multistart(), optimr does not
# fix is to add badval=1e8 to the control list, see optimizer() in djm_optimization.R

data(dynamics)
lt = diff(c(tempos$note_onset,61))
source("my_model_dynamics.R")
makeRegistry("my-mazurka-dynamics", packages=c('dpf','optimr'), 
             source = "my_model_dynamics.R",
             seed = 20181109)
batchMap(optimizer, as.list(select(dynamics, -meas_num, -note_onset, -beat)),
         more.args = list(lt=lt))
submitJobs(resources = list(ppn=1, nodes=1, memory='16gb', walltime='24:00:00'))
