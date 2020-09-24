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
# fix is to add badval=1e8 to the control list, see optimizer() in my_model.R

data(tempos)
lt = diff(c(tempos$note_onset,61))
source("my_model-rev1.R")
source("prior-checks.R")

pfuns = list(
  original = logprior,
  less_var = lessvar,
  igs = logprior_igs,
  uvars = logprior_uvar,
  uprobs = logprior_uprobs)

makeRegistry("mazurka-prior-checks", packages=c('dpf','optimr'), 
             source = c("my_model-rev1.R", "prior-checks.R"),
             seed = 20181109)
batchMap(optimizer, priorfun = pfuns,
         more.args = list(lt=lt, perf=tempos["Richter_1976"]))
#submitJobs(resources = list(ppn=1, nodes=1, memory='16gb', walltime='24:00:00'))
