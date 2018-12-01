library(dpf)
library(batchtools)
library(dplyr)
library(optimr) # cran version
# library(optimrx) # dev version (installed on my laptop)
# Note: optimrx correctly handles inability to find an optimum within multistart(), optimr does not
# fix is to add badval=1e8 to the control list, see optimizer() in djm_optimization.R

data(tempos)
lt = diff(c(tempos$note_onset,61))
source("mazurka_code/optimization_new_prior.R")
makeRegistry("mazurka3", packages=c('dpf','optimr'), 
             source = "mazurka_code/optimization_new_prior.R",
             seed = 20180918)
batchMap(optimizer, as.list(select(tempos, -meas_num, -note_onset, -beat)),
         more.args = list(lt=lt))
submitJobs(resources = list(ppn=1, nodes=1, memory='32gb', walltime='24:00:00'))
