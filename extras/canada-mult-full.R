library(dpf)
library(batchtools)
library(dplyr)
library(optimr) # cran version
# library(optimrx) # dev version (installed on my laptop)
# Note: optimrx correctly handles inability to find an optimum within multistart(), optimr does not
# fix is to add badval=1e8 to the control list, see optimizer() in my_model.R

data(tempos)
lt = rep(1, nrow(tempos))
source("mazurka-code/multiplicative_model.R")
source("mazurka-code/prior-checks-mult.R")

pfuns = list(
  original = logprior,
  less_var = lessvar,
  igs = logprior_igs,
  uvars = logprior_uvar,
  uprobs = logprior_uprobs)

makeRegistry("mazurka-mult-data", packages=c('dpf','optimr'), 
             source = c("mazurka-code/multiplicative_model.R", "mazurka-code/prior-checks-mult.R"),
             seed = 20181109)
batchMap(optimizer, as.list(select(tempos, -meas_num, -note_onset, -beat)),
         more.args = list(lt=lt, priorfun=logprior))


makeRegistry("mazurka-mult-priors", packages=c('dpf','optimr'), 
             source = c("mazurka-code/multiplicative_model.R", "mazurka-code/prior-checks-mult.R"),
             seed = 20181109)
batchMap(optimizer, priorfun = pfuns,
         more.args = list(lt=lt, perf=tempos[["Richter_1976"]]))

