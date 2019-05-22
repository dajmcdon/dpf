# Discrete Particle Filtering

This package is intended to estimate Markov switching state-space models using maximum likelihood _fast_.

This repo contains both the `R` package and the manuscript itself.

Our implementation of the Markov switching model allows the parameter matrices of a linear Gaussian state-space model to depend on previous discrete switch states. In other words, the hidden Markov (switch) state determines which parameter
matrices govern the evolution of the system. If the  matrices  depend on a vector of unknown parameters $\theta$, one needs to jointly estimate the switch path $S_1,\ldots,S_n$, the continuous hidden states $X_1,\ldots,X_n$, and the parameter vector.


## Why this package?

Most available packages on CRAN are for state-space models only (not switching) or for hidden Markov models only (without continuous hidden states). Without switch states, the Kalman filter [@Kalman1960] returns the likelihood of a parameter vector $\theta$, and the Kalman smoother gives the estimates of the hidden continuous states conditional on all the observed data. In the case of an HMM, the Viterbi algorithm is analogous to the Kalman smoother, producing the most likely path conditional on all the data.

However, when combined, computing the likelihood for $\theta$ or producing the most likely path becomes exponentially hard: given 2 possible values for $S_i$, there are $2^n$ possible sequences $(S_1,\ldots,S_n)$. Maximizing the likelihood for $\theta$ then means evaluating the Kalman filter at each of these paths _for each proposed value of $\theta$_.

## Contribution

The `dpf` package provides a lightweight, fast, greedy method for calculating the likelihood of $\theta$. It is implemented in C++ using `Rcpp` and takes advantage of an accurate algorithm for finding the most likely state sequence proposed in [@FearnheadClifford2003].

We provide 3 main functions which are useful for generic switching state space models:

1. `kalman()` computes the standard Kalman filter and smoother given a collection of parameter matrices which are allowed to be time dependent. It returns the likelihood, means and variances for the filter and smoother distributions, and predicted values for the observations.

2. `beamSearch()` implements the Greedy HMM algorithm. For any $N$, it attempts to compute the $N$ most likely sequences $(S_1,\ldots,S_n)$ given a collection of parameter matrices which are time and state dependent. The result is this collection of paths as well as probabilities associated with each of the $N$ paths. Using the path with the highest weight in `kalman()` then corresponds to a frequentist evaluation of a particular parameter $\theta$ while sampling proportional to the weights allows for approximate Bayesian inference (ignoring the $s^N-N$ paths which presumably have negligible likelihood).

3. `getLogLike()` computes only the components of the Kalman filter necessary to evaluate the likelihood.

For a generic switching state-space model whose evolution matrices depend on an unknown vector $\theta$, one would use `beamSearch()` and `getLogLike()` combined with any standard optimization method to produce an ML estimate $\widehat{\theta}$. Then, for $\widehat{\theta}$, one would use `beamSearch()` and `kalman()` to return the ML estimates for $(S_1,\ldots,S_n)$ and $(X_1,\ldots,X_n)$.


## Switching models for musical tempo decisions

The remaining functions in `dpf-package` are created for our manuscript which proposes a switching model for the decisions made by classical music performers. Our model essentially imagines that musicians have an intentional tempo for a piece of music which is partially determined by a Markov chain on 4 states: constant tempo, acceleration, deceleration, and stress. While our working manuscript provides the details, we here illustrate the `beamSearch()` function on a recording of Chopin's Mazurka Op. 68 No. 3 made by the famous Russian pianist Nina Milkina in 1970. We estimated $\theta$ offline. The `musicModel()` function creates the appropriate parameter matrices in the above equation given $\theta$. 

These data are included with the package.


## Installation

To install and simultaneously build the vignette run

```
devtools::install_github(‘dajmcdon/dpf’, build_opts = c("--no-resave-data", "--no-manual"))
```

