#' dpf: Discrete particle filtering
#'
#' The dpf package provides three categories of important functions:
#' Kalman filters, a greedy discrete particle filter (beam search), and 
#' functions for applying these to musical tempo analysis.
#' 
#' For a concise description, see the package vignette.
#' 
#' The most recent version of this package is available on github, and installable
#' with
#' 
#' \code{devtools::install_github('dajmcdon/dpf', build_vignettes=TRUE)}
#' 
#' @section Kalman filters
#'
#' \code{\link{kalman}}
#' \code{\link{getLogLike}} 
#'
#' @section Beam search
#' 
#' \code{\link{beamSearch}}
#' 
#' @section Music tempo analysis
#' 
#' \code{\link{musicModel}}
#' \code{\link{convert10to4}}
#' 
#' 
#'
#' @docType package
#' @name dpf
NULL

#' @useDynLib dpf
#' @import ggplot2
#' @importFrom usedist dist_make
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats rgamma
NULL

#' Dynamics of 46 Chopin piano recordings
#'
#' A dataset containing the dynamics of 46 recordings of Chopin's Mazurka
#' Op. 68, No. 3
#'
#' @format A data frame with 231 rows and 49 variables:
#' \describe{
#'   \item{meas_num}{The number of the measure in the score}
#'   \item{beat}{The beat (within the measure) at which the note occurs}
#'   \item{note_onset}{The relative onset time of the note. This is calculated as (beat-meas_num)/3 + meas_num.}
#'   \item{remaining columns}{The loudness (in dB) of each note as played in the particular performance. Column names give the name of the performer and year of the recording.}
#' }
#' 
#' @source \url{http://www.mazurka.org.uk}
"dynamics"

#' Tempos of 46 Chopin piano recordings
#'
#' A dataset containing the relative tempos of 46 recordings of Chopin's Mazurka
#' Op. 68, No. 3
#'
#' @format A data frame with 231 rows and 49 variables:
#' \describe{
#'   \item{meas_num}{The number of the measure in the score}
#'   \item{beat}{The beat (within the measure) at which the note occurs}
#'   \item{note_onset}{The relative onset time of the note. This is calculated as (beat-meas_num)/3 + meas_num.}
#'   \item{remaining columns}{The instantaneous tempo (in beats per minute) of each note as played in the particular performance. Column names give the name of the performer and year of the recording. These are calculated from the raw data by dividing the actual duration of the note by the notated duration. That is, a note which should take up one beat in the score and lasts for .5 seconds is recorded as 120 bpm.}
#' }
#' 
#' @source \url{http://www.mazurka.org.uk}
"tempos"

#' Recordings
#'
#' A dataset containing the performers and years of performance of 46 recordings of Chopin's Mazurka Op. 68, No. 3
#'
#' @format A data frame with 46 rows and 2 variables:
#' \describe{
#'   \item{performer}{The last name of the pianist}
#'   \item{year}{The year in which the recording was made}
#' }
#' 
#' @source \url{http://www.mazurka.org.uk}
"recordings"

