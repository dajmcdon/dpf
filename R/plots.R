#' Plots the observed tempo and estimated discrete states given parameters.
#'
#' @param performance The named performance, e.g. 'Richter_1976'. See \code{names(tempos)}.
#' @param params A vector of parameters of length 14. No checks are performed.
#' @param y A vector of tempos.
#' @param onset A vector of note onset times.
#' @param particleNumber Number of particles for \code{beamSearch()}. Default is 200.
#' @param initialMean Mean for the first note in constant tempo. Length 2. Default is (132,0).
#' @param initialVariance Variance for the first note in constant tempo. Length 2. Default is (400,0).
#'
#' @export
#' @examples
#' params = c(426.69980736, 136.33213703, -11.84256691, -34.82234559, 
#'           439.37886221, 1, 1, 0.84916635, 0.04611644, 0.74119571, 
#'           0.43966082, 0.02116317, 0.24513563, 0.17253254)
#' data(tempos)
#' y = tempos[,'Richter_1976']
#' onset = tempos$note_onset
#' plotStates('Richter_1976', params, y, onset)
plotStates <- function(performance, params, y, onset,
                       particleNumber = 200, initialMean = c(132,0), 
                       initialVariance = c(400,10)){
    if(is.list(params)) params = unlist(params)
    y = matrix(y, nrow = 1)
    lt = diff(c(onset, 61))
    mats = musicModel(lt, params[1], params[2:4], params[5:7], params[8:14], 
                      initialMean, initialVariance)
    bs = beamSearch(mats$a0, mats$P0, c(1,rep(0,9)), mats$dt, mats$ct, 
                    mats$Tt, mats$Zt, mats$HHt, mats$GGt, y, 
                    mats$transMat, particleNumber)
    bestpath = bs$paths[which.max(bs$weights),]
    kal = kalman(mats, bestpath, y)
    df = data.frame(
      measure = onset, 
      tempo = c(y), 
      inferred = c(kal$ests), 
      state = factor(
        convert10to4(bestpath), 
        labels=c('constant tempo', 'decelerating','accelerating','stress')
        )
    )
    ggplot2::ggplot(df) + 
      ggplot2::geom_rect(
        data=data.frame(xmin = 33, xmax = 45, ymin = -Inf, ymax = Inf),
        mapping=ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
        fill = 'gray90', color = 'gray90') +
      ggplot2::geom_line(ggplot2::aes(x=measure, y=tempo), color='black') +
      ggplot2::geom_point(ggplot2::aes(x=measure, y=inferred, color=state)) +
      ggplot2::scale_color_brewer(palette = 'Spectral') +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = 'bottom', 
                     legend.title = ggplot2::element_blank()) +
      ggplot2::ggtitle(performance)
}
