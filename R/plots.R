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


#' Plots the Music Model states and Kalman smoothed estimates
#'
#' @param theta parameters in same order
#' 
#' @param states vector of states at each time period (not converted)
#' 
#' @param y vector of tempos or dynamics of a single performance
#' 
#' @param onset vector of onsets.
#' 
#' @param model If using the dynamics model put "dynamics".  If using the tempo model put "tempo".
#' 
#' @param priormean Prior Mean array
#' 
#' @param priorvar Prior Variance array
#' 
#' @param title Change the title of the plot
#'
#' @return A plot showing the estimated values of the kalman smoother along with the states.
#' 
#' @examples
#'theta = c(426.69980736, 136.33213703, -11.84256691, -34.82234559, 
#'          439.37886221, 1, 1, 0.84916635, 0.04611644, 0.74119571, 
#'          0.43966082, 0.02116317, 0.24513563, 0.17253254)
#'y = matrix(tempos[,'Richter_1976'], 1)
#'lt = diff(c(tempos$note_onset, 61))
#'pmats = musicModel(lt, theta[1], theta[2:4], theta[5:7], theta[8:14], 
#'                   c(132,0), c(400,10)) # prior means and variances on X_1
#'beam = with(pmats, beamSearch(a0, P0, c(1,0,0,0,0,0,0,0,0,0), dt, ct, Tt, Zt,
#'                              HHt, GGt, y, transMat, 200))
#'states = beam$paths[which.max(beam$weights),]
#'
#'plotStates2(theta,states,y,tempos$note_onset,model="tempo",c(132,0),c(400,10),title="Richter 1976")
#'
#'
#'theta = c(20,10,
#'          12,13,14,
#'          21,22,23,
#'          .3,.2,.1,.5)
#'y = matrix(dynamics[,'Richter_1976'], 1)
#'lt = diff(c(dynamics$note_onset, 61))
#'
#'pmats = musicModeldynamics(lt, theta[1], theta[2], theta[3:5], theta[6:8], theta[9:12], c(20,0,0), c(25,10,0)) 
#'beam = with(pmats, beamSearch(a0, P0, c(1,0,0,0), dt, ct, Tt, Zt,
#'                              HHt, GGt, y, transMat, 200))
#'states = beam$paths[which.max(beam$weights),]
#'kal <- kalman(pmats,states,y)
#'plotStates2(theta,states,y,dynamics$note_onset,model="dynamics",c(20,0,0), c(25,10,0),title="Richter 1976")
#' 
#' @export
plotStates2 <- function(theta, states, y, onset, model, priormean, priorvar,title="Music Plot"){
  lt = diff(c(onset,61))
  
  if(model=="dynamics"){
    mats = musicModeldynamics(lt, theta[1], theta[2:4], theta[5:7], theta[8:11],
                              priormean, priorvar)
    kal <- kalman(pmats,states,y)
    df = data.frame(
      measure = onset, 
      dynamics = c(y), 
      inferred = c(kal$ests), 
      state = factor(
        convert4to2dynamics(states),
        levels=c(1,2),
        labels=c('New Value', 'Smooth Progression')
      )
    )
    myplot <- ggplot2::ggplot(df) + 
      ggplot2::geom_rect(
        data=data.frame(xmin = 33, xmax = 45, ymin = -Inf, ymax = Inf),
        mapping=ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
        fill = 'gray90', color = 'gray90') +
      ggplot2::geom_line(ggplot2::aes(x=measure, y=dynamics), color='black') +
      ggplot2::geom_point(ggplot2::aes(x=measure, y=inferred, color=state)) +
      ggplot2::scale_color_viridis_d() +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = 'bottom', 
                     legend.title = ggplot2::element_blank()) +
      ggplot2::ggtitle(title)
  }
  
  if(model=="tempo"){
    mats = musicModel(lt, theta[1], theta[2:4], theta[5:7], theta[8:14], 
                      priormean, priorvar)
    kal <- kalman(pmats,states,y)
    df = data.frame(
      measure = onset, 
      tempo = c(y), 
      inferred = c(kal$ests), 
      state = factor(
        convert10to4(states),
        levels=c(1,2,3,4),
        labels=c('constant tempo', 'decelerating','accelerating','stress')
      )
    )
    myplot <- ggplot2::ggplot(df) + 
      ggplot2::geom_rect(
        data=data.frame(xmin = 33, xmax = 45, ymin = -Inf, ymax = Inf),
        mapping=ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
        fill = 'gray90', color = 'gray90') +
      ggplot2::geom_line(ggplot2::aes(x=measure, y=tempo), color='black') +
      ggplot2::geom_point(ggplot2::aes(x=measure, y=inferred, color=state)) +
      ggplot2::scale_color_brewer(palette = 'Paired') +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = 'bottom', 
                     legend.title = ggplot2::element_blank()) +
      ggplot2::ggtitle(title)
  }
  return(myplot)
}
