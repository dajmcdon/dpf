#' @export
plotStates <- function(performance, performer, onset, params, particleNumber = 200, initialMean = c(132,0), initialVariance = c(400,10)){
    if(is.list(params)) params = unlist(params)
    lt = diff(c(onset, 61))
    y = matrix(performance, nrow = 1)
    mats = yupengMats(lt, params[1], params[2:4], params[5:7], params[8:11], initialMean, initialVariance)
    bs = beamSearch(mats$a0, mats$P0, c(1,0,0,0,0,0,0,0), mats$dt, mats$ct, mats$Tt, mats$Zt,
                    mats$HHt, mats$GGt, y, mats$transMat, particleNumber)
    bestpath = bs$paths[which.max(bs$weights),]
    kal = kalman(mats, bestpath, y)
    df = data.frame(measure = tempos$note_onset, tempo = performance, 
                    inferred = c(kal$ests), state = convert8to4(bestpath))
    ggplot2::ggplot(df, ggplot2::aes(x=measure, y=tempo)) + ggplot2::ylim(0,max(df$tempo)) +
        ggplot2::geom_rect(ggplot2::aes(xmin = 33, xmax = 45, ymin = 0, ymax = max(df$tempo),
                      fill = 'gray90', color = 'gray90'))+
        ggplot2::geom_line(ggplot2::aes(y=df$tempo), color='black')+
        ggplot2::geom_point(ggplot2::aes(y=df$inferred), 
                   color = c("blue","red","green","orange")[df$state])+
        #scale_color_manual(values = c("blue","red","green","orange"))+
        ggplot2::theme(legend.position = 'none', legend.title = ggplot2::element_blank())+
        ggplot2::ggtitle(performer)
}