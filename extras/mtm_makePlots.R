library(ggplot2)
library(dpf)
library(gridExtra)
data("tempos")
load("extras/mazurka1results.Rdata")
#similar and different plots
performances = rownames(pvec_ml)
similarPerformances = c("Hatto_1988","Chiu_1999",
                        "Kapell_1951","Rubinstein_1966")
similar = which(performances %in% similarPerformances)
similarPlots = vector("list", 4)
differentPerformances = c("Rubinstein_1939","Ohlsson_1999",
                          "Fliere_1977","Fou_1978")
different = which(performances %in% differentPerformances)
differentPlots = vector("list", 4)
for(i in 1:4){
    s = similar[i]
    similarPlots[[i]] = plotStates(tempos[,s+3], similarPerformances[i], tempos$note_onset, unlist(pvec_ml[s,]))
    d = different[i]
    differentPlots[[i]] = plotStates(tempos[,d+3], differentPerformances[i], tempos$note_onset, unlist(pvec_ml[d,]))
}
ggsave("extras/similar_performances.pdf", marrangeGrob(similarPlots, nrow = 2, ncol = 2))
ggsave("extras/different_performances.pdf", marrangeGrob(differentPlots, nrow = 2, ncol = 2))
#raw tempo heatmap
library(heatmaply)
mat = t(tempos[,-c(1:3)])
heatmaply(as.matrix(dist(mat)),k_row=8,k_col=8,symm=TRUE,
          labCol=rep(NA,nrow(mat)), file='extras/raw_tempo_heatmap.pdf')
