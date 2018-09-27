library(ggplot2)
library(dpf)
library(gridExtra)
data("tempos")
load("extras/mazurka3results.Rdata")
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
          labCol=rep(NA,nrow(mat)), file='extras/raw_tempo_heatmap.pdf', height = 800, width = 1280)
#every performane plot
plots = vector("list", 4)
for(i in 1:nrow(pvec_ml)){
    plots[[i]] = plotStates(tempos[,i+3], rownames(pvec_ml)[i], tempos$note_onset, unlist(pvec_ml[i,]))
}
ggsave("extras/all_performances3.pdf", marrangeGrob(plots, nrow = 2, ncol = 2))
#density plots
library(ggplot2)
library(dpf)
library(gridExtra)
load("manuscript/mazurka2results.Rdata")
load('extras/ClusterLabels.Rdata')
pvec = pvec_ml[-14,]#bad_perf which was removed before clustering
pvec$clust = as.factor(clusts)
plots = vector("list", ncol(pvec)-1)
i = 1
for(p in colnames(pvec)){
    if(p != "clust"){
        local({
            p <- p
            plots[[i]] <<- ggplot(pvec, aes(x = eval(parse(text = p)), fill = clust))+geom_density(alpha = 0.5)+xlab(p)
        })
        i = i + 1
    }
}
ggsave("extras/cluster_density.pdf", marrangeGrob(plots, nrow = 2, ncol = 2, newpage = FALSE))
