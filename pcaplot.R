# PCA
pcaplot <- function(counts, colorvec, dotvec){
  storelist <- list()
  countsfilt <- counts[rowSums(counts) >= 10,]
  tcountsfilt = t(countsfilt)
  tcountsfilt = data.frame(tcountsfilt)
  tcountsfilt["Color"] <-  rownames(tcountsfilt)
  dim(tcountsfilt)
  
  dfx <- tcountsfilt[,-c(dim(tcountsfilt)[2])]
  PC <- prcomp(dfx, scale. = T)
  PCi <- data.frame(PC$x,Color=tcountsfilt$Color)
  percentage <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
  percentage <- paste( colnames(PCi), "(", paste( as.character(percentage), "%", ")", sep="") )
  theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
  
  out <- ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
    theme + xlab(percentage[1]) + ylab(percentage[2])+
    geom_point(size=4,alpha=1,aes(shape=Color))+
    scale_color_manual(values = colorvec)+
    scale_shape_manual(values=dotvec) + theme_bw()
  storelist[["pca"]] <- out
  storelist[["countsfilt"]] <- countsfilt
  storelist[["transposed_countsfilt"]] <- tcountsfilt
  return(storelist)
}
