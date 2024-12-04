#' Plot Heterogeneity-seq Classifier Results
#'
#' This plotting functions creates AUC Scatter plots to visualize Classifier Results.
#' @param table Table from HetSeq using the Classifier method.
#' @param highlights A list of gene vectors to highlight in the plot.
#' @param highlights.color A vector of colors for gene highlights.
#' @param highlights.cutoff Set a cutoff at which gene highlights should be set.
#' @param label.cutoff Set a cutoff at which gene highlights will be labeled. Values > 1 will lead to no genes being labeled.
#' @param density.n Set granularity of 2d density color.
#' @param point.scale Set point size.
#' @param xlab Set label of the x-axis.
#' @param ylab Set label of the y-axis.
#' @param linetype Set the linetype of the baseline AUC line.
#' @return ggplot object
#' @examples 
#' res=hetseq(method="classify", data, trajectories, score.name = "score", quantiles = c(0.25,0.75), compareGroups = c("Low", "High"), posClass=posClass)
#' plot.classify(res, highlights=list(c("GAPDH", "MYC", "ISG15")))
#' @export
PlotClassify <- function(table, highlights=NULL, highlights.color=Seurat::DiscretePalette(36), highlights.cutoff=NULL, label.cutoff=1.1, density.n=500, point.scale=0.5, xlab="AUC", ylab=bquote(log[2]~FC~(`0h`)), linetype="dashed"){
  plot <- ggplot(table, aes(x=AUC, y=LFC,label=Gene, color=density2d(x=AUC,y=LFC, n = density.n)))+geom_point_rast(scale=point.scale)+
    scale_color_viridis_c()+xlab(xlab)+ylab(ylab)+theme_cowplot()+
    geom_vline(xintercept = table["baseline",]$AUC, linetype=linetype)+geom_hline(yintercept = 0)+
    theme(legend.position="none")
  
  for(i in c(1:length(highlights))){
    plot <- plot + geom_point_rast(data=table[table$Gene%in%highlights[[i]] & table$AUC > highlights.cutoff,],aes(x=AUC,y=LFC),color=highlights.color[i], scale=point.scale)
    plot <- plot + geom_label_repel(data=table[table$Gene%in%highlights[[i]] & table$AUC > label.cutoff,], aes(x=AUC,y=LFC,label=Gene),color=highlights.color[i])
  }
  
  plot
}

#' Plot Heterogeneity-seq DoubleML Results
#'
#' This plotting functions creates a Vulcano Plot to visualize DoubleML Results.
#' @param table Table from the Hetseq using the doubleML method.
#' @param highlights A list of gene vectors to highlight in the plot.
#' @param highlights.color A vector of colors for gene highlights.
#' @param highlights.cutoff Set a cutoff at which gene highlights should be set.
#' @param label.cutoff Set a cutoff at which gene highlights will be labeled. Values > 1 will lead to no genes being labeled.
#' @param density.n Set granularity of 2d density color.
#' @param point.scale Set point size.
#' @param xlab Set label of the x-axis.
#' @param ylab Set label of the y-axis.
#' @param linetype Set the linetype of the p-value and estimate cutoff line.
#' @return ggplot object
#' @examples 
#' hetseq(method="doubleML", data, trajectories, score.name = "score", quantiles = c(0.25,0.75), compareGroups = c("Low", "High"))
#' plot.doubleml(res, highlights=list(c("GAPDH", "MYC", "ISG15")))
#' @export
PlotDoubleML <- function(table, highlights=NULL, p.cutoff = 0.05, est.cutoff=NULL, highlights.color=Seurat::DiscretePalette(36),highlight.cutoff=NULL, label.cutoff=1.1, density.n=500, point.scale=0.5, xlab="Estimate", ylab=bquote("-" ~ log[10] ~ FDR), linetype="dashed"){
  plot<- ggplot(table, aes(x=Estimate, y=-log10(p.adj),color=density2d(Estimate,-log10(p.adj), n=density.n)))+geom_point_rast(scale=point.scale)+
    scale_color_viridis_c()+xlab(xlab)+ylab(ylab)+theme_cowplot()+
    geom_hline(yintercept = -log10(p.cutoff), linetype=linetype)+
    theme(legend.position="none")
  if(!is.null(est.cutoff)){
    plot<-plot+geom_vline(c(-est.cutoff,est.cutoff), linetype=linetype)
  }
  for(i in c(1:length(highlights))){
    plot <- plot + geom_point_rast(data=table[table$Gene%in%highlights[[i]] & table$p.adj > p.cutoff,],aes(x=Estimate,y=p.adj),color=highlights.color[i], scale=point.scale)
    plot <- plot + geom_label_repel(data=table[table$Gene%in%highlights[[i]] & table$p.adj > label.cutoff,], aes(x=Estimate,y=p.adj,label=Gene),color=highlights.color[i])
  }
  
  plot
}
  