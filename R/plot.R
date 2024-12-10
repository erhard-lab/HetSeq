library(ggrastr)
library(cowplot)
library(ggrepel)
#' Plot Heterogeneity-seq Classifier Results
#'
#' This plotting functions creates AUC Scatter plots to visualize Classifier Results.
#' @param table Table from HetSeq using the Classifier method.
#' @param highlights A vector of genes to highlight in the plot.
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
PlotClassify <- function(table, highlights=NULL, highlights.color=NULL, highlights.cutoff=NULL, label.cutoff=1.1, density.n=500, point.scale=0.5, xlab="AUC", ylab=bquote(log[2]~FC~(`0h`)), linetype="dashed"){
  plot <- ggplot2::ggplot(table, aes(x=AUC, y=LFC,label=Gene, color=grandR::density2d(x=AUC,y=LFC, n = density.n)))+ggrastr::geom_point_rast(scale=point.scale)+
    viridis::scale_color_viridis_c()+xlab(xlab)+ylab(ylab)+cowplot::theme_cowplot()+
    geom_vline(xintercept = table["baseline",]$AUC, linetype=linetype)+geom_hline(yintercept = 0)+
    theme(legend.position="none")
  
  if(length(highlights)!=length(highlights.color)){
    warning("Number of highlights did not match number of highlight colors. Using default color scale.")
    highlights.color=scales::hue_pal()(length(highlights))
  }
  
  if(length(highlights)!=0){
    plot <- plot + ggrastr::geom_point_rast(data=table[table$Gene%in%highlights & table$AUC > highlights.cutoff,],aes(x=AUC,y=LFC),color=highlights.color, scale=point.scale)
    plot <- plot + ggrepel::geom_label_repel(data=table[table$Gene%in%highlights & table$AUC > label.cutoff,], aes(x=AUC,y=LFC,label=Gene),color=highlights.color)
  }
  plot
}

#' Plot Heterogeneity-seq DoubleML Results
#'
#' This plotting functions creates a Vulcano Plot to visualize DoubleML Results.
#' @param table Table from the Hetseq using the doubleML method.
#' @param p.cutoff Adds a dashed horizontal line at the given adjusted p-value cutoff.
#' @param est.cutoff Adds two dashed vertical lines (+/-) at the given estimate cutoff.
#' @param highlights A vector of genes to highlight in the plot.
#' @param highlights.color A vector of colors for gene highlights.
#' @param highlight.p.cutoff Set a adjusted p-value cutoff below which genes will be highlighted. Defaults to the p.cutoff.
#' @param highlight.est.cutoff Set an estimate cutoff above/below which genes will be highlighted. Defaults to 0.
#' @param label.p.cutoff Set a cutoff at which gene highlights will be labeled. Defaults to the highlight.p.cutoff.
#' @param label.est.cutoff Set a cutoff at which gene highlights will be labeled. Defaults to the highlight.est.cutoff.
#' @param label.repulsion Represents the force parameter of the ggrepel::geom_label_repel() function. Higher values reduce label overlap.
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
PlotDoubleML <- function(table, highlights=NULL, p.cutoff = 0.05, est.cutoff=NULL, highlights.color=NULL, highlight.p.cutoff=NULL, highlight.est.cutoff=NULL, label.p.cutoff=NULL, label.est.cutoff=NULL, label.repulsion = 1, density.n=500, point.scale=0.5, xlab="Estimate", ylab=bquote("-" ~ log[10] ~ FDR), linetype="dashed"){
  plot<- ggplot2::ggplot(table, aes(x=Estimate, y=-log10(p.adj),color=grandR::density2d(Estimate,-log10(p.adj), n=density.n)))+ggrastr::geom_point_rast(scale=point.scale)+
    viridis::scale_color_viridis_c()+xlab(xlab)+ylab(ylab)+theme_cowplot()+
    geom_hline(yintercept = -log10(p.cutoff), linetype=linetype)+
    theme(legend.position="none")
  
  if(length(highlights)!=length(highlights.color)){
    warning("Number of highlights did not match number of highlight colors. Using default color scale.")
    highlights.color=scales::hue_pal()(length(highlights))
  }
  if(is.null(highlight.p.cutoff)){
    highlight.p.cutoff=p.cutoff
  }
  if(is.null(highlight.est.cutoff)){
    highlight.est.cutoff=0
  }
  if(is.null(label.p.cutoff)){
    label.p.cutoff=highlight.p.cutoff
  }
  if(is.null(label.est.cutoff)){
    label.est.cutoff=highlight.est.cutoff
  }
  if(!is.null(est.cutoff)){
    plot<-plot+geom_vline(c(-est.cutoff,est.cutoff), linetype=linetype)
  }
  if(length(highlights)!=0){
    plot <- plot + ggrastr::geom_point_rast(data=table[table$Gene%in%highlights & table$p.adj < highlight.p.cutoff & abs(table$Estimate) > highlight.est.cutoff,],aes(x=Estimate,y=-log10(p.adj)),color=highlights.color, scale=point.scale)
    plot <- plot + ggrepel::geom_label_repel(data=table[table$Gene%in%highlights & table$p.adj < label.p.cutoff & abs(table$Estimate) > label.est.cutoff,], aes(x=Estimate,y=-log10(p.adj),label=Gene),color=highlights.color, force = label.repulsion)
  }
  
  plot
}
  