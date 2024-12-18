#' Plot Heterogeneity-seq Classifier Results
#'
#' This plotting functions creates AUC Scatter plots to visualize Classifier Results.
#' @param table Table from HetSeq using the Classifier method.
#' @param highlights A vector of genes to highlight in the plot.
#' @param highlights.color A vector of colors for gene highlights.
#' @param plot.baseline Inserts a vertical line at the baseline AUC value (= AUC of the classifier with basefeatures but no further gene info).
#' @param auc.cutoff Inserts a vertical line at the cutoff AUC value.
#' @param density.color Color data points by density. Default is TRUE.
#' @param density.n Set granularity of 2d density color.
#' @param point.scale Set point size.
#' @param xlab Set label of the x-axis.
#' @param ylab Set label of the y-axis.
#' @param linetype Set the linetype of the baseline AUC line.
#' @return ggplot object
#' @export
PlotClassify <- function(table, highlights=NULL, highlights.color=NULL, auc.cutoff = NULL, plot.baseline=TRUE, density.color = TRUE, density.n=500, point.scale=0.5, xlab="AUC", ylab=bquote(log[2]~FC~(`0h`)), linetype="dashed"){
  AUC<-LFC<-Gene<-NULL
  
  if(density.color){
    plot <- ggplot(table, aes(x=AUC, y=LFC,label=Gene, color=grandR::density2d(x=AUC,y=LFC, n = density.n)))
  } else {
    plot <- ggplot(table, aes(x=AUC, y=LFC,label=Gene))
  }
  plot <-plot+ggrastr::geom_point_rast(scale=point.scale)+
    scale_color_viridis_c()+xlab(xlab)+ylab(ylab)+cowplot::theme_cowplot()+
    geom_vline(xintercept = table["baseline",]$AUC, linetype=linetype)+geom_hline(yintercept = 0)+
    theme(legend.position="none")
  
  if(length(highlights)!=length(highlights.color)){
    warning("Number of highlights did not match number of highlight colors. Using default color scale.")
    highlights.color=scales::hue_pal()(length(highlights))
  }
  
  if(is.null(auc.cutoff)){
    plot <- plot + geom_vline(xintercept = auc.cutoff, linetype=linetype)
  }
  
  if(length(highlights)!=0){
    plot <- plot + ggrastr::geom_point_rast(data=table[table$Gene%in%highlights,],aes(x=AUC,y=LFC),color=highlights.color, scale=point.scale)
    plot <- plot + ggrepel::geom_label_repel(data=table[table$Gene%in%highlights,], aes(x=AUC,y=LFC,label=Gene),color=highlights.color)
  }
  if(plot.baseline){
    if(is.na(table["baseline",]$Gene)){
      warning("baseline info not in the table.")
    } else {
      plot <- plot + geom_vline(xintercept = table["baseline",]$AUC, linetype=linetype)
    }
  }
  plot
}

#' Plot Heterogeneity-seq DoubleML Results
#'
#' This plotting functions creates a Vulcano Plot to visualize DoubleML Results.
#' @param table Table from the Hetseq using the doubleML method.
#' @param highlights A vector of genes to highlight in the plot.
#' @param p.cutoff Adds a dashed horizontal line at the given adjusted p-value cutoff.
#' @param est.cutoff Adds two dashed vertical lines (+/-) at the given estimate cutoff.
#' @param highlights.color A vector of colors for gene highlights.
#' @param label.repulsion Represents the force parameter of the ggrepel::geom_label_repel() function. Higher values reduce label overlap.
#' @param density.color Color data points by density. Default is TRUE.
#' @param density.n Set granularity of 2d density color.
#' @param point.scale Set point size.
#' @param xlab Set label of the x-axis.
#' @param ylab Set label of the y-axis.
#' @param linetype Set the linetype of the p-value and estimate cutoff line.
#' @return ggplot object
#' @export
PlotDoubleML <- function(table, highlights=NULL, p.cutoff = 0.05, est.cutoff=NULL, highlights.color=NULL, label.repulsion = 1, density.color=TRUE, density.n=500, point.scale=0.5, xlab="Estimate", ylab=bquote("-" ~ log[10] ~ FDR), linetype="dashed"){
  Estimate<-p.adj<-Gene<-NULL
  
  if(density.color){
    plot <- ggplot(table, aes(x=Estimate, y=-log10(p.adj),color=grandR::density2d(Estimate,-log10(p.adj), n=density.n)))
  } else {
    plot <- ggplot(table, aes(x=Estimate, y=-log10(p.adj)))
  }
                
  plot <- plot+ggrastr::geom_point_rast(scale=point.scale)+
    scale_color_viridis_c()+xlab(xlab)+ylab(ylab)+cowplot::theme_cowplot()+
    geom_hline(yintercept = -log10(p.cutoff), linetype=linetype)+
    theme(legend.position="none")
  
  if(length(highlights)!=length(highlights.color)){
    warning("Number of highlighted genes did not match number of highlight colors. Using default color scale.")
    highlights.color=scales::hue_pal()(length(highlights))
  }

  if(!is.null(est.cutoff)){
    plot<-plot+geom_vline(c(-est.cutoff,est.cutoff), linetype=linetype)
  }
  if(!is.null(p.cutoff)){
    plot<-plot+geom_hline(-log10(p.cutoff), linetype=linetype)
  }
  if(length(highlights)!=0){
    plot <- plot + ggrastr::geom_point_rast(data=table[table$Gene%in%highlights,],aes(x=Estimate,y=-log10(p.adj)),color=highlights.color, scale=point.scale)
    plot <- plot + ggrepel::geom_label_repel(data=table[table$Gene%in%highlights,], aes(x=Estimate,y=-log10(p.adj),label=Gene),color=highlights.color, force = label.repulsion)
  }
  
  plot
}
  