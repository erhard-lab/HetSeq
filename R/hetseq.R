#' @import ggplot2
#' @import stats
#' @import Seurat
NULL
#> NULL



unregister_pp <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


#' Heterogeneity-seq wrapper
#'
#' Wrapper function for all Heterogeneity-seq functions.
#' 
#' Heterogeneity-seq uses intercellular heterogeneity to identify modulators of cellular response to perturbations. Three approaches to identify these factors are currently available:
#' 
#' - test: Differential Gene Expression testing between cells from two response groups.
#'
#' - classify: Predicting cellular outcome by expression of single genes (+informative features) reveals genes with strong predictive capabilities (high AUC) as potential pathway modulators.
#'
#' - doubleML: A strict predictor that utilizes Causal Inference (doubleML) to distinguish causal factors from simply correlating genes. Genes with high estimated effects on the outcome and high significance are identified as potential pathway modulators.
#' 
#' @param method The method to run Heterogeneity-seq. Calls hetseq.test, hetseq.classify or hetseq.doubleML.
#' @param ... Parameters given to the chosen Hetseq method. See respective help pages.
#' @return Table of Heterogeneity-seq results.
#' @examples 
#' # Testing for differential gene expression: 
#' labels = c("Low","Middle","High")
#' group = cut(data$score[trajectories[,ncol(trajectories)]],breaks=c(-10,0,5,10),labels=labels)
#' tab = FetchData(data,vars = genes,cells = trajectories[,1],layer = "data")
#' mat = t(as.matrix(tab))[,trajectories[,1]]
#' hetseq(method="test", mat, group=="High",group=="Low")
#' 
#' # Using the Classifier to predict cellular outcomes by expression of single genes.
#' # Genes with strong predictive capabilities (high AUC) are candidate genes for pathway modulators.
#' hetseq(method="classify", data, trajectories, score.name = "score")
#' 
#' # Using a Classifier and Causal Inference (doubleML) to predict cellular outcomes by expression of single genes.
#' # Effect of genes on the outcome are estimated to distinguish genes with causal effects from correlating genes.
#' hetseq(method="doubleML", data, trajectories, score.name = "score")
#' @export
Hetseq = function(method=c("test", "classify", "doubleML"), ...){
    if (method == "test") {
      return(HetseqTest(...))
    } else if (method == "classify") {
      return(HetseqClassify(...))
    } else if (method == "doubleML") {
      return(HetseqDoubleML(...))
    } else {
      stop("Invalid hetseq method specified. Please use 'test', 'classify', or 'doubleML'.")
    }
}

#' Heterogeneity-seq: Testing for differential gene expression
#'
#' Testing for differential gene expression in predecessors of treated cells
#' @param mat Gene expression matrix of control cells
#' @param A Vector of cells (columns) in positive class
#' @param B Vector of cells (columns) in negative class
#' @return Table of log2FC, p-values and adjusted p-values of differentially expressed genes
#' @examples 
#' labels = c("Low", "Middle", "High")
#' group = cut(data$score[trajectories[,ncol(trajectories)]],breaks=c(-10,0,5,10),labels=labels)
#' tab = FetchData(data,vars = genes,cells = trajectories[,1],layer = "data")
#' mat = t(as.matrix(tab))[,trajectories[,1]]
#' hetseq(method="test", mat, group=="High",group=="Low")
#' 
#' @export
HetseqTest = function(mat,A,B) {
  ps=sapply(1:nrow(mat),function(i) wilcox.test(mat[i,A],mat[i,B])$p.value)
  eff=sapply(1:nrow(mat),function(i) (log1p(mean(expm1(mat[i,A])))-log1p(mean(expm1(mat[i,B]))))/log(2))
  expr=sapply(1:nrow(mat),function(i) log1p(mean(expm1(mat[i,A|B]))))
  re = data.frame(Expr=expr,LFC=eff,P=ps, Q = p.adjust(ps,method="BH"))
  rownames(re) = rownames(mat)
  re
}

#' Heterogeneity-seq: Classifying cellular response by gene expression values
#'
#' Classifying the cellular response of control cells using single gene expression (+ informative features) to identify features with the strongest predictive capabilities.
#' @param object Seurat object
#' @param trajectories Matrix of cell-cell trajectories. Columns represent time points, rows represent trajectories of connected cells over time points.
#' @param score.group A named vector of response groups. Names represent cells, the values represent the score groups. If no score.group is set, use score.name and quantiles parameters must be set to define score groups.
#' @param score.name The name of a numeric Seurat meta data column, which will be used to calculate score groups. Only used if no score.group is given.
#' @param quantiles Thresholds of the score.name meta data to define 3 response groups. Low, Middle, High.
#' @param compareGroups Which score groups to test. Default: Low vs. High
#' @param posClass Define the positive Class for classification.
#' @param basefeatures Additional informative features to include in the classification. Must be meta data available in the Seurat object.
#' @param genes Vector of genes to test.
#' @param assay The name of the Seurat assay to perform Heterogeneity-seq on. If NULL, the default assay will be used.
#' @param split Set a training-test data split. Must be in [0,1]
#' @param kernel The kernel for the SVM. linear, polynomial, radial or sigmoid. Default: radial. 
#' @param cross Number of cross-validations.
#' @param num_cores The number of cores used in parallel processing.
#' @return Table of log2FC and AUC values for each gene and an additional AUC value for the baseline features.
#' @examples 
#' # Using meta data in the Seurat object called "score" to automatically define response groups
#' HetseqClassify(data, trajectories, score.name = "score")
#' 
#' # Use a manually defined grouping by handing a named vector to score.group.
#' # data$infectionState contains the infection state of all cells and all time points. 
#' # Get the infectionState of all cells from the last time point in the trajectories.
#' groups <- data$infectionState[trajectories[,ncol(trajectories)],]
#' # Adapt compareGroups to fit the low and high group names of the manual definition.
#' HetseqClassify(data, trajectories, score.group = groups, compareGroups = c("Bystander", "High"))
#' @export
HetseqClassify<-function(object, trajectories, score.group = NULL, score.name=NULL,  quantiles = c(0.25,0.75), compareGroups = c("Low", "High"), posClass=NULL, basefeatures=NULL, genes=NULL, assay=NULL, split=NULL, kernel="radial",  cross=10, num_cores = 1){
  #library(parallel)
  #library(doParallel)
  #library(e1071)
  #library(pROC)
  
  unregister_pp()
  
  if(!is.null(split) && (split > 1 | split < 0)) stop("split must be either NULL or in [0,1]")
  
  if(is.null(score.group)){
    if(is.null(score.name)){
      stop("Either provide a score.group or a score.name (Seurat metadata) on which to calculate the score.group")
    }
    br1=quantile(object[[score.name]][trajectories[,length(colnames(trajectories))],], probs = quantiles)[1]
    br2=quantile(object[[score.name]][trajectories[,length(colnames(trajectories))],], probs = quantiles)[2]
    score.group = cut(object[[score.name]][trajectories[,length(colnames(trajectories))],],breaks=c(-10,br1,br2,10),labels=c("Low","Middle","High"))
    object$Trajectory = "none"
    for(x in c(1:length(colnames(trajectories)))){
      object$Trajectory[trajectories[,x]] = as.character(score.group)
    }
    object$Trajectory = factor(object$Trajectory,levels=c("none",levels(score.group)))
  } else {
    object$Trajectory = "none"
    for(x in c(1:length(colnames(trajectories)))){
      object$Trajectory[trajectories[,x]] = as.character(score.group)
    }
    object$Trajectory <- as.factor(object$Trajectory)
  }
  object <- subset(object, cells = trajectories[,1], subset = Trajectory %in% compareGroups)
  object$Trajectory <- droplevels(object$Trajectory)
  print(table(object$Trajectory))
  
  
  if(is.null(basefeatures)){
    basefeatures <- c("uninformative")
    object$uninformative <- 1
  }
  
  if(is.null(split)){
    train_idx <- Cells(object)
    test_idx <-  Cells(object)
  } else {
    train_idx <- sample(1:ncol(object), split * ncol(object))
    test_idx <- setdiff(c(1:ncol(object)), train_idx)
  }

  # Register parallel backend
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  if(is.null(assay)){
    assay=DefaultAssay(object)
  }
  DefaultAssay(object) <- assay
  if(is.null(genes)) genes <- rownames(object)
  df_traj<-object@meta.data[train_idx,"Trajectory"]
  if(is.null(posClass)) posClass = names(table(df_traj))[2]
  y_numeric <- as.numeric(df_traj == posClass)
  model <- e1071::svm(x = FetchData(object,cells=train_idx, vars=basefeatures), y = y_numeric, kernel = kernel, cross = cross)
  predictions <- predict(model, newdata = FetchData(object, cells=test_idx, vars=basefeatures))
  df_traj<-object@meta.data[test_idx,]$Trajectory
  y_numeric <- as.numeric(df_traj == posClass)
  auc_baseline=pROC::auc(pROC::roc(y_numeric, predictions))[[1]]
  
  
  improvements <- foreach::foreach(gene = genes, .combine = c) %dopar% {
    library(Seurat)
    library(caret)
    library(foreach)
    library(doParallel)
    library(e1071)
    library(pROC)
    
    # Add gene expression data to features
    features <- c(gene, basefeatures)
    df_traj<-object@meta.data[train_idx,"Trajectory"]
    y_numeric <- as.numeric(df_traj == posClass)
    tab_features <- FetchData(object, cells=train_idx, vars=features)
    model <- svm(x = tab_features, y = y_numeric, kernel = kernel, cross = cross)
    
    # Evaluate model performance with added gene
    df_test_traj<-object@meta.data[test_idx, "Trajectory"]
    y_numeric_test <- as.numeric(df_test_traj == posClass)
    tab_features <- FetchData(object, cells=test_idx, vars=features)
    predictions <- predict(model, newdata = tab_features)
    auc=auc(roc(y_numeric_test, predictions))[[1]]
    
    names(auc) <- gene
    auc
  }
  
  # Stop parallel backend
  parallel::stopCluster(cl)
  unregister_pp()
  
  # Combine results
  improvements <- sort(improvements, decreasing = TRUE)
  t<-as.data.frame(improvements)
  t$Gene <- rownames(t)
  colnames(t) <- c("AUC", "Gene")
  
  # calculate LFCs
  mat=t(as.matrix(FetchData(object,vars = genes,layer = "data")))
  mat <- mat[t$Gene,]
  lfcs=sapply(1:nrow(mat),function(i) (log1p(mean(expm1(mat[i,object$Trajectory==compareGroups[2]])))-log1p(mean(expm1(mat[i,object$Trajectory==compareGroups[1]]))))/log(2))
  t$LFC <- lfcs
  t <- t[t$LFC != 0,]
  
  #Add AUC from baseline features
  t <- rbind("baseline"=c(auc_baseline, "baseline", 0), t)
  t$AUC <- as.numeric(t$AUC)
  t$LFC <- as.numeric(t$LFC)
  t
}

#' Heterogeneity-seq: Classifying cellular response by gene expression values including causal inference by doubleML
#'
#' Classifying the cellular response of control cells using single gene expression (+ informative features) to identify features with the strongest predictive capabilities and applying causal inference by a doubleML approach.
#' @param object Seurat object
#' @param trajectories Matrix of cell-cell trajectories. Columns represent time points, rows represent trajectories of connected cells over time points.
#' @param score.group A named vector of response groups. Names represent cells, the values represent the score groups. If no score.group is set, use score.name and quantiles parameters must be set to define score groups.
#' @param score.name The name of a numeric Seurat meta data column, which will be used to calculate score groups. Only used if no score.group is given.
#' @param quantiles Thresholds of the score.name meta data to define 3 response groups. Low, Middle, High.
#' @param compareGroups Which score groups to test. Default: Low vs. High
#' @param posClass Define the positive Class for classification.
#' @param basefeatures Additional informative features to include in the classification. Must be meta data available in the Seurat object.
#' @param background A set of genes that will be considered as potential confounding factors in the doubleML analysis. Must contain all genes set in the genes parameter. By default, all genes are used.
#' @param genes Vector of genes to test.
#' @param assay The name of the Seurat assay to perform Heterogeneity-seq on. If NULL, the default assay will be used.
#' @param split Set a training-test data split. Must be in [0,1]
#' @param kernel The kernel for the SVM. linear, polynomial, radial or sigmoid. Default: radial. 
#' @param cross Number of cross-validations.
#' @param num_cores The number of cores used in parallel processing.
#' @return Table of log2FC and AUC values for each gene and an additional AUC value for the baseline features.
#' @examples 
#' # Using meta data in the Seurat object called "score" to automatically define response groups
#' HetseqDoubleML(data, trajectories, score.name = "score")
#' 
#' # Use a manually defined grouping by handing a named vector to score.group.
#' # data$infectionState contains the infection state of all cells and all time points. 
#' # Get the infectionState of all cells from the last time point in the trajectories.
#' groups <- data$infectionState[trajectories[,ncol(trajectories)],]
#' # Adapt compareGroups to fit the low and high group names of the manual definition.
#' HetseqDoubleML(data, trajectories, score.group = groups, compareGroups = c("Bystander", "High"))
#' @export
HetseqDoubleML <- function(object, trajectories, score.group = NULL, score.name=NULL,  quantiles = c(0.25,0.75), compareGroups = c("Low", "High"), posClass=NULL, basefeatures=NULL, genes=NULL, background=NULL, assay=NULL, split=NULL, kernel="radial", cross=10, num_cores=1){
  #library(parallel)
  #library(doParallel)
  #library(e1071)
  #library(pROC)
  #library(DoubleML)
  #library(mlr3)
  #library(mlr3learners)
  unregister_pp()
  
  if(!is.null(split) && (split > 1 | split < 0)) stop("split must be either NULL or in [0,1]")
  
  if(is.null(score.group)){
    if(is.null(score.name)){
      stop("Either provide a score.group or a score.name (Seurat metadata) on which to calculate the score.group")
    }
    br1=quantile(object[[score.name]][trajectories[,length(colnames(trajectories))],], probs = quantiles)[1]
    br2=quantile(object[[score.name]][trajectories[,length(colnames(trajectories))],], probs = quantiles)[2]
    score.group = cut(object[[score.name]][trajectories[,length(colnames(trajectories))],],breaks=c(-10,br1,br2,10),labels=c("Low","Middle","High"))
    object$Trajectory = "none"
    for(x in c(1:length(colnames(trajectories)))){
      object$Trajectory[trajectories[,x]] = as.character(score.group)
    }
    object$Trajectory = factor(object$Trajectory,levels=c("none",levels(score.group)))
  } else {
    object$Trajectory = "none"
    for(x in c(1:length(colnames(trajectories)))){
      object$Trajectory[trajectories[,x]] = as.character(score.group)
    }
    object$Trajectory <- as.factor(object$Trajectory)
  }
  object <- subset(object, cells = trajectories[,1], subset = Trajectory %in% compareGroups)
  object$Trajectory <- droplevels(object$Trajectory)
  print(table(object$Trajectory))
  
  
  if(is.null(basefeatures)){
    basefeatures <- c("uninformative")
    object$uninformative <- 1
  }
  
  if(is.null(split)){
    train_idx <- Cells(object)
    test_idx <-  Cells(object)
  } else {
    train_idx <- sample(1:ncol(object), split * ncol(object))
    test_idx <- setdiff(c(1:ncol(object)), train_idx)
  }
  
  # Register parallel backend
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  
  if(is.null(assay)){
    assay=DefaultAssay(object)
  }
  DefaultAssay(object) <- assay
  if(is.null(background)) background <- rownames(object)
  if(is.null(genes)) genes <- rownames(object)
  if(!all(genes %in% background)){
    stop("Error: Every element in genes must also be in background")
  }
  df_traj<-object@meta.data[train_idx,"Trajectory"]
  if(is.null(posClass)) posClass = names(table(df_traj))[2]
  
  X = data.frame(FetchData(object,cells=train_idx, vars=c(background, basefeatures)))
  
  renamed_genes <- colnames(data.frame(FetchData(object,cells=train_idx, vars=genes)))
  
  improvements <- foreach::foreach(gene = renamed_genes, .combine = rbind) %dopar% {
    library(Seurat)
    library(caret)
    library(foreach)
    library(doParallel)
    library(e1071)
    library(pROC)
    library(DoubleML)
    library(mlr3)
    library(mlr3learners)
    # Add gene expression data to features
    df_traj<-object@meta.data[train_idx,"Trajectory"]
    y_numeric <- as.numeric(df_traj == posClass)
    X$outcome <- y_numeric
    
    dml <- DoubleMLData$new(data = X, y_col = "outcome", d_cols = gene)
    learner_lasso <- lrn("regr.cv_glmnet")
    learner_gam <- lrn("regr.ranger")
    learner_svm = lrn("regr.svm")
    
    # Initialize the DoubleMLPLR model (Partial Linear Regression)
    dml_plr <- DoubleMLPLR$new(dml, ml_g=learner_gam, ml_m = learner_gam, n_folds = cross)
    dml_plr$fit()
    
    est=dml_plr$coef[gene]
    pval=dml_plr$pval
    
    names(est) <- "Estimate"
    names(pval) <- "pval"
    names(gene) <- "Gene"
    c(est,pval,gene)
  }
  print("Iterations done")
  
  # Stop parallel backend
  parallel::stopCluster(cl)
  unregister_pp()
  
  # Combine results
  t<-as.data.frame(improvements)
  t$oldGeneNames <- genes

  t$Estimate <- as.numeric(t$Estimate)
  t$pval <- as.numeric(t$pval)
  t$p.adj <- p.adjust(t$pval, method="BH")
  t <- t[order(t$Estimate, decreasing = TRUE),]
  
  t
}