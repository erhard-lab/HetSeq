unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

distmat=function(prev.t,next.t, prevAssay, nextAssay, gene_subset = NULL) {
  if(is.null(gene_subset)){
    gene_subset = c(rownames(prev.t),rownames(next.t))
  }
  DefaultAssay(prev.t)=prevAssay
  DefaultAssay(next.t)=nextAssay
  prev.t=FindVariableFeatures(NormalizeData(prev.t))
  next.t=FindVariableFeatures(NormalizeData(next.t))
  features <- intersect(SelectIntegrationFeatures(list(prev.t,next.t)), gene_subset)
  anchors <- FindIntegrationAnchors(list(prev.t,next.t), anchor.features = features)
  combined <- IntegrateData(anchors)
  
  DefaultAssay(combined) <- "integrated"
  combined <- ScaleData(combined, verbose = FALSE)
  combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
  combined <- RunUMAP(combined, reduction = "pca", dims = 1:15)
  combined <- FindNeighbors(combined, reduction = "pca", dims = 1:15)
  #m=combined@assays$integrated@data
  m=t(combined@reductions$pca@cell.embeddings[,1:15])
  t1=as.matrix(m[,1:ncol(prev.t)])
  t2=as.matrix(m[,(1+ncol(prev.t)):ncol(m)])
  
  sqeuc=function(a,b) sum((a-b)^2)
  allcomb = apply(expand.grid(1:ncol(t1),1:ncol(t2)),1,function(ii) sqeuc(t1[,ii[1]],t2[,ii[2]]))
  re=matrix(allcomb,nrow=ncol(t1))
  rownames(re)=colnames(t1)
  colnames(re)=colnames(t2)
  re
}

compute.hetseq = function(mat,A,B) {
  ps=sapply(1:nrow(mat),function(i) wilcox.test(mat[i,A],mat[i,B])$p.value)
  eff=sapply(1:nrow(mat),function(i) (log1p(mean(expm1(mat[i,A])))-log1p(mean(expm1(mat[i,B]))))/log(2))
  expr=sapply(1:nrow(mat),function(i) log1p(mean(expm1(mat[i,A|B]))))
  re = data.frame(Expr=expr,LFC=eff,P=ps, Q = p.adjust(ps,method="BH"))
  rownames(re) = rownames(mat)
  re
}

hetseq_classify<-function(object, target, basefeatures, genes=NULL, assay=NULL, split, kernel, posClass=NULL, cross=10){
  library(parallel)
  library(doParallel)
  library(e1071)
  library(pROC)
  
  unregister_dopar()
  if(is.null(split)){
    train_idx <- Cells(object)
    test_idx <-  Cells(object)
  } else {
    train_idx <- sample(1:ncol(object), split * ncol(object))
    test_idx <- setdiff(c(1:ncol(object)), train_idx)
  }
  num_cores <- 30
  
  # Register parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  if(is.null(assay)){
    assay=DefaultAssay(object)
  }
  DefaultAssay(object) <- assay
  if(is.null(genes)) genes <- rownames(object)
  df_traj<-object@meta.data[train_idx,target]
  if(is.null(posClass)) posClass = names(table(df_traj))[2]
  y_numeric <- as.numeric(df_traj == posClass)
  model <- svm(x = FetchData(object,cells=train_idx, vars=basefeatures), y = y_numeric, kernel = kernel, cross = cross)
  predictions <- predict(model, newdata = FetchData(object, cells=test_idx, vars=basefeatures))
  df_traj<-object@meta.data[test_idx,]$Trajectory
  y_numeric <- as.numeric(df_traj == posClass)
  auc_baseline=auc(roc(y_numeric, predictions))[[1]]
  
  
  improvements <- foreach(gene = genes, .combine = c) %dopar% {
    library(Seurat)
    library(caret)
    library(foreach)
    library(doParallel)
    library(e1071)
    library(pROC)
    
    # Add gene expression data to features
    features <- c(gene, basefeatures)
    df_traj<-object@meta.data[train_idx,target]
    y_numeric <- as.numeric(df_traj == posClass)
    tab_features <- FetchData(object, cells=train_idx, vars=features)
    model <- svm(x = tab_features, y = y_numeric, kernel = kernel, cross = cross)
    # Evaluate model performance with added gene
    df_test_traj<-object@meta.data[test_idx, target]
    y_numeric_test <- as.numeric(df_test_traj == posClass)
    tab_features <- FetchData(object, cells=test_idx, vars=features)
    predictions <- predict(model, newdata = tab_features)
    auc=auc(roc(y_numeric_test, predictions))[[1]]
    #feature_importances <- t(model$coefs) %*% model$SV
    #auc <- paste0(auc,"_",feature_importances[1])
    names(auc) <- gene
    auc
  }
  
  # Stop parallel backend
  stopCluster(cl)
  unregister_dopar()
  
  # Combine results
  improvements <- sort(improvements, decreasing = TRUE)
  improvements <- c("baseline"=auc_baseline, improvements)
  t<-as.data.frame(improvements)
  t$Gene <- rownames(t)
  colnames(t) <- c("AUC", "Gene")
  t
}


HetseqClassify<-function(object, trajectories, score.group = NULL, score.name=NULL,  quantiles = c(0.25,0.75), compareGroups = c("Low", "High"), basefeatures=NULL, genes=NULL, assay=NULL, split, kernel, posClass=NULL, cross=10){
  library(parallel)
  library(doParallel)
  library(e1071)
  library(pROC)
  
  unregister_dopar()
  
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
  num_cores <- 30
  
  # Register parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  if(is.null(assay)){
    assay=DefaultAssay(object)
  }
  DefaultAssay(object) <- assay
  if(is.null(genes)) genes <- rownames(object)
  df_traj<-object@meta.data[train_idx,"Trajectory"]
  if(is.null(posClass)) posClass = names(table(df_traj))[2]
  y_numeric <- as.numeric(df_traj == posClass)
  model <- svm(x = FetchData(object,cells=train_idx, vars=basefeatures), y = y_numeric, kernel = kernel, cross = cross)
  predictions <- predict(model, newdata = FetchData(object, cells=test_idx, vars=basefeatures))
  df_traj<-object@meta.data[test_idx,]$Trajectory
  y_numeric <- as.numeric(df_traj == posClass)
  auc_baseline=auc(roc(y_numeric, predictions))[[1]]
  
  
  improvements <- foreach(gene = genes, .combine = c) %dopar% {
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
    #feature_importances <- t(model$coefs) %*% model$SV
    #auc <- paste0(auc,"_",feature_importances[1])
    names(auc) <- gene
    auc
  }
  
  # Stop parallel backend
  stopCluster(cl)
  unregister_dopar()
  
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

library(lpSolve)
library(ggplot2)

prune=function(D.list,top.n = 10) {
  if(length(D.list)==1){
    D=D.list[[1]]
    re=lapply(D.list,function(D) {
      prune = matrix(0,nrow=nrow(D),ncol=ncol(D))
      for (i in 1:nrow(D)) prune[i,rank(D[i,])>top.n] = prune[i,rank(D[i,])>top.n]+1
      for (i in 1:ncol(D)) prune[rank(D[,i])>top.n,i] = prune[rank(D[,i])>top.n,i]+1
      D[prune==2]=Inf
      D
    })
    keep = apply(re[[1]],1,function(v) sum(is.finite(v))>0)
    re[[1]] = re[[1]][keep,]
    re
  } else {
    re=lapply(D.list,function(D) {
      prune = matrix(0,nrow=nrow(D),ncol=ncol(D))
      for (i in 1:nrow(D)) prune[i,rank(D[i,])>top.n] = prune[i,rank(D[i,])>top.n]+1
      for (i in 1:ncol(D)) prune[rank(D[,i])>top.n,i] = prune[rank(D[,i])>top.n,i]+1
      D[prune==2]=Inf
      D
    })
    keep = apply(re[[1]],1,function(v) sum(is.finite(v))>0)
    re[[1]] = re[[1]][keep,]
    for (i in 2:length(re)) {
      keep = apply(re[[i-1]],2,function(v) sum(is.finite(v))>0) & apply(re[[i]],1,function(v) sum(is.finite(v))>0)
      re[[i-1]] = re[[i-1]][,keep]
      re[[i]] = re[[i]][keep,]
    }
    keep = apply(re[[length(re)]],2,function(v) sum(is.finite(v))>0)
    re[[length(re)]] = re[[length(re)]][,keep]
    re
  }
}

mincostflow = function(D.list, verbose=TRUE) {
  
  if (length(D.list)>1) for (i in 2:length(D.list)) if (ncol(D.list[[i-1]])!=nrow(D.list[[i]]) || !all(colnames(D.list[[i-1]])==rownames(D.list[[i]]))) stop("Dimensions do not match!")
  groups = c(nrow(D.list[[1]]),sapply(D.list,ncol))
  dist.edges = sapply(D.list,function(D) sum(is.finite(D)))
  num.edges = sum(dist.edges) + sum(groups)
  
  if (verbose) cat(sprintf("Computing max flow...\n"))
  
  ren=function(D,i) {rownames(D)=paste(i,rownames(D));colnames(D)=paste(i,colnames(D));D}
  g=do.call("rbind",lapply(1:length(D.list),function(di) reshape2::melt(ren(D.list[[di]],di))))
  g=g[is.finite(g$value),]
  g$value=1
  colnames(g) <- c("from", "to", "capacity")
  g = rbind(g,data.frame(from="source",to=paste(1,rownames(D.list[[1]])),capacity=1))
  if (length(D.list)>1) for (i in 2:length(D.list)) g = rbind(g,data.frame(from=paste(i-1,rownames(D.list[[i]])),to=paste(i,rownames(D.list[[i]])),capacity=1))
  g = rbind(g,data.frame(from=paste(length(D.list),colnames(D.list[[length(D.list)]])),to="sink",capacity=1))
  g <- igraph::graph_from_data_frame(g)
  total = igraph::max_flow(g, source = "source", target = "sink")$value
  
  if (verbose) cat(sprintf("max flow=%d, bottleneck=%d\n",total,min(groups)))
  
  #total = min(groups)
  num.constraints = sum(groups) + 2 + sum(sapply(D.list,function(D) nrow(D)+ncol(D)))
  # the first edges are the finite dist edges of D1, then D2, ... , Dn
  
  if (verbose) cat(sprintf("Allocate %dx%d matrix...\n",num.constraints,num.edges))
  lhs = matrix(NA,ncol=num.edges,nrow = num.constraints)
  
  
  # construct contraints enforcing capacity of 1 for each cell
  capacities = sum(groups)
  if (verbose) cat(sprintf("Constructing capacity constraints for %d cells...\n",capacities))
  lhs[1:capacities,] = cbind(matrix(0,nrow=capacities,ncol=sum(dist.edges)),diag(capacities))
  dir = rep("<=",capacities)
  rhs = rep(1,capacities)
  
  cpad = rep(0,sum(dist.edges))
  # add source constraints
  lhs[capacities+1,] = c(cpad,rep(1,groups[1]),rep(0,capacities-groups[1]))
  dir = c(dir,"==")
  rhs = c(rhs, total)
  # add sink constraints
  lhs[capacities+2,] = c(cpad,rep(0,capacities-groups[length(groups)]),rep(1,groups[length(groups)]))
  dir = c(dir,"==")
  rhs = c(rhs, total)
  
  # add flow constraints
  if (verbose) cat(sprintf("Constructing flow constraints and costs for %d edges...\n",sum(dist.edges)))
  cpad = rep(0,sum(groups))
  bpad = c()
  cpos = 1
  ncstr = capacities+3
  for (D in D.list) {
    # add contraints for left group
    for (r in 1:nrow(D)) {
      rr=matrix(0,nrow=nrow(D),ncol=ncol(D))
      rr[r,] = 1
      rr=as.vector(rr)[is.finite(D)]
      row = c(bpad,rr)
      row = c(row,rep(0,sum(dist.edges)-length(row)))
      rest = rep(0,num.edges-length(row))
      rest[cpos] = -1
      lhs[ncstr,] = c(row,rest)
      ncstr=ncstr+1
      cpos=cpos+1
    }
    # add contraints for right group
    for (c in 1:ncol(D)) {
      cc=matrix(0,nrow=nrow(D),ncol=ncol(D))
      cc[,c] = 1
      cc=as.vector(cc)[is.finite(D)]
      row = c(bpad,cc)
      row = c(row,rep(0,sum(dist.edges)-length(row)))
      rest = rep(0,num.edges-length(row))
      rest[cpos] = -1
      lhs[ncstr,] = c(row,rest)
      ncstr=ncstr+1
      cpos=cpos+1
    }
    cpos = cpos-ncol(D)
    bpad = c(bpad,rep(0,sum(is.finite(D))))
  }
  dir = c(dir,rep("==",nrow(lhs)-length(dir)))
  rhs = c(rhs,rep(0,nrow(lhs)-length(rhs)))
  
  obj = c(do.call("c",lapply(D.list,function(D) as.vector(D))), rep(0,sum(groups)))
  obj = obj[is.finite(obj)]
  #print(rbind(cbind(lhs,rhs),c(obj,NA)))
  
  if (verbose) cat(sprintf("Optimizing %d variables with %d constraints...\n",ncol(lhs),nrow(lhs)))
  solution <- lpSolve::lp(
    direction = 'min',
    objective.in = obj,
    const.mat = lhs,
    const.dir = dir,
    const.rhs = rhs)
  
  mergeorder = function(a,b) {
    re=merge(a,b,by='C')
    re=re[,c(2:ncol(a),1,(ncol(a)+1):ncol(re))]
    re
  }
  if (verbose) cat(sprintf("Constructing %d trajectories...\n",total))
  df=NULL
  sol = solution$solution
  pos = 1
  for (di in 1:length(D.list)) {
    D=D.list[[di]]
    edges = dist.edges[di]
    mm = rep(0,nrow(D)*ncol(D))
    mm[as.vector(is.finite(D))]=sol[pos:(pos+edges-1)]
    pos = pos+edges
    
    mm=matrix(mm,ncol=ncol(D))
    rownames(mm) = rownames(D)
    colnames(mm) = colnames(D)
    mm = reshape2::melt(mm)
    mm = mm[mm$value!=0,]
    dfn = setNames(data.frame(as.character(mm$Var1),as.character(mm$Var2)),c("C",""))
    df = if (is.null(df)) dfn else mergeorder(df,dfn)
    names(df)=c(rep("",ncol(df)-1),"C")
  }
  names(df)=paste0("G",1:(length(D.list)+1))
  df
}


randomizeTrajectories <- function(trajectories, percentage){
  nrow = dim(trajectories)[1]
  ncol = dim(trajectories)[2]
  count = round(nrow*percentage)
  randRows = sample.int(nrow, count)
  for(i in 1:count){
    if(i%%2==0){
      next
    } else if(i==count){
      next
    }
    
    c1 = trajectories[randRows[i], ncol]
    trajectories[randRows[i], ncol] <- trajectories[randRows[i+1], ncol]
    trajectories[randRows[i+1], ncol] <- c1
  }
  trajectories
}

percentageCorrectTrajectories <- function(traj, noisyTraj){
  rownames(noisyTraj) <- noisyTraj[,1]
  sum(noisyTraj[traj[,1],dim(traj)[2]] == traj[,dim(traj)[2]])/dim(traj)[1]
}

HetSeq_causalTree <- function(object, trajectories, score.group = NULL, score.name=NULL,  quantiles = c(0.25,0.75), compareGroups = c("Low", "High"), basefeatures=NULL, genes=NULL, assay=NULL, returnforest=FALSE, split, kernel, posClass=NULL, cross=10, num_cores=1){
  library(parallel)
  library(doParallel)
  library(e1071)
  library(pROC)
  library(grf)
  
  unregister_dopar()
  
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
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  if(is.null(assay)){
    assay=DefaultAssay(object)
  }
  DefaultAssay(object) <- assay
  if(is.null(genes)) genes <- rownames(object)
  df_traj<-object@meta.data[train_idx,"Trajectory"]
  if(is.null(posClass)) posClass = names(table(df_traj))[2]
  
  
  X <- as.matrix(FetchData(object, cells = train_idx, vars = genes))
  Y <- as.numeric(object@meta.data[train_idx, "Trajectory"] == posClass)
  if(!returnforest){
    improvements <- foreach(gene = genes, .combine = rbind, .packages = c("grf")) %dopar% {
      W <- X[, gene]
      X_other <- X[, setdiff(colnames(X), gene)]
      causal_forest <- causal_forest(X_other, Y, W)
      
      # Estimate the average treatment effect
      ate <- average_treatment_effect(causal_forest, target.sample = "all")
      out <- c(gene, ate)
      names(out) <- c("Gene", "ATE" ,"stderr")
      out
    }
    # Stop parallel backend
    stopCluster(cl)
    unregister_dopar()
    
    t<-as.data.frame(improvements)
    t$ATE <- as.numeric(t$ATE)
    t$stderr <- as.numeric(t$stderr)
    t <- t[order(t$ATE, decreasing = TRUE),]
    t
  } else {
    improvements <- foreach(gene = genes, .packages = c("grf")) %dopar% {
      W <- X[, gene]
      X_other <- X[, setdiff(colnames(X), gene)]
      causal_forest <- causal_forest(X_other, Y, W)
      causal_forest
    }
  }
}
HetSeq_doubleML <- function(object, trajectories, score.group = NULL, score.name=NULL,  quantiles = c(0.25,0.75), compareGroups = c("Low", "High"), basefeatures=NULL, genes=NULL, background=NULL, assay=NULL, split, kernel, posClass=NULL, cross=10, num_cores=1){
  library(parallel)
  library(doParallel)
  library(e1071)
  library(pROC)
  library(DoubleML)
  library(mlr3)
  library(mlr3learners)
  unregister_dopar()
  
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
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
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
  # y_numeric <- as.numeric(df_traj == posClass)
  # model <- svm(x = FetchData(object,cells=train_idx, vars=basefeatures), y = y_numeric, kernel = kernel, cross = cross)
  # predictions <- predict(model, newdata = FetchData(object, cells=test_idx, vars=basefeatures))
  # df_traj<-object@meta.data[test_idx,]$Trajectory
  # y_numeric <- as.numeric(df_traj == posClass)
  # auc_baseline=auc(roc(y_numeric, predictions))[[1]]
  
  X = data.frame(FetchData(object,cells=train_idx, vars=c(background, basefeatures)))


  
  renamed_genes <- colnames(data.frame(FetchData(object,cells=train_idx, vars=genes)))

  improvements <- foreach(gene = renamed_genes, .combine = rbind) %dopar% {
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

    #model <- svm(x = FetchData(object,cells=train_idx, vars=basefeatures), y = y_numeric, kernel = kernel, cross = cross)
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
  stopCluster(cl)
  unregister_dopar()
  
  # Combine results
  t<-as.data.frame(improvements)
  t$oldGeneNames <- genes
  #t$Gene <- rownames(t)
  
  #colnames(t) <- c("Estimate", "Gene", "geneRename")
  t$Estimate <- as.numeric(t$Estimate)
  t$pval <- as.numeric(t$pval)
  t <- t[order(t$Estimate, decreasing = TRUE),]
  
  t
}
