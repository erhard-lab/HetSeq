# Heterogeneity-Seq

Heterogeneity-seq is a concept recently developed in our preprint on [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.10.28.620481v1). In short, tracing back perturbed cells to the most likely pre-perturbation control cells enables exploiting intercellular heterogeneity to identify factors that modulate cellular responses to perturbations. The HetSeq functions can be applied to calculate metabolic labeling based trajectories in a scSLAM-seq time course data set and subsequently identify pathway modulators with three different approaches. In general, these HetSeq approaches can also be applied in other types of data sets and trajectory calculation methods.

The following tutorial considers the analysis of a scSLAM-seq time course data set using the Seurat package.


## Installation

	install.packages("devtools")
	devtools::install_github("https://github.com/erhard-lab/HetSeq/")
	library(HetSeq)

Installing via install_github will also ensure that the following packages will be installed alongside the HetSeq package:

```
caret,doParallel,DoubleML,e1071,foreach,ggplot2,grf,igraph,lpSolve,mlr3,mlr3learners,parallel,pROC,reshape2,Seurat
```

## Calculating previous RNA profiles (optional)


When creating trajectories based on labeling information, we can connect cells across time points by calculating the distances between expression profiles of cells. If labeling information is available (e.g. in a scSLAM-seq experiment), we can calculate distances between the pre-existing (old/unlabeled) RNA of cells from a time point with the total RNA profiles of cells from a previous time point, as this old RNA profile represents RNA from before the onset of labeling and should therefore be similar to the overall expression of cells from the previous time point. However, we can even go a step further and take half-lives and RNA degradation into account to calculate a more accurate representation of a cells previous RNA profile. For this, we can make use of the [grandR](https://github.com/erhard-lab/grandR) package:


	library(Seurat)
	library(grandR)
	data #Seurat object with oldRNA and newRNA assays
	data_grandR #grandR object of the data set
	labeling_dur = 2

	SetParallel(cores = 16) 
	data_grandR<-ComputeSteadyStateHalfLives(data_grandR, as.analysis = F)
	HL <- GetMatrix(data_grandR, "HL")
	prevRNA = (data@assays$oldRNA$counts * exp(labeling_dur * log(2)/pmin(pmax(HL, 0.25), 24)))[rownames(data),colnames(data)]
	data@assays$prevRNA <- CreateAssayObject(counts = prevRNA, key = "prevrna_")


This would generate previous RNA profiles for a labeling duration of 2 hours. For a complete example on how to load your data set as a grandR object and how to calculate half-lives, check out the [grandR vignettes](https://grandr.erhard-lab.de/articles/getting-started.html).


## Building cellular trajectories

To connect cells treated cells back to pre-treated control cells, we need to build cellular trajectories. This segment will demonstrate how to generate these trajectories based on metabolic labeling information using the HetSeq package. Cell-cell trajectories built any other way can also be utilized, but need to be in the following matrix format, where every row represents one trajectory from the last time point back to the first (control sample) time point:


	| Timepoint1 | Timepoint2 | Timepoint3 | Timepoint4 | Timepoint5 |
	|------------|------------|------------|------------|------------|
	| Cell1      | Cell2      | Cell3      | Cell4      | Cell5      |
	| Cell6      | Cell7      | Cell8      | Cell9      | Cell10     |
	| Cell11     | Cell12     | Cell13     | Cell14     | Cell15     |
	| Cell16     | Cell17     | Cell18     | Cell19     | Cell20     |
	| Cell21     | Cell22     | Cell23     | Cell24     | Cell25     |


Consider a Seurat object, where the time course information is saved in the treatmentTime meta data and we have total RNA (RNA) and previous RNA (prevRNA) assays:

	treatment.list <- SplitObject(subset(data, subset = treatmentTime %in% c("0h", "2h", "4h")), split.by = "treatmentTime")
	D.list=list(
	  distmat(treatment.list[["0h"]],treatment.list[["2h"]], "RNA", "prevRNA"),
	  distmat(treatment.list[["2h"]],treatment.list[["4h"]], "RNA", "prevRNA")
	)
	D.list = prune(D.list)
	trajectories = mincostflow(D.list)

First, we create distance matrices between cells of two time points, using the total RNA profiles of the previous time point and the previous RNA profiles of the succeeding time point. Alternatively, if previous RNA profiles were not calculated, one could also use old RNA profiles instead. The list of distance matrices is then pruned to retain only the top 10 connections for every cell. From this, in a min-cost max-flow approach, the overall best set of trajectories is generated as a matrix (see above) and can now be used for the HetSeq functions.

## Predicting cellular outcomes using an SVM


## Identifying causal pathway modulators with a Causal Inference approach

