## _Trajectory inference across conditions: Differential expression and differential progression_
## Pre-requisites
Software:
* Basic knowledge of R syntax
* Familiarity with single-cell RNA-sequencing
* Familiarity with the SingleCellExperiment class
Background reading:
* The textbook "Orchestrating Single-Cell Analysis with Bioconductor" is a great reference for single-cell analysis using Bioconductor packages.
* [Slingshot paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4772-0)
* [tradeSeq paper](https://www.nature.com/articles/s41467-020-14766-3)
## Learning goals
* Grasp the complexity of analyzing large scRNA-seq datasets with the goal of extracting relevant biological information.
* Understand the concepts of differential progression and differential expression along a trajectory path.
* Learn how to analyze single-cell RNA-seq data using Bioconductor packages.
* Import and explore large scRNA-seq datasets.
* Understand the challenges of trajectory inference.
* Compose analysis pipeline that allows interpretation of complex scRNA-seq datasets.
* Assess the added complexity of handling multiple conditions in these dynamic systems and how it influences the analysis pipeline.
* Discuss how the analysis pipeline can incorporate this change and evaluate it.
## Dataset being used: [epithelial-to-mesenchymal transition](https://www.nature.com/articles/s41588-019-0489-5) 
**Studying EMT (epithelial to mesenchymal transition)** _One of many examples of the type of dataset being worked with in TI_
<img src="https://github.com/AlicenJoyHenning/TrajectoryInference/assets/129797527/7c9444de-fc7b-4f6c-a3fc-5c7b127dbcad" />
Epithelial cells are taken and allowed to grow. As they divide 
and differentiate they turn into mesenchymal cells, 
from which samples were taken both from the inner 
and outer regions of the plate. Inner sample was found to be more epithelial cells while outer cells 
were more mesenchymal cells. 
_But_ our experimental design of interest duplicates this, looking at both control and treatment continuums of cells. Here we are interested in comparing cell development these conditions. 
## Upstream workflow 
The upstream workflow is adjacent to what I have done before in Seurat  involving normalization with ```SCTransform```, batch integration of different conditions in Seurat with ```IntegrateData```, and identification of trajectory structure with ```Slingshot```.
* **Integration**: representing the datasets in the same space
* **Trajectory Inference**: identify the primary axis of EMT

For the upstream workflow there is conversion between Seurat objects and SingleCellExperiment objects 

## Workshop focus
### Differential Topology
But these are not the main focus of this workshop and are choices made that are debatable. This workflow begins with differential topology (looking for differences in large structure of the conditions). For example: 

Single trajectory with 2 conditions that are evenly mixed throughout, or on the right we see two conditions that have diverged and end up in different places- where not appropriate to fit one trajectory to both conditions. Here it may not be appropriate to fit a branching trajectory where an early cell, already by being from condition B, will be prediposed to one branch. As opposed to a branching trajectory, this would be two separate trajectories. 

![image](https://github.com/AlicenJoyHenning/TrajectoryInference/assets/129797527/e6405098-21ae-461e-aa4f-23ad818fd28f)

### Differential expression within & between conditions 
Differential expression in the context of continuous psuedotime with a package called **tradeSeq**, well suited to handle differential conditions problem. This allows us to discover genes that are differentially expressed. This allows us to focus on 2 main questions, 

(1) **Differential expresison within each condition**: Is a gene's expression varying significantly along psuedotime? _Are these genes exhibiting dynamic expression over psuedotime in either condition?_ 
```associationTest```

(2) **Differential expression between conditions**: Are the expression patterns along psuedotime different between the conditions? _Are genes exhibiting different patterns over psuedotime between the two conditions?_ 
```conditionTest```

This is very similar to the types of questions used to be asked of branching trajectories where we ask _is a gene following the same pattern in this branch vs another branch?_, but now we're asking in the context of conditions, _is a gene following a pattern within a condition?_

## Workshop notes 
1. In addition to the packages, you must install the package that was created to carry out this tutorial:
  ```R
  install.packages("remotes")
  remotes::install_github("kstreet13/bioc2020trajectories")
  library(bioc2020trajectories)
  ```
2. Also, and essentially, you must clone the [GitHub repo](https://github.com/kstreet13/bioc2020trajectories)
     
  ![image](https://github.com/AlicenJoyHenning/TrajectoryInference/assets/129797527/48a88c44-51b8-4a87-8029-73353550a25c)

~ATM I'm having trouble with doing the preprocessing so will just work through the workflow stuff using the preprocessed object: 
```R
data("sce", package = "bioc2020trajectories")
```
> **Differential Topology**
> Looking at neighbourhoods around each cell and asking does the composition of cells in that neighbourhood match the overall composition of the condition? Within each local neighbourhood around a given cell, we want to see if there is an even divide between the two conditions. 
> ![image](https://github.com/AlicenJoyHenning/TrajectoryInference/assets/129797527/266d6e3f-ca42-4900-82e8-8f3a3734f5b2)
>
> _kbet_: method for assessing normalizing methods is very similar 



