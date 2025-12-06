# Assessing difference between different NMF and cNMF methods 


Gene programs inferred from single-cell genomic data (scRNASeq., scATACseq., multi-omics and Perturb-seq.) are useful in discovering contextual biological mechanisms. These programs can be viewed as data-driven hypotheses of gene interactions. We aim to implement a flexible framework to evaluate the plausibility of programs inferred by computational methods. 

We break down the accessment into vanilla NMF methods and cNMF methods. The former focuses on basic metric benchmarking and stability evaluation. The latter is broken down into themes such as goodness if fit (ability to explain the data), co-regulation, mechanistic interactions etc. Under each theme, multiple evaluation tasks are conceptualised and implemented using appropriate statistical tests.


## cNMF benchmarking

Understand the difference between different cNMF methods. Both Jupyter Notebook version and Slurm version are avalible to run inference, evaluation, and plotting seperately and all-together. 

### Inference

Versions of cNMF:

1. CPU powered cNMF
  * Individual NMF inference using: [sklearn.decomposition.non_negative_factorization](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.non_negative_factorization.html)
  * consensus NMF using [sk-cNMF](https://github.com/EngreitzLab/sk_cNMF) which is a slightly modified version from the [Orginal cNMF](https://github.com/dylkot/cNMF/tree/main) with more flexiblity to choose solver and loss function. 
  * Solver choice: multiplicative update, coordinate descent
  * Mode: batch

2.  GPU powered cNMF
  * Individual NMF inference using: [NMF-Torch](https://github.com/lilab-bcb/nmf-torch)
  * consensus NMF using: [torch-cNMF](https://github.com/ymo6/torch_based_cNMF) 
  * Solver choice: multiplicative update, hierarchical alternative least square
  * Mode: batch, mini-batch (online) 


### Evaluation
Basic metric evaluation:
  * Speed
  * Memory usage

 Statistical evaluation:
  * Reconstructive error
  * Stability with silhouette score
  * Euclidean distance clustermap 
  * Correlation clustermap
  * Top 300 gene overlap clustermap

 Biologcoal metrics: 
  * Goodness of fit 
  * Variation across category levels
  * Gene-set enrichment
  * Motif enrichment
  * Trait enrichment 
  * Perturbation sensitivity

### Plotting
K-selection plots: 
  * Stability &Error
  * GO/Genesets/Trait enrichment
  * perturbation sensitivity
  * explained variances
  * program dot plot by conditions 

Compare model plots (with same K): 
  * clustermap and boxplots for shared gene
  * GO/Genesets/Trait enrichment
  * perturbation sensitivity
  * coefficient of variance

Program QC plots
  * program UMAP
  * program violin plot
  * program loading correlations
  * top GO term plot
  * top loading genes
  * volcano plot + dot plot + waterfall plot + bar plot for regulated programs per condition of cells

Perturbed-gene plots
  * gene UMAP
  * guide UMAP
  * gene dotplot
  * gene loading correlations
  * top-loading programs
  * volcano plot + dot plot + waterfall plot + bar plot for regulated programs per condition of cells
  * Heatmap plot for regulator expression in conditions  

Excel summarization: 
  * Integrate the mdata + evaluation results information together   

## NMF benchmarking

Understand the difference between different NMF methods. 

Running replicates of different NMF methods, currently tested:
* [sklearn.decomposition.non_negative_factorization](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.non_negative_factorization.html)
* [NMF-Torch](https://github.com/lilab-bcb/nmf-torch)
* [pytorch-NMF](https://github.com/yoyolicoris/pytorch-NMF/tree/master) -> no longer considered
* [pyDNMFk](https://github.com/lanl/pyDNMFk) > no longer considered
