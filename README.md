# Assessing difference between different NMF methods 


Gene programs inferred from single-cell genomic data (scRNASeq., scATACseq., multi-omics and Perturb-seq.) are useful in discovering contextual biological mechanisms. These programs can be viewed as data-driven hypotheses of gene interactions. We aim to implement a flexible framework to evaluate the plausibility of programs inferred by computational methods. 

We break down the accessment into vanilla NMF methods and cNMF methods. The former focuses on basic metric benchmarking and stability evaluation. The latter is broken down into themes such as goodness if fit (ability to explain the data), co-regulation, mechanistic interactions etc. Under each theme, multiple evaluation tasks are conceptualised and implemented using appropriate statistical tests.


## cNMF benchmarking

Understand the difference between different cNMF methods. Both Jupyter Notebook version and Slurm version are avalible to run inference, evaluation, and plotting seperately and all-together. 

### Inference
Versions of cNMF:

1. CPU powered cNMF
  * NMF inference using: [sklearn.decomposition.non_negative_factorization](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.non_negative_factorization.html)
  * Solver choice: multiplicative update, coordinate descent
  * Mode: batch

2.  GPU powered cNMF
  * NMF inference using: [NMF-Torch](https://github.com/lilab-bcb/nmf-torch)
  * Solver choice: multiplicative update, hierarchical alternative least square
  * Mode: batch, mini-batch (online) 


### Evaluation
#### Basic metric evaluation:
  * Speed
  * Memory usage

#### Statisitical evaluation:
  * Reconstructive error
  * Stability with silhouette score
  * Euclidean distance clustermap 
  * Correlation clustermap
  * Top 300 gene overlap clustermap

#### Biologcoal metrics: 

| Criterion    | Implementation | External resource | Interpretation | Caveats |
| -------- | ------- | -------- | ------- | ------- |
| Goodness of fit  | [Explained variance](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.explained_variance_score.html) per program | None | A program explaining more variance in the data might represent dominant biological variation. | Technical variation might be the highest source of variance (e.g. batch effects). |
| Variation across category levels | [Kruskall-Wallis non-parametric ANOVA](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) + Dunn's posthoc test | None | If program scores are variable between batch levels then the component likely is modelling technical noise. Alternatively, if program scores are variable between a biological category like cell-type or condition then the program is likely modelling a biological process specific to the category. | If batches are confounded with biological conditions, then the relative contribution of technical and biological variation cannot be decomposed. |
| Gene-set enrichment | [GSEA](https://gseapy.readthedocs.io/en/latest/introduction.html) using program x feature scores | MsigDB, Enrichr | If a program is significantly associated with a gene-set then it could explain the biological process the program represents | |
| Motif enrichment | Pearson correlation of motif counts per gene (promoter or enchancer) and program x gene scores | HOCOMOCO v12 | If genes with high contributions to a program are also enriched with same enhancer/promoter motifs they could be co-regulated | A biological pathway could involve genes with different regulation but still contribute to a common function | 
| Trait enrichment | [Fisher's exact test](https://en.wikipedia.org/wiki/Fisher%27s_exact_test) | OpenTargets database | If a program is significantly associated with a trait then it could explain the biological process the program represents | |
| Perturbation sensitivity | Mann-Whitney U test of program scores between perturbed cells and non-targeted/reference cells | Perturbation data | Cell x programs score distribution shifts greater than expected due to the direct effect of perturbation on genes in the program could indicate hierarchical relationships b/w genes in the program | Expression of genes upstream of the perturbed gene are unlikely to be affected | 


### Plotting
Visual representations of the evaluation metrics. 
