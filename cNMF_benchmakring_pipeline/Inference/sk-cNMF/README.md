To run sk-cNMF, run 

pip install git+https://github.com/EngreitzLab/sk_cNMF.git


| field | meaning |
|---------------|---------------------------|
| counts_fn | Path to input counts matrix. |
| output_directory | Output directory path |
| run_name | Name for the analysis run |
| k | Number of components to factorize |
| seed | A number to set seed for reproducibility |
| n_iter | Number of iterations for factorization. |
| loss | Loss function choices: "frobenius", "kullback-leibler" |
| init | Number of NMF runs (recommend 100 for actual data analysis, and 10 for testing the pipeline) |
| algo | Solver choices: "cd", "mu"<br>- 'cd' is a Coordinate Descent solver that uses Fast Hierarchical Alternating Least Squares (Fast HALS).<br>- 'mu' is a Multiplicative Update solver. |
| total_workers | Number of processes to run in parallel |
| numhvgenes | Number of highly variable genes to use |
| max_NMF_iter | Maximum number of iterations for NMF solver |
| sel_thresh | Threshold for filtering NMF runs during consensus step |

