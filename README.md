# marbles
MultivARiate Bayesian Layered variablE Selection

R package for "[Mohammed, S.](shariq-mohammed.github.io), Kurtek, S., Bharath, K., Rao, A. and Baladandayuthapani, V., 2021. Tumor Radiogenomics with Bayesian Layered Variable Selection." <!--[https://arxiv.org/abs/2106.10941](https://arxiv.org/abs/2106.10941)"-->

Code to run the model in (Mohammed et.al, 2021) on the imaging and genomic data in gliomas using the package `marbles`.

### Contents

1. Accessing low grade glioma (LGG) data from the package
2. Fit the multivariate Bayesian layered variable selection model for LGG data
3. Plot radiogenomic associations in LGG

### Accessing low grade glioma data from the package

Installing the package
```
# install the package (devtools package needed)
if(!require(devtools)) install.packages("devtools")
devtools::install_github('author/marbles')
```

Load the package
```
library(marbles)
```

#### Example

The complete description for each of the data files can be accessed in **R** or **RStudio** as follows:

```
?cancer_gene_exp_LGG
?mri_data
?pdfs
?pc_scores
?group_id
```

Data can be accessed directly using the variable names once the package is loaded in the environment. For example, the structure and dimensions of the `cancer_gene_exp_LGG` object can be checked as follows
```
# structure
str(cancer_gene_exp_LGG)

# dimensions
dim(cancer_gene_exp_LGG)
```

### Fit the multivariate Bayesian layered variable selection model for LGG data

Run the sequential estimation model for different values of $v_0$.
```
# sequence of values for v_0
v0.seq = seq(from = 0.001, by = 0.001, length.out = 40)
```

Estimate the model parameters with `pc_scores` as the multivariate response and `cancer_gene_exp_LGG` as the predictors.
```
res = seq_est_model_selection(Y = pc_scores, X = cancer_gene_exp_LGG,
                              g_id = group_id, nCores = 4, v0seq = v0.seq)
```

### Plot radiogenomic associations in LGG

Create a disc plot using MR imaging sequence names and gene names.
```
# imaging sequence names
im.seq = c('FLAIR', 'T1', 'T1Gd', 'T2')

# gene names
X.colnames = colnames(cancer_gene_exp_LGG)

disc_plot(res$seq_res$w, im.seq, X.colnames, layer = T)
```
