# Processed Data

Processed data are too big to be stored on Github but download links are provided below:

- ANOVA results for kidney [kidney_anova_table.csv](https://www.dropbox.com/s/1fvh4mzngmgvxfn/kidney_anova_table.csv?dl=0)
- ANOVA results for heart [heart_anova_table.csv](https://www.dropbox.com/s/a47z5f76feq7p3g/heart_anova_table.csv?dl=0)
- PCA results for kidney: [kidney_pca_loadings.csv](https://www.dropbox.com/s/cbdxdxs5qfw42q4/kidney_pca_loadings.csv?dl=0) and [kidney_pca_predict.csv](https://www.dropbox.com/s/f3fr33xuh50t1cy/kidney_pca_predict.csv?dl=0)
- PCA results for heart: [heart_pca_loadings.csv](https://www.dropbox.com/s/q7a4lr0nc455sa7/heart_pca_loadings.csv?dl=0) and [heart_pca_predict.csv](https://www.dropbox.com/s/s59ib1fm5bkv8ek/heart_pca_predict.csv?dl=0)

### ANOVA output:

The columns in CSV table with results are as follows:

- `id`, `gene_id`, `symbol`, `chr`, `start`, `end`, `strand`, `biotype` = gene/protein annotation.
- `p.X_Y.Z` is p-value for a test of a submodel `lm(X ~ Y + Z)` vs `lm(X ~ Z)`. Generation factor is included in all models even if not mentioned explicitely.
- `r.X_Y.Z` is a partial correlation coefficient  . Generation factor is always included in `Z` covars even if not mentioned explicitely.
- `p.X_Int.Z` is p-value for a test of an interaction submodel `lm(X ~ Age * Sex + Z)` vs `lm(X ~ Age + Sex + Z)`.
- `lod` columns
- `average` columns

Only genes with both mRNA and protein expression are used.

### PCA output:

For top 10 PCA, both for mRNA and proteins, the output tables contain loadings (=eigenvectors) and predictions (sample-level coefficients).

Currently, only genes with both mRNA and protein expression with no missing values are used.
