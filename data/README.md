# Data

RNA and Protein expression data collected from  [Diversity Outbred](https://www.jax.org/strain/009376) mice and distributed as `.RData` files. 

 1) Kidney: [DO188b_kidney_noprobs.RData](https://www.dropbox.com/s/67mhez9xqfmsjww/DO188b_kidney_noprobs.RData?dl=0)
 
 2) Heart: [DO189_heart_v2_noprobs.RData](https://www.dropbox.com/s/bv6dvyyk37kv73l/DO189_heart_v2_noprobs.RData?dl=0)

The objects in .RData files are as follows:

- `expr.mrna` / `expr.protein` = a matrix with rankZ-tranformed normalized expression values (rows = samples, columns = genes)
- `annot.mrna` / `annot.protein` = an annotation of genes, each row is one gene (Ensembl id, MGI symbol, genome position, nearest marker, biotype)
- `annot.samples` = an annotation of samples, each row is one sample
- `N` = a list, elements are a number of genes, a number of proteins, a number of gene-protein pairs and a number of gene-protein pairs with no NA values (some protein expressions are NA)

Genotype and haplotype reconstructions are available upon request.
