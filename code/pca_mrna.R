# PCA - mRNA only
# 
# Principal component analysis for mRNA/Protein 
#
# Usage:
#
# Rscript pca_mrna.R input.RData output_loadings.csv output_predict.csv
# 
# output_loadings.csv =   gene-level PCA coefficients
# output_predict.csv  = sample-level PCA coefficients

# Options -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# for kidney
# args = list("../kidney2/R/DO188b_kidney_noprobs.RData", 
#              "../kidney2/kidney_pca_loadings_mrna.csv", 
#              "../kidney2/kidney_pca_predict_mrna.csv")

# two arguments expected (one input and two output files)
stopifnot(length(args)==3)

input.file = args[[1]]
stopifnot(file.exists(input.file))
output.file1 = args[[2]] # loadings (gene-level predictions)
output.file2 = args[[3]] # predict (sample-level predictions)

load(input.file)

# PCA ---------------------------------------------------------------------

# PCA calculation
pca.mrna <- prcomp(scale(as.matrix(expr.mrna)))

# how much is explained by top10 PCAs?
print(pca.mrna$sdev[1:10]^2 / sum(pca.mrna$sdev^2))

# Output loadings of top 10 PCAs------------------------------------

loadings.mrna <- pca.mrna$rotation[,1:10]
colnames(loadings.mrna) <- paste(colnames(loadings.mrna), "mrna", sep="_")

annot.cols <- c("id",	"symbol",	"chr",	"start",	"end",
                "strand",	"biotype", "duplicated")
output <- cbind(annot.mrna[, annot.cols], loadings.mrna)
output <- subset(output, !duplicated)

write.csv(output, file=output.file1, row.names=FALSE)


# Output sample predictions for top 10 PCAs -------------------------------

predict.mrna <- predict(pca.mrna)[,1:10]
colnames(predict.mrna) <- paste(colnames(predict.mrna), "mrna", sep="_")

output2 <- cbind(annot.samples, predict.mrna)
write.csv(output2, file=output.file2, row.names=FALSE)
