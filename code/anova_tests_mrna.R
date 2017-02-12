# ANOVA - mRNA only
# 
# Calculates averages and fit several ANOVA models.
# (special mRNA version because most of mRNA do not 
# have the corresponding protein, derived from anova_tests.R)
#
# Usage:
#
# Rscript anova_tests_mrna.R input.RData output.csv

library(broom)
# library(ppcor) # library for partial correlation coefs

# Options -----------------------------------------------------------------

# args <- commandArgs(trailingOnly = TRUE)

# for kidney
args = list("../kidney2/R/DO188b_kidney_noprobs.RData", "../kidney2/kidney_anova_mrna_table.csv")
# for heart
# args = list("../heart2/DO189_heart_v2_noprobs.RData", "../heart2/heart_anova_mrna_table.csv")

# two arguments expected (input and output files)
stopifnot(length(args)==2)

input.file = args[[1]]
stopifnot(file.exists(input.file))
output.file = args[[2]]

load(input.file)


# Averages per group -------------------------------------------------------
## calculates mean expression for each group
## (for example mRNA average expression of Age==6 samples)

averages.mrna <- NULL
vars <- c("Age", "Sex", "Generation")
for (v in vars)
  for (u in sort(unique(annot.samples[,v]))) {
    sel <- annot.samples[,v] == u
    new.average <- colMeans(expr.mrna[sel,], na.rm=TRUE)
    averages.mrna <- cbind(averages.mrna, new.average)
    colnames(averages.mrna)[ncol(averages.mrna)] <- paste("average", "mRNA", v, u, sep=".")
  }



# ANOVA (p-values and normalized coefs) -----------------------------------

## test for dependence between Age/Sex and x (x is mRNA/Prot expression)
anova_tests_2 <- function(x) {
  # full model without interaction
  tmp.full <- lm(x ~ Age + Sex + Generation, data=annot.samples)
  lm.full <- tidy(tmp.full)

  pvalues <- c(subset(lm.full, term=="Age")$p.value,
               subset(lm.full, term=="Sex")$p.value)

  # partial correlation coefficient equals normalized beta * CONST
  # see http://stats.stackexchange.com/questions/76815/multiple-regression-or-partial-correlation-coefficient-and-relations-between-th
  
  pres <- !is.na(x) # must be calculated on the same observations
  tmp.age  <- lm(Age ~ Sex + Generation, data=annot.samples[pres,])
  tmp.sex  <- lm(Sex ~ Age + Generation, data=annot.samples[pres,])
  sigma.full <- sd(tmp.full$resid, na.rm = TRUE)
  sigma.age <- sd(tmp.age$resid, na.rm = TRUE)
  sigma.sex <- sd(tmp.sex$resid, na.rm = TRUE)
    
  coefs <- c(subset(lm.full, term=="Age")$estimate * sigma.age,
             subset(lm.full, term=="Sex")$estimate * sigma.sex) / sigma.full
  
  return(c(pvalues,coefs))
}

result.table <- matrix(NA, ncol(expr.mrna), 4)
colnames(result.table) <- c("p.mRNA_Age.Sex", "p.mRNA_Sex.Age",
                            "r.mRNA_Age.Sex", "r.mRNA_Sex.Age")
                            
                    
## normalize verything to mean=0, sd=1
## for easy calculation of partial correlation coefficients
for (i in 1:ncol(expr.mrna)) {
  expr.mrna[,i] <- scale(expr.mrna[,i])
}
annot.samples$Age <- scale(annot.samples$Age)  
annot.samples$Sex <- scale(as.numeric(factor(annot.samples$Sex)))

# ANOVA for genes
print("Testing for dependence between mRNA and Age/Sex...")
for (i in 1:ncol(expr.mrna)) {
  if (i %% 100 == 0) print(i)
  result.table[i, 1:4] <- anova_tests_2(expr.mrna[,i])
}



# Output ------------------------------------------------------------------

annot.cols <- c("id",	"symbol",	"chr",	"start",	"end",
                "strand",	"biotype", "duplicated")

output <- cbind(annot.mrna[, annot.cols], 
                result.table,
                round(averages.mrna,4))
# remove mRNAs that are there twice because of pairing
output <- subset(output, !duplicated)

write.csv(output, output.file, row.names=FALSE)
