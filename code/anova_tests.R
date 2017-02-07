# ANOVA
# 
# Calculates averages and fit several ANOVA models.
#
# Usage:
#
# Rscript anova_tests.R input.RData output.csv

library(broom)
# library(ppcor) # library for partial correlation coefs

# Options -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# for kidney
# args = list("../kidney2/R/DO188b_kidney_noprobs.RData", "../kidney2/kidney_anova_table.csv")
# for heart
# args = list("../heart2/DO189_heart_v2_noprobs.RData", "../heart2/heart_anova_table.csv")

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
    new.average <- colMeans(expr.mrna[sel,1:N[["pairs"]]], na.rm=TRUE)
    averages.mrna <- cbind(averages.mrna, new.average)
    colnames(averages.mrna)[ncol(averages.mrna)] <- paste("average", "mRNA", v, u, sep=".")
  }

averages.protein <- NULL
vars <- c("Age", "Sex", "Generation")
for (v in vars)
  for (u in sort(unique(annot.samples[,v]))) {
    sel <- annot.samples[,v] == u
    new.average <- colMeans(expr.protein[sel,1:N[["pairs"]]], na.rm=TRUE)
    averages.protein <- cbind(averages.protein, new.average)
    colnames(averages.protein)[ncol(averages.protein)] <- paste("average", "Prot", v, u, sep=".")
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

## test for dependence between Age/Sex and x (x is mRNA/Prot expression)
## conditioned on y (y is Prot/mRNA expression)
anova_tests_3 <- function(x,y) {
  # full model without interaction
  tmp.full <- lm(x ~ Age + Sex + Generation + y, data=annot.samples)
  lm.full <- tidy(tmp.full)
  # the model without Age
  tmp.noage <- lm(x ~ Sex + Generation + y, data=annot.samples)
  lm.sex <- tidy(tmp.noage)
  # the model without Sex
  tmp.nosex <- lm(x ~ Age + Generation + y, data=annot.samples)
  lm.age <- tidy(tmp.nosex)
  
  pvalues <- c(subset(lm.full, term=="Age")$p.value,
               subset(lm.full, term=="Sex")$p.value,
               subset(lm.full, term=="y")$p.value,
               subset(lm.sex, term=="y")$p.value,
               subset(lm.age, term=="y")$p.value)
  
  # partial correlation coefficient equals normalized beta * CONST
  # see http://stats.stackexchange.com/questions/76815/multiple-regression-or-partial-correlation-coefficient-and-relations-between-th
  
  pres <- !is.na(x) & !is.na(y) # must be calculated on the same observations
  tmp.age  <- lm(Age ~ Sex + Generation + y[pres], data=annot.samples[pres,])
  tmp.sex  <- lm(Sex ~ Age + Generation + y[pres], data=annot.samples[pres,])
  tmp.y    <- lm(y[pres] ~ Sex + Age + Generation, data=annot.samples[pres,])
  sigma.full <- sd(tmp.full$resid, na.rm = TRUE)
  sigma.age <- sd(tmp.age$resid, na.rm = TRUE)
  sigma.sex <- sd(tmp.sex$resid, na.rm = TRUE)
  sigma.y <- sd(tmp.y$resid, na.rm = TRUE)
  tmp.noage.y <- lm(y[pres] ~ Sex + Generation, data=annot.samples[pres,])
  tmp.nosex.y <- lm(y[pres] ~ Age + Generation, data=annot.samples[pres,])
  sigma.noage <- sd(tmp.noage$resid, na.rm = TRUE)
  sigma.noage.y <- sd(tmp.noage.y$resid, na.rm = TRUE)
  sigma.nosex <- sd(tmp.nosex$resid, na.rm = TRUE)
  sigma.nosex.y <- sd(tmp.nosex.y$resid, na.rm = TRUE)
  
  coefs <- c(subset(lm.full, term=="Age")$estimate * sigma.age / sigma.full,
             subset(lm.full, term=="Sex")$estimate * sigma.sex / sigma.full,
             subset(lm.full, term=="y")$estimate * sigma.y / sigma.full,
             subset(lm.sex, term=="y")$estimate * sigma.noage.y /sigma.noage,
             subset(lm.age, term=="y")$estimate * sigma.nosex.y / sigma.nosex)
  
  return(c(pvalues,coefs))
}

## test for interaction between Age and Sex
anova_tests_int <- function(x,y) {
  # full model without interaction
  lm.x <- tidy(lm(x ~ Age * Sex + Generation, data=annot.samples))
  lm.y <- tidy(lm(y ~ Age * Sex + Generation, data=annot.samples))
  lm.x.y <- tidy(lm(x ~ Age * Sex + Generation + y, data=annot.samples))
  lm.y.x <- tidy(lm(y ~ Age * Sex + Generation + x, data=annot.samples))
  
  pvalues <- c(subset(lm.x, term=="Age:SexM")$p.value,
               subset(lm.y, term=="Age:SexM")$p.value,
               subset(lm.x.y, term=="Age:SexM")$p.value,
               subset(lm.y.x, term=="Age:SexM")$p.value)

  return(pvalues)
}

result.table <- matrix(NA, N[["pairs"]], 32)
colnames(result.table) <- c("p.mRNA_Age.Sex", "p.mRNA_Sex.Age",
                            "r.mRNA_Age.Sex", "r.mRNA_Sex.Age",
                            "p.Prot_Age.Sex", "p.Prot_Sex.Age",
                            "r.Prot_Age.Sex", "r.Prot_Sex.Age",
                            "p.mRNA_Age.SexProt", "p.mRNA_Sex.AgeProt",
                            "p.mRNA_Prot.SexAge", "p.mRNA_Prot.Sex",
                            "p.mRNA_Prot.Age",
                            "r.mRNA_Age.SexProt", "r.mRNA_Sex.AgeProt",
                            "r.mRNA_Prot.SexAge", "r.mRNA_Prot.Sex",
                            "r.mRNA_Prot.Age",
                            "p.Prot_Age.SexmRNA", "p.Prot_Sex.AgemRNA",
                            "p.Prot_mRNA.SexAge", "p.Prot_mRNA.Sex",
                            "p.Prot_mRNA.Age",
                            "r.Prot_Age.SexmRNA", "r.Prot_Sex.AgemRNA",
                            "r.Prot_mRNA.SexAge", "r.Prot_mRNA.Sex",
                            "r.Prot_mRNA.Age",
                            "p.mRNA_Interaction", "p.Prot_Interaction",
                            "p.mRNA_Interaction.Prot", "p.Prot_Interaction.mRNA")    

# interaction tests
print("Testing for interaction between Age and Sex...")
for (i in 1:N[["pairs"]]) {
  if (i %% 100 == 0) print(i)
  result.table[i, 29:32] <- anova_tests_int(expr.mrna[,i], expr.protein[,i])
}

## normalize verything to mean=0, sd=1
## for easy calculation of partial correlation coefficients
for (i in 1:N[["pairs"]]) {
  expr.mrna[,i] <- scale(expr.mrna[,i])
  expr.protein[,i] <- scale(expr.protein[,i])
}
annot.samples$Age <- scale(annot.samples$Age)  
annot.samples$Sex <- scale(as.numeric(factor(annot.samples$Sex)))

# ANOVA for genes
print("Testing for dependence between mRNA and Age/Sex...")
for (i in 1:N[["pairs"]]) {
  if (i %% 100 == 0) print(i)
  result.table[i, 1:4] <- anova_tests_2(expr.mrna[,i])
}

# ANOVA for proteins
print("Testing for dependence between mRNA and Age/Sex...")
for (i in 1:N[["pairs"]]) {
  if (i %% 100 == 0) print(i)
  result.table[i, 5:8] <- anova_tests_2(expr.protein[,i])
}

# ANOVA for mRNA | Protein expression
print("Testing for dependence between mRNA and Age/Sex conditioned on Protein...")
for (i in 1:N[["pairs"]]) {
  if (i %% 100 == 0) print(i)
  result.table[i, 9:18] <- anova_tests_3(expr.mrna[,i], expr.protein[,i])
}

# ANOVA for mRNA | Protein expression
print("Testing for dependence between mRNA and Age/Sex conditioned on Protein...")
for (i in 1:N[["pairs"]]) {
  if (i %% 100 == 0) print(i)
  result.table[i, 19:28] <- anova_tests_3(expr.protein[,i], expr.mrna[,i])
}

# reorder columns - p-values first, corelation coefs second

pcols <- grep("^p", colnames(result.table))
rcols <- grep("^r", colnames(result.table))

result.table <- result.table[,c(pcols,rcols)]


# LOD (Log-Likelihood) statistics -----------------------------------------

## test of submodel: x ~ test + y + covar vs x ~ y + covar
## xy = c(x,y), test = age or sex
## to be run with parApply
question.loddrop <- function(x, y, test, covar) {
  
  # fit models
  lm.m1 <- lm(x ~ test + y + covar)
  lm.m0 <- lm(x ~ y + covar)
  lm.m3 <- lm(x ~ test + covar)
  lm.m2 <- lm(x ~ covar)
  
  # 
  lod <- c((logLik(lm.m1)-logLik(lm.m0))/log(10), (logLik(lm.m3)-logLik(lm.m2)) / log(10))
  return(lod)
}

lods <- matrix(NA, N[["pairs"]], 4)
colnames(lods) <- c("lod.Prot_Age.mRNASex", "lod.Prot_Age.Sex",
                    "lod.Prot_Sex.mRNAAge", "lod.Prot_Sex.Age")

test.age <- model.matrix(~ Age, data=annot.samples)[,-1]
covar.age<- model.matrix(~ Sex + Generation, data=annot.samples)[,-1] 
test.sex <- model.matrix(~ Sex, data=annot.samples)[,-1]
covar.sex<- model.matrix(~ Age + Generation, data=annot.samples)[,-1] 

# LOD statistics
print("Is Age/Sex effect mediated by RNA...")
for (i in 1:N[["pairs"]]) {
  if (i %% 100 == 0) print(i)
  lods[i, 1:2] <- question.loddrop(expr.protein[,i], expr.mrna[,i],
                                   test=test.age,covar=covar.age)
  lods[i, 3:4] <- question.loddrop(expr.protein[,i], expr.mrna[,i],
                                   test=test.sex,covar=covar.sex)
}

# Output ------------------------------------------------------------------

annot.cols <- c("id",	"gene_id",	"symbol",	"chr",	"start",	"end",
                "strand",	"biotype")

output <- cbind(annot.protein[1:N[["pairs"]], annot.cols], 
                result.table, 5,
                round(lods,4),
                round(averages.mrna,4), 
                round(averages.protein,4))

write.csv(output, output.file, row.names=FALSE)
