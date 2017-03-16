# Age/Sex GLMNET Predictions
# 
# Use the glmnet package to predict Age and Sex.
#
# Usage:
#
# Rscript predict_glmnet.R input.RData output_predictions.csv

# Options -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# for kidney
# args = list("../kidney2/R/DO188b_kidney_noprobs.RData", 
#              "../kidney2/kidney_predict_glmnet.csv")

# two arguments expected (one input and two output files)
stopifnot(length(args)==2)

input.file = args[[1]]
stopifnot(file.exists(input.file))
output.file = args[[2]] 

load(input.file)

# Glmnet ------------------------------------------------------------------

library(glmnet)
library(dplyr)

# convert sex to numeric
annot.samples$Sex <- as.numeric(factor(annot.samples$Sex))

# which numeric columns from annot.samples are we going to predict
to.be.predicted <- c("Age", "Sex")

# other covariates to use besides Generation
covars <- list(Age = model.matrix(~0+Sex, data=annot.samples),
               Sex = model.matrix(~0+Age, data=annot.samples))

# level of prediction
levels <- c("prot", "mrna", "all_mrna") # all_mrna - to be added later

# factors to be used
xs = list(prot = expr.protein[,1:N[["complete"]]],
         mrna = expr.mrna[,1:N[["complete"]]],
         all_mrna = expr.mrna)  

output.list <- list()
  
for (p in to.be.predicted)
  for (l in levels) {
    
    print(paste(p,l))
    
    # true values
    y = annot.samples[, p]
    # predictors
    x = xs[[l]]
    # other predictors
    covar <- cbind(covars[[p]], 
                   model.matrix(~0+Generation, data=annot.samples))
    
    X <- cbind(x, covar)
    
    # object for the predicted values
    pred <- rep(NA, nrow(annot.samples))
      
    # let predict each sample from all the others
    for (i in 1:nrow(annot.samples)) {
      if (i%%50==0) print(i)
      cvfit = cv.glmnet(X[-i,], y[-i])
      cf <- coef(cvfit, s = "lambda.min")
      pred[i] <- colSums(cf[-1] * cbind(X[i,]))+ cf[1]  
    }
    
    output.name <- paste(p, l, sep=".") # column name
    output.list[[output.name]] <- pred
  }


# Output results ----------------------------------------------------------

# write output into a file
output <- cbind(annot.samples, as.data.frame(output.list))
write.csv(output, file=output.file, row.names=FALSE)
