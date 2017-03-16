all: anova pca

anova: ../kidney2/kidney_anova_table.csv \
       ../heart2/heart_anova_table.csv \
       ../kidney2/kidney_anova_mrna_table.csv \
       ../heart2/heart_anova_mrna_table.csv

pca: ../kidney2/kidney_pca_loadings.csv \
     ../heart2/heart_pca_loadings.csv \
     ../kidney2/kidney_pca_loadings_mrna.csv \
     ../heart2/heart_pca_loadings_mrna.csv

predictions: glmnet

glmnet: ../kidney2/kidney_predict_glmnet.csv \
        ../heart2/heart_predict_glmnet.csv

../kidney2/kidney_anova_table.csv: code/anova_tests.R ../kidney2/R/DO188b_kidney_noprobs.RData
	Rscript code/anova_tests.R ../kidney2/R/DO188b_kidney_noprobs.RData \
	  ../kidney2/kidney_anova_table.csv

../heart2/heart_anova_table.csv: code/anova_tests.R ../heart2/DO189_heart_v2_noprobs.RData
	Rscript code/anova_tests.R ../heart2/DO189_heart_v2_noprobs.RData \
	  ../heart2/heart_anova_table.csv
	  
../kidney2/kidney_anova_mrna_table.csv: code/anova_tests_mrna.R ../kidney2/R/DO188b_kidney_noprobs.RData
	Rscript code/anova_tests_mrna.R ../kidney2/R/DO188b_kidney_noprobs.RData \
	  ../kidney2/kidney_anova_mrna_table.csv

../heart2/heart_anova_mrna_table.csv: code/anova_tests_mrna.R ../heart2/DO189_heart_v2_noprobs.RData
	Rscript code/anova_tests_mrna.R ../heart2/DO189_heart_v2_noprobs.RData \
	  ../heart2/heart_anova_mrna_table.csv

../kidney2/kidney_pca_loadings.csv: code/pca.R ../kidney2/R/DO188b_kidney_noprobs.RData
	Rscript code/pca.R ../kidney2/R/DO188b_kidney_noprobs.RData \
    ../kidney2/kidney_pca_loadings.csv ../kidney2/kidney_pca_predict.csv

../kidney2/kidney_pca_predict.csv: code/pca.R ../kidney2/R/DO188b_kidney_noprobs.RData
	Rscript code/pca.R ../kidney2/R/DO188b_kidney_noprobs.RData \
    ../kidney2/kidney_pca_loadings.csv ../kidney2/kidney_pca_predict.csv

../heart2/heart_pca_loadings.csv: code/pca.R ../heart2/DO189_heart_v2_noprobs.RData
	Rscript code/pca.R ../heart2/DO189_heart_v2_noprobs.RData \
    ../heart2/heart_pca_loadings.csv ../heart2/heart_pca_predict.csv

../heart2/heart_pca_predict.csv: code/pca.R ../heart2/DO189_heart_v2_noprobs.RData
	Rscript code/pca.R ../heart2/DO189_heart_v2_noprobs.RData \
    ../heart2/heart_pca_loadings.csv ../heart2/heart_pca_predict.csv
    
../kidney2/kidney_pca_loadings_mrna.csv: code/pca_mrna.R ../kidney2/R/DO188b_kidney_noprobs.RData
	Rscript code/pca_mrna.R ../kidney2/R/DO188b_kidney_noprobs.RData \
    ../kidney2/kidney_pca_loadings_mrna.csv ../kidney2/kidney_pca_predict_mrna.csv

../kidney2/kidney_pca_predict_mrna.csv: code/pca_mrna.R ../kidney2/R/DO188b_kidney_noprobs.RData
	Rscript code/pca_mrna.R ../kidney2/R/DO188b_kidney_noprobs.RData \
    ../kidney2/kidney_pca_loadings_mrna.csv ../kidney2/kidney_pca_predict_mrna.csv

../heart2/heart_pca_loadings_mrna.csv: code/pca_mrna.R ../heart2/DO189_heart_v2_noprobs.RData
	Rscript code/pca_mrna.R ../heart2/DO189_heart_v2_noprobs.RData \
    ../heart2/heart_pca_loadings_mrna.csv ../heart2/heart_pca_predict_mrna.csv

../heart2/heart_pca_predict_mrna.csv: code/pca_mrna.R ../heart2/DO189_heart_v2_noprobs.RData
	Rscript code/pca_mrna.R ../heart2/DO189_heart_v2_noprobs.RData \
    ../heart2/heart_pca_loadings_mrna.csv ../heart2/heart_pca_predict_mrna.csv

../kidney2/kidney_predict_glmnet.csv: code/predict_glmnet.R ../kidney2/R/DO188b_kidney_noprobs.RData
	Rscript code/predict_glmnet.R ../kidney2/R/DO188b_kidney_noprobs.RData ../kidney2/kidney_predict_glmnet.csv

../heart2/heart_predict_glmnet.csv: code/predict_glmnet.R ../heart2/DO189_heart_v2_noprobs.RData
	Rscript code/predict_glmnet.R ../heart2/DO189_heart_v2_noprobs.RData ../heart2/heart_predict_glmnet.csv