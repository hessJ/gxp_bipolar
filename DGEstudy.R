
# load these packages (install if needed)
require(affy)
require(oligo)
require(plyr)
require(ggplot2)
require(data.table)
require(gcrma)
library(limma)
require(pd.hg.u133.plus.2)
require(pd.hg.u133a)
require(pd.hg.u95av2)
require(CellMix)
require(AnnotationDbi)
require(sva)
require(org.Hs.eg.db)
require(hgu133a.db)
require(hgu133plus2.db)
require(illuminaHumanv4.db)
require(hgu95av2.db)
require(hugene10stprobeset.db)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)


datExprFull = data.table(collapse)
study_id = unique(datExprFull$FACTOR_studyID)

datExprFull$FACTOR_dx = factor(datExprFull$FACTOR_dx, levels=c("BD","CT"))

save_results = list()
for( i in 1:length(study_id)){
  
  cat("\nStudy-wise differential expression analysis:",i)
  expr.tmp = datExprFull[datExprFull$FACTOR_studyID %in% study_id[[i]]]
  
  N = table(expr.tmp$FACTOR_dx)
  Nca = N[[1]]
  Nco = N[[2]]
  cat("\n   Detected", Nca,"cases and",Nco,"controls")
  
  # calculate missingness column-wise
  colmisrate = colSums(is.na(expr.tmp[,!grepl("FACTOR_",colnames(expr.tmp)),with=F]))/nrow(expr.tmp) # proportion of missingness
  colexclude = colmisrate[colmisrate > 0.0001] # missingness filter 
  
  expr.tmp = expr.tmp[,!colnames(expr.tmp) %in% names(colexclude),with = F]
  
  # statistical analysis (via linear model)
  y = expr.tmp[,colnames(expr.tmp) %in% names(colmisrate),with=F]
  x = expr.tmp[,!colnames(expr.tmp) %in% names(colmisrate),with = F]
  
  # identify columns in y that are non-numeric
  class_check = lapply(y,class)
  class_check = unlist(class_check)
  wrong_class = which(class_check != "numeric")
  
  # add wrong class to demographics
  if(length(wrong_class)>0){
    x = data.table(x, y[,colnames(y) %in% names(wrong_class),with=F]);
    y = y[,!colnames(y) %in% names(wrong_class),with=F]}
  
  # remove columns with missingness ( > 50%)
  colmis = colSums(is.na(x)) >= nrow(x)*.25
  colmis = which(colmis == TRUE)
  if(length(colmis) > 0){x = x[,!colnames(x) %in% names(colmis),with=F]}
  
  # remove rows with some missingness
  rowmis = rowSums(is.na(x))
  rowmis = which(rowmis > 0)
  if(length(rowmis) > 0){x = x[-rowmis]; y = y[-rowmis]}
  
  # Model matrix (basic)
  predictors = x[,grepl("dx|sex|age|psychosis|medicated|race", colnames(x)),with=F]
  
  counts = lapply(predictors,table)
  count_class = unlist(lapply(counts,length))
  count_class = count_class[count_class <= 1]
  
  if(length(count_class) > 0){predictors = predictors[,!colnames(predictors) %in% names(count_class),with=F]}
  
  
  # remove genes with low variance
  var_filter = lapply(y, sd)
  low_var_filter = var_filter[is.na(var_filter) | var_filter < .001]
  
  if(length(low_var_filter) > 0){y = y[,!colnames(y) %in% names(low_var_filter),with=F]}
  
  N = table(x$FACTOR_dx)
  Nca = N[[1]]
  Nco = N[[2]]
  cat("\n   Kept", Nca,"cases and",Nco,"controls")
  
  
  # CELLMIX SECTION
  # estimate leukocyte abundance
  cat("\n   Estimating leukocyte abundances with CellMix...")
  exprs = as.data.frame(y) # extract normalized gene expression intensities for subjects
  # colnames(exprs) = gsub("ENTREZID_","",colnames(exprs))
  # convert = select(hgu133a.db, keys = colnames(exprs), keytype="ENTREZID", columns="PROBEID")
  convert = select(hgu133a.db, keys = colnames(exprs), keytype="SYMBOL", columns="PROBEID")
  
  minval = min(apply(exprs, 2, function(x) min(x, na.rm = T))) # find minimum value of full matrix
  exprs = exprs + abs(minval) # add |minimum value| to achieve non-negative matrix
  exprs = exprs[,!colSums(is.na(exprs)) > 0]
  
  exprs = data.frame(SYMBOL = colnames(exprs), t(exprs))
  exprs = merge(convert, exprs, by="SYMBOL")
  exprs = exprs[,!colnames(exprs) %in% "SYMBOL"]
  exprs = exprs[!duplicated(exprs$PROBEID), ]
  exprs = exprs[!is.na(exprs$PROBEID), ]
  rownames(exprs) = exprs$PROBEID
  
  exprs = ExpressionSet(as.matrix(exprs[,-1])) # convert to ExpressionSet object
  
  res = gedBlood(exprs, verbose = T, normalize = T) # run gedBlood algorithm (non-negative matrix factorization)
  
  wb.coef = coef(res)
  wb.coef = as.data.frame(t(wb.coef))
  wb.coef = wb.coef[,!colSums(wb.coef == 0) > 0.98*nrow(wb.coef)] # remove cell types with high missingness
  
  ## PCA reduction of peripheral leukocyte proportions
  design =  model.matrix( ~ -1 + ., predictors)
  
  leuk.pca = prcomp(wb.coef, scale = T, center = T)
  leuk.pca.sd = (leuk.pca$sdev^2)/sum(leuk.pca$sdev^2)
  
  leuk.pc = as.data.frame(leuk.pca$x)
  leuk.pc = leuk.pc[gsub("X","",rownames(leuk.pc)) %in% rownames(design), ]
  
  lmod = lm(as.matrix(leuk.pc) ~ design)
  lmod = summary(lmod)
  coefs = lapply(lmod, function(x) broom::tidy(x$coefficients))
  names(coefs) = colnames(leuk.pc)
  coefs = ldply(coefs)
  coefs = coefs[grepl("dx", coefs$.rownames), ]
  
  keep_cell_pc = coefs$.id[coefs$Pr...t.. < .05] # retain significant leukocyte factors
  
  if(length(keep_cell_pc) > 0){
    ltrap = data.frame(leuk.pc[,colnames(leuk.pc)%in%keep_cell_pc])
    colnames(ltrap) = keep_cell_pc
    predictors = data.frame(predictors, ltrap)
  }
  
  ## Surrogate variable analysis - default method 
  exprs = ExpressionSet(as.matrix(t(y)))
  mod = model.matrix(~ ., data=predictors) # model with known factors and covariates
  mod0 = model.matrix(~1,data=predictors) # intercept only model
  n.sv = num.sv(exprs(exprs),mod,method="be") # number of significant surrogate variables.. Using "be" works well with smaller samples
  
  if(n.sv > 0){
    svobj = sva(exprs(exprs),mod,mod0,n.sv=n.sv)
    svdf = data.frame(NULL)
    svdf = as.data.frame(svobj$sv)
    colnames(svdf) = paste("SV",1:ncol(svdf), sep = "")
    
    predictors = data.frame(predictors, svdf)
  }
  
  # final design matrix for differential expression analysis
  design = model.matrix( ~ -1 + ., predictors)
  
  # fit the linear model (log2 expression as response variable)
  lmFit = lm(as.matrix(y) ~  design)
  stdLmFit = lm(as.matrix(scale(y)) ~  design)
  # extract summary statistics
  summary = summary(lmFit)
  stdSummary = summary(stdLmFit)
  # format summary statistics into a table
  tidy_table = lapply(summary, function(x) broom::tidy(x))
  names(tidy_table) = colnames(y)
  tidy_stdTable = lapply(stdSummary, function(x) broom::tidy(x))
  names(tidy_stdTable) = colnames(y)
  # squash into a big table
  big_table = ldply(tidy_table)
  big_table$N_cases = Nca
  big_table$N_controls = Nco
  big_table$N = nrow(x)
  stdTable = ldply(tidy_stdTable)
  big_table$BetaStd = stdTable$estimate
  big_table$SEStd = stdTable$std.error
  
  nGenes = length(unique(big_table$.id))
  
  names(big_table)[names(big_table) %in% ".id"] = "GeneSymbol"
  save_results[[i]]  = big_table
  names(save_results)[[i]] = study_id[[i]]
  
}


## merge study-wise sum stats into table
mergestats = ldply(save_results)
names(mergestats)[names(mergestats) %in% ".id"] = "studyID"
mergestats$term = gsub("design", "", mergestats$term)

# check that t-values are equal from log2/se and z/se
# tval = mergestats$estimate[mergestats$term %in% "dxBD"]/mergestats$std.error[mergestats$term %in% "dxBD"]
# std.tval = mergestats$BetaStd[mergestats$term %in% "dxBD"]/mergestats$SEStd[mergestats$term %in% "dxBD"]

graph_df = mergestats
graph_df$term[graph_df$term %in% "psychosisYes"] = "Psychosis (Yes/No)"
graph_df$term[graph_df$term %in% "sexM"] = "Sex (Male/Female)"
graph_df$term[graph_df$term %in% "tobaccoYes"] = "Tobacco (Yes/No)"
graph_df$term[graph_df$term %in% "dxBD"] = "Diagnosis (BD v. CT)"
graph_df$term[graph_df$term %in% "racenotEUR"] = "European (Yes/No)"
graph_df$term[graph_df$term %in% "age"] = "Age (years)"
non_sv = graph_df[!grepl("SV",graph_df$term), ]
sv_df = graph_df[grepl("SV",graph_df$term), ]
sv_df$term = factor(sv_df$term,levels=paste("SV",1:20,sep=""))

combn = ldply(list(non_sv,sv_df))
combn$term = factor(combn$term, levels=unique(combn$term))

g = ggplot(combn[!grepl("Intercept", ignore.case=T, combn$term), ], aes(x = term, y = -log10(p.value), fill = factor(term))) +
  ylab(expression(paste("Differential expression, -log"[10],"(P-value)"))) + 
  facet_wrap(~studyID, ncol = 4) +
  xlab(NULL) + 
  geom_violin() + 
  guides(fill = FALSE) + 
  theme_bw() +
  geom_hline(yintercept = -log10(.05), col = "red", linetype ="dashed", lwd = 0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_boxplot(width=0.05, outlier.shape = NA, fill = "lightgrey", col = "black")

g
