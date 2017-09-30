
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
require(metafor)



genes = unique(mergestats$GeneSymbol)

res_save = list()
loo_save = list()
for( i in 1:length(genes)){
  
  cat("\rMeta-analysis:",i)
  sub = mergestats[mergestats$GeneSymbol %in% genes[[i]], ]
  sub = sub[grepl("dx", sub$term), ]
  
  sub_combn = data.frame(estimate = c(sub$estimate),
                         BetaStd = c(sub$BetaStd),
                         SE = c(sub$std.error),
                         SEstd = c(sub$SEStd),
                         N = c(sub$N),
                         N_case = c(sub$N_cases),
                         N_control = c(sub$N_controls),
                         studyID = sub$studyID)
  
  if(nrow(sub_combn) < 2) next
  
  
  if(nrow(sub) <= 1) next
  
  direction = sign(sub_combn$estimate)
  direction = ifelse(direction == 1, "+", "-")
  direction = paste(direction, collapse= "", sep = "")
  
  res = NA
  
  res = try(suppressWarnings(metafor::rma(yi = sub_combn$estimate,
                                          sei = sub_combn$SE, 
                                          verbose = FALSE,
                                          weighted = T, 
                                          method = "DL")))
  
  res.std = try(suppressWarnings(metafor::rma(yi = sub_combn$BetaStd,
                                              sei = sub_combn$SEstd, 
                                              verbose = FALSE,
                                              weighted = T, 
                                              method = "DL")))
  
  # try to obtain results with fewer input studies
  if( unique(grepl("Error",res)) == T ) {
    
    sub_combn = sub_combn[sub_combn$N > 30, ];
    
    res = try(suppressWarnings(metafor::rma(yi = sub_combn$estimate,
                                            sei = sub_combn$SE, 
                                            verbose = FALSE,
                                            weighted = TRUE, 
                                            method = "FE")));
    
    res.std = try(suppressWarnings(metafor::rma(yi = sub_combn$BetaStd,
                                                sei = sub_combn$SEstd, 
                                                verbose = FALSE,
                                                weighted = TRUE, 
                                                method = "FE")))
    
  } 
  
  if( unique(grepl("Error",res)) == T ) next
  
  if(nrow(sub_combn) > 1){
    loo_tmp = leave1out(res)
    loo_save[[i]] = data.frame(studyID = sub_combn$studyID, loo_tmp, GeneSymbol = genes[[i]])
  }
  
  res = data.frame(GeneSymbol = sub$GeneSymbol[[1]],
                   Direction = direction, 
                   Nstudy = nrow(sub_combn),
                   Nsample = sum(sub_combn$N),
                   Ncase = sum(sub_combn$N_case),
                   Ncontrol = sum(sub_combn$N_control),
                   Log2FC = res$beta, 
                   SE = res$se, 
                   BetaStd = res.std$beta,
                   SEstd = res.std$se,
                   P = res$pval,
                   HetQ = res$QE,
                   HetP = res$QEp,
                   study_kept = paste(sub_combn$studyID, collapse =  ", ", sep = ""))
  
  res_save[[i]] = res
}

res_df = ldply(res_save)
