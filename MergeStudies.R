
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


# Combine normalized data with factors 

norm_files = list.files(path = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/expr_covs/",full.names=T)
norm_files

read = lapply(norm_files, function(x) as.data.frame(fread(x,h=T)))



for( i in 1:length(read)){
  read[[i]]$FACTOR_studyID = basename(norm_files[[i]])
  tmp = read[[i]]
  names(tmp)[names(tmp) %in% "FACTOR_Gender"] = "FACTOR_sex"
  names(tmp)[names(tmp) %in% "FACTOR_medication"] = "FACTOR_medicated"
  names(tmp)[names(tmp) %in% "FACTOR_ethnicity"] = "FACTOR_race"
  tmp = tmp[,!colnames(tmp) %in% c("FACTOR_sampleIDs.1")]
  read[[i]] = data.table(tmp)
}

collapse = ldply(read) # collapse matrices into one data frame

# harmonize factor levels for categorical data 
collapse$FACTOR_dx[collapse$FACTOR_dx %in% "BP"] = "BD"
collapse = collapse[collapse$FACTOR_dx %in% c("BD","CT"), ] # retain BD and CT subjects for analysis
collapse$FACTOR_medicated[collapse$FACTOR_medicated %in% ""] = NA
collapse$FACTOR_medicated[collapse$FACTOR_medicated %in% "yes"] = "Yes"
collapse$FACTOR_medicated[collapse$FACTOR_medicated %in% "no"] = "No"
collapse$FACTOR_sex[collapse$FACTOR_sex %in% "Male"] = "M"
collapse$FACTOR_sex[collapse$FACTOR_sex %in% "Female"] = "F"

table(collapse$FACTOR_race) # checks out
table(collapse$FACTOR_sex) # checks out
table(collapse$FACTOR_medicated) # checks out

