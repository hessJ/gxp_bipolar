
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

# Load Normalized expression data from each study
# Add sample characteristics, covariates, and other factors
# Predict leukocyte abundances with CellMix package
# Identify significnat latent (unmeasured) variables with surrogate variable analysis ('sva') package

data_list = list.files(path = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/normalized_data/",full.names = T,pattern='.txt')
sample_factors = list.files(path = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/sample_factors/",full.names = T)


# read huex probe table
library(biomaRt)
ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
# Anno <- getBM(attributes=c("hgnc_symbol"),filters=c("affy_huex_1_0_st_v2"),mart=ensembl)



# Witt et al. (GSE46416)
cdf = fread("/Users/jonathanhess/Google Drive/mac_storage/TWAS/bd_mega/data/blood/E-GEOD-46416_mania-euthymia/GPL11028-11233.txt",h=T)
cdf = cdf[!cdf$gene_symbol %in% ""]
cdf = cdf[,colnames(cdf) %in% c("ID","gene_symbol"),with=F]
convert = AnnotationDbi::select(org.Hs.eg.db,keys=cdf$gene_symbol,keytype="SYMBOL",columns="ENTREZID")
colnames(convert) = c("gene_symbol","ENTREZID")
convert = convert[match(cdf$gene_symbol, convert$gene_symbol), ]
cdf$ENTREZID = convert$ENTREZID
cdf = cdf[!is.na(cdf$ENTREZID), ]
cdf = cdf[,c(1:2),with=F]
colnames(cdf) = c("PROBEID","SYMBOL")

sf = fread("/Users/jonathanhess/Google Drive/mac_storage/TWAS/bd_mega/data/blood/sample_factors//E-GEOD-46416.sdrf.txt", h = T)
df = fread("/Users/jonathanhess/Google Drive/mac_storage/TWAS/bd_mega/data/blood/normalized_data//E-GEOD-46416_mania-euthymia_normalizedProbes.txt",h=T)

d = strsplit(colnames(df), "_")
colnames(df) = unlist(lapply(id, function(x) x[[1]]))

sf = sf[,grepl("Source Name|Sample_title", colnames(sf)),with=F]
sf$studyID = "Witt et al"
sf$sampleIDs = gsub(" 1", "", sf$`Source Name`)
sf$dx[grepl("Bipolar", sf$`Comment [Sample_title]`)] = "BD"
sf$dx[grepl("Control", sf$`Comment [Sample_title]`)] = "CT"
sf$age = unlist(lapply(strsplit(sf$`Comment [Sample_title]` , ", "),function(x) x[[2]]))
sf$age = gsub("1_DE|2_DE", "", sf$age)

sf=data.frame(sampleIDs = sf$sampleIDs,
              dx = sf$dx,
              age = sf$age,
              studyID = sf$studyID)
colnames(sf) = paste("FACTOR_",colnames(sf),sep="")

m = as.data.frame(df)
mdf = merge(cdf,m,by="PROBEID")

datExpr0 = aggregate(mdf[,-c(1:2),with=F], by=list(mdf$SYMBOL), function(x) median(x,na.rm=T))
rownames(datExpr0) = datExpr0[,1]
datExpr0 = as.data.frame(t(datExpr0[,-1]))
datExpr0 = data.frame(FACTOR_sampleIDs = rownames(datExpr0), datExpr0) # expression df (Rows = subjects, Cols = genes)
datFull = merge(sf, datExpr0, by="FACTOR_sampleIDs") # Combine with factor variables

fwrite(datFull,
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/expr_covs/bd_witt_qc-normalized.txt",
       quote = F, 
       row.names = F,
       sep = "\t")



# Tsuang et al.
sf = fread("/Users/jonathanhess/Google Drive/mac_storage/TWAS/bd_mega/data/blood/sample_factors//Tsuang2005-SZ_BD_Control.csv",h=T)
sf = sf[,colnames(sf) %in% c("sampleIDs","dx","gender","age","ethnicity","medication","studyID"),with=F]
names(sf)[names(sf) %in% "gender"] = "sex"
sf$dx[sf$dx %in% "bd"]  = "BD"
sf$dx[sf$dx %in% "control"] = "CT"
colnames(sf) = paste("FACTOR_",colnames(sf),sep = "")

df = fread("/Users/jonathanhess/Google Drive/mac_storage/TWAS/bd_mega/data/blood/normalized_data/MingTsuang_PMID15645418_normalizedProbes.txt",h=T)
id = lapply(strsplit(colnames(df)[-1], "[.]"), function(x) x[[2]])
colnames(df)[-1] = unlist(id)
df = as.data.frame(df)

convert = select(hgu133plus2.db, keys=df$PROBEID, keytype="PROBEID",columns="SYMBOL")
df = merge(convert, df, by="PROBEID")
df = data.table(df)

datExpr0 = aggregate(df[,-c(1:2),with=F], by=list(df$SYMBOL), function(x) median(x,na.rm=T))
rownames(datExpr0) = datExpr0[,1]
datExpr0 = as.data.frame(t(datExpr0[,-1]))
datExpr0 = data.frame(FACTOR_sampleIDs = rownames(datExpr0), datExpr0) # expression df (Rows = subjects, Cols = genes)

datFull = merge(sf, datExpr0, by="FACTOR_sampleIDs") # Combine with factor variables

fwrite(datFull,
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/expr_covs/bd_tsuang_qc-normalized.txt",
       quote = F, 
       row.names = F,
       sep = "\t")


# Savitz et al.
sf = fread("/Users/jonathanhess/Google Drive/mac_storage/TWAS/bd_mega/data/blood/GSE39653/gene_expression_omnibus_IDs.csv",h=T)
sf = sf[,grepl("sampleID|subjectID|Gender|age|dx|medicated|race|studyID", colnames(sf)),with=F]
colnames(sf) = paste("FACTOR_",colnames(sf),se)


df = fread("/Users/jonathanhess/Google Drive/mac_storage/TWAS/bd_mega/data/blood/normalized_data//GSE39653_normalizedProbes.txt",h=T)
colnames(df)[-1] = sf$FACTOR_sampleIDs

convert = select(illuminaHumanv4.db, keys=df$PROBEID,keytype="PROBEID",columns="SYMBOL")
convert = convert[!is.na(convert$SYMBOL), ]

df = merge(convert,df,by="PROBEID")
df = data.table(df)

datExpr0 = aggregate(df[,-c(1:2),with=F], by=list(df$SYMBOL), function(x) median(x,na.rm=T))
rownames(datExpr0) = datExpr0[,1]
datExpr0 = as.data.frame(t(datExpr0[,-1]))
datExpr0 = data.frame(FACTOR_sampleIDs = rownames(datExpr0), datExpr0) # expression df (Rows = subjects, Cols = genes)
a = datExpr0

ilmn = fread("/Users/jonathanhess/Google Drive/mac_storage/TWAS/bd_mega/data/blood/GSE39653/GSE39653_ILMN.txt",h=T)
colnames(ilmn)[1] = "PROBEID"
ilmn = merge(convert,ilmn,by="PROBEID")
ilmn = data.table(ilmn)

datExpr0 = aggregate(ilmn[,-c(1:2),with=F], by=list(ilmn$SYMBOL), function(x) median(x,na.rm=T))
rownames(datExpr0) = datExpr0[,1]
datExpr0 = as.data.frame(t(datExpr0[,-1]))
datExpr0 = data.frame(FACTOR_sampleIDs = rownames(datExpr0), datExpr0)


datFull = merge(sf, datExpr0, by="FACTOR_sampleIDs") # Combine with factor variables

fwrite(datFull,
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/expr_covs/bd_savitz_qc-normalized.txt",
       quote = F, 
       row.names = F,
       sep = "\t")


# Padmos et al.
sf = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/sample_factors/E-MEXP-1275.sdrf.txt",h=T)
colnames(sf)[1] = "sampleIDs"
sf = sf[,!colnames(sf) %in% "cellType",with=F]
sf$age = as.integer(sf$age)

colnames(sf) = paste("FACTOR_",colnames(sf),sep="")

df = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/normalized_data/E-MEXP-1275_normalizedProbes.txt",h=T)

convert = AnnotationDbi::select(hgu95av2.db,keys=df$PROBEID,keytype="PROBEID",columns="SYMBOL")

df = merge(convert,df,by="PROBEID")

df = df[!is.na(df$SYMBOL), ]
df = data.table(df)

datExpr0 = aggregate(df[,-c(1:2),with=F], by=list(df$SYMBOL), function(x) median(x,na.rm=T))
rownames(datExpr0) = datExpr0[,1]
datExpr0 = as.data.frame(t(datExpr0[,-1]))
datExpr0 = data.frame(FACTOR_sampleIDs = rownames(datExpr0), datExpr0)


datFull = merge(sf, datExpr0, by= "FACTOR_sampleIDs") # Combine with factor variables

fwrite(datFull,
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/expr_covs/bd_padmos_qc-normalized.txt",
       quote = F, 
       row.names = F,
       sep = "\t")



# Clelland et al.
sf = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/sample_factors/GSE46449_sampleInfo.csv",h=T)
sf$replicate = unlist(lapply(strsplit(sf$samples, " "), function(x) x[[3]]))

df = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/normalized_data/GSE46449_normalizedProbes.txt",h=T)

convert = AnnotationDbi::select(hgu133plus2.db, key=df$PROBEID,keytype = "PROBEID", columns="SYMBOL")

df = merge(convert,df,by="PROBEID")

df = data.table(df)

datExpr0 = aggregate(df[,-c(1:2),with=F], by=list(df$SYMBOL), function(x) median(x,na.rm=T))
rownames(datExpr0) = datExpr0[,1]
datExpr0 = as.data.frame(t(datExpr0[,-1]))
datExpr0 = data.frame(FACTOR_sampleIDs = rownames(datExpr0), datExpr0)

# --> mean biological replicates

rep_list = split(sf, sf$replicate)

rep_save = list()
for(y in 1:length(rep_list)){
  repname = rep_list[[y]]$sampleIDs
  repname = paste(repname, collapse="|",sep="")
  repExpr = datExpr0[grepl(repname, datExpr0$FACTOR_sampleIDs), ]
  colmeans = colMeans(repExpr[,-1])
  colExpr = data.frame(FACTOR_sampleIDs = repExpr$FACTOR_sampleIDs[[1]], t(colmeans))
  rownames(colExpr) = colExpr[,1]
  colExpr = colExpr[,-1]
  rep_save[[y]] = colExpr
}
dat = rbindlist(rep_save)

sampleids = lapply(rep_list, function(x) x$sampleIDs[[1]])

df = data.table(FACTOR_sampleIDs = unlist(sampleids), dat)


sf = sf[,colnames(sf) %in% c("dx","sampleIDs","medication","age","sex"),with=F]
colnames(sf) = paste("FACTOR_",colnames(sf),sep="")

datFull = merge(sf, df, by= "FACTOR_sampleIDs") # Combine with factor variables

fwrite(datFull,
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/expr_covs/bd_clelland_qc-normalized.txt",
       quote = F, 
       row.names = F,
       sep = "\t")


# Bousman et al.,
sf = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/GSE18312/GSE18321_factors.txt",h=T)
colnames(sf) = paste("FACTOR_",colnames(sf),sep="")

df = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/normalized_data/GSE18312_normalizedProbes.txt",h=T)
colnames(df) = gsub(".CEL", "", colnames(df))


cdf = fread("~/Downloads/GPL5175-3188.txt",h=T)
probemap = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/normalized_data/GSE18312_PROBEID.txt",h=T)
colnames(probemap) = c("TSID", "PROBEID")

cdf = cdf[!cdf$GB_LIST %in% ""]
cdf = cdf[!cdf$gene_assignment %in% "---"]

cdf = cdf[,c("ID", "gene_assignment"),with=F]
colnames(cdf) = c("PROBEID","gene_assignment")
cdf = merge(cdf, probemap,by="PROBEID")
cdf = cdf[!duplicated(cdf$PROBEID)]
split = strsplit(cdf$gene_assignment, " // ")
names(split) = cdf$TSID

cdf_df = ldply(split, cbind)
colnames(cdf_df) = c("PROBEID","MAP")

genes = genes(keepStandardChromosomes(TxDb.Hsapiens.UCSC.hg19.knownGene))
genes$SYMBOL = select(org.Hs.eg.db, keys=as.character(genes$gene_id),keytype="ENTREZID", columns = "SYMBOL")$SYMBOL

cdf_df = cdf_df[cdf_df$MAP %in% genes$SYMBOL, ]
colnames(cdf_df) = c("PROBEID","SYMBOL")

df = merge(cdf_df, df,by="PROBEID")
df = data.table(df)

datExpr0 = aggregate(df[,-c(1:2),with=F], by=list(df$SYMBOL), function(x) median(x,na.rm=T))
rownames(datExpr0) = datExpr0[,1]
datExpr0 = as.data.frame(t(datExpr0[,-1]))
datExpr0 = data.frame(FACTOR_sampleIDs = rownames(datExpr0), datExpr0)

datFull = merge(sf, datExpr0, by= "FACTOR_sampleIDs") # Combine with factor variables

fwrite(datFull,
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/expr_covs/bd_bousman_qc-normalized.txt",
       quote = F, 
       row.names = F,
       sep = "\t")


# Beech et al.
sf = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/sample_factors/GSE23848_sampleInfo.csv",h=T)
colnames(sf) = paste("FACTOR_",colnames(sf),sep="")

df = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/normalized_data/GSE23848_normalizedProbes.txt",h=T)

cdf = fread("~/Downloads/GPL6106-11578.txt",h=T)
cdf = cdf[cdf$ID %in% df$PROBEID]
cdf = cdf[!cdf$Symbol %in% ""]
cdf = cdf[,c("ID","Symbol"),with=F]
colnames(cdf) = c("PROBEID","SYMBOL")
convert = select(org.Hs.eg.db,keys=cdf$SYMBOL,keytype="SYMBOL",columns="ENTREZID")
convert = convert[!is.na(convert$ENTREZID),]
cdf = merge(cdf,convert,by="SYMBOL")

df = merge(cdf[,c("SYMBOL","PROBEID")], df, by="PROBEID")

df = data.table(df)

datExpr0 = aggregate(df[,-c(1:2),with=F], by=list(df$SYMBOL), function(x) median(x,na.rm=T))
rownames(datExpr0) = datExpr0[,1]
datExpr0 = as.data.frame(t(datExpr0[,-1]))


datFull = data.frame(sf, datExpr0)


fwrite(datFull,
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/expr_covs/bd_beech_qc-normalized.txt",
       quote = F, 
       row.names = F,
       sep = "\t")

