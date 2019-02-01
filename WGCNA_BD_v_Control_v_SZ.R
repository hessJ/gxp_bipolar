
# load packages
require(data.table)
require(WGCNA)
require(plyr)
require(ggplot2)
require(igraph)
require(CellMix)
require(hgu133a.db)

# rm(list=ls())

# manhattan plot
require(org.Hs.eg.db)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes = as.data.frame(genes)
genes$SYMBOL = select(org.Hs.eg.db, 
                      keys=as.character(genes$gene_id), 
                      keytype="ENTREZID", 
                      columns="SYMBOL")$SYMBOL

convert = genes
colnames(convert) = c("CHR", "START", "END", "width", "strand", "ENTREZID", "GeneSymbol")
autosomes = convert[convert$CHR %in% c(paste("chr",1:22,sep="")),]
sexGenes = convert[convert$CHR %in% c("chrX", "chrY"), ]

## import BD-controls data
# mydata = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_7studiescombined_zscaled.txt", sep= "\t",h = T, stringsAsFactors = F)
# convert column names to gene symbols
# convert = AnnotationDbi::select(org.Hs.eg.db,keys=as.character(gsub("ENTREZID_", "", colnames(mydata)[grepl("ENTREZID",colnames(mydata))])), keytype="ENTREZID",columns="SYMBOL")
# convert$ENTREZID = paste("ENTREZID_",convert$ENTREZID,sep="")
# replace ENTREZ ids with HGNC symbols
# colnames(mydata)[grepl("ENTREZID",colnames(mydata))] = convert$SYMBOL

mydata = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/expr_covs/BD_allstudies_standardized.txt",h=T)
# mydata = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_combat_7studies.txt",h=T)
genes = colnames(mydata)[!grepl("FACTOR", colnames(mydata))]
colnames(mydata) = gsub("FACTOR_","",colnames(mydata))
mydata$sex = ifelse(grepl('witt', ignore.case = T, mydata$studyBatch) == T, 'male', mydata$sex)

## import SZ-controls data
szdata = fread("~/Google Drive/mac_storage/TWAS/scz_mega/blood_final_4.csv", header = T, sep = ",", stringsAsFactors = F)
names(szdata)[names(szdata) %in% "dx_status"] = "dx"
names(szdata)[names(szdata) %in% "study"] = "studyID"
names(szdata)[names(szdata) %in% "match"] = "sampleIDs"

# find columns that are in both SZ and BD data sets
# int = intersect(colnames(szdata)[-c(1:12)], colnames(mydata)[-c(1:12)])
int = intersect(colnames(szdata), colnames(mydata))


## combine using intersected columns 
datExpr0 = mydata
datExpr1 = szdata

datList = list(datExpr0, datExpr1)

datCombine = ldply(datList)
datCombine = datCombine[,colnames(datCombine) %in% int] # note that this includes genes and some shared covariates 

datCombine = datCombine[!is.na(datCombine$sex), ]

# add subject ID to rownames
rownames(datCombine) = datCombine$sampleIDs

# grab columns containing microarray data
datExpr = datCombine[,6:ncol(datCombine)]

datExpr = datExpr[,!colSums(is.na(datExpr))>nrow(datExpr)*.10]


keepSamples = datCombine[,c(1:5)]
keepSamples = keepSamples[which(gsg$goodSamples == TRUE),]

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK


keepGenes = colnames(datExpr)



cat("\n****** Building WGCNA modules for data set:******\n")
phenos = as.data.frame(keepSamples)
rownames(phenos) = keepSamples$sampleIDs
phenos = phenos[phenos$sampleIDs %in% rownames(datCombine), ]
phenos[phenos %in% ""] = NA
phenos = phenos[!rowSums(is.na(phenos))>0,]
phenos$sex = ifelse(phenos$sex %in% c("F","f","female"), "female", phenos$sex)
phenos$sex = ifelse(phenos$sex %in% c("m","M","male"), "male", phenos$sex)
phenos$dx = ifelse(phenos$dx %in% c("CT","control"), "control", phenos$dx)
rownames(phenos) = phenos$sampleIDs

table(phenos$dx, phenos$sex)

rownames(datExpr) = datCombine$sampleIDs
datExpr = datExpr[,colnames(datExpr) %in% keepGenes]
datExpr = datExpr[,order(colnames(datExpr))]
# datExpr = datExpr[rownames(datExpr) %in% rownames(phenos), ]

# Surrogate variable analysis - default method
nonMissingExpr = datExpr[,!colSums(is.na(datExpr))>1]
nonMissingExpr = nonMissingExpr[rownames(nonMissingExpr) %in% rownames(phenos),]
exprs = ExpressionSet(as.matrix(t(nonMissingExpr)))
mod = model.matrix(~ dx + sex + age, data=phenos) # model with known factors and covariates
mod0 = model.matrix(~1,data=phenos) # intercept only model

require(sva)
n.sv = num.sv(exprs(exprs), mod, method="leek", B = 100) # number of significant surrogate variables.. Using "be" works well with smaller samples

if(n.sv > 0){
  svobj = sva(exprs(exprs),mod,mod0,n.sv=n.sv)
  svdf = data.frame(NULL)
  svdf = as.data.frame(svobj$sv)
  colnames(svdf) = paste("SV",1:ncol(svdf), sep = "")
  phenos = data.frame(svdf, phenos)
}

table(phenos$dx, phenos$sex)

exprs = as.data.frame(datExpr) # extract normalized gene expression intensities for subjects
exprs = exprs[,colSums(is.na(exprs))==0]

exprs = exprs

# colnames(exprs) = gsub("ENTREZID_","",colnames(exprs))
# convert = select(hgu133a.db, keys = colnames(exprs), keytype="ENTREZID", columns="PROBEID")
convert = AnnotationDbi::select(hgu133a.db, keys = colnames(exprs), keytype="SYMBOL", columns="PROBEID")

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

res = gedBlood(exprs, verbose = TRUE, normalize = TRUE) # run gedBlood algorithm (non-negative matrix factorization)

wb.coef = coef(res)
wb.coef = as.data.frame(t(wb.coef))

rownames(wb.coef) = gsub("X", "", rownames(wb.coef))
wb.coef = wb.coef[rownames(wb.coef) %in% phenos$sampleIDs, ]

phenos$dx = relevel(as.factor(phenos$dx), ref = 'control')
fit = lm(as.matrix(wb.coef) ~ phenos$dx)
coefs = summary(fit)
coefs = lapply(coefs, function(x) broom::tidy(x$coefficients))
coefs = ldply(coefs)
coefs = coefs[grepl("phenos", coefs$.rownames),]
coefs = coefs[coefs$Pr...t.. < 0.05, ]

boxplot(wb.coef$`Tc act` ~ phenos$dx, outline = FALSE)

# correlation of blood cell concentrations to surrogate variables
cor_res = psych::corr.test(wb.coef, phenos[,colnames(phenos) %in% c("SV1","SV2","SV3","SV4")], adjust = 'none')
cor_p = cor_res$p
cor_r = cor_res$r
cor_r[cor_p]

require(corrplot)

cor_r[which(cor_p > 0.05)] = NA

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/Cor_cellmix_sv.png",res=300,units="in",height=5,width=8)
corrplot(t(cor_r), 
         p.mat = t(cor_p), sig.level = 0.05, 
         insig = 'blank',
         method = "color", 
         outline = T,  
         tl.srt = 90,
         addgrid.col = "darkgray",
         rect.col = "black", 
         rect.lwd = 5,
         cl.pos = "b",
         tl.col = "black",
         tl.cex = 0.75, 
         cl.cex = 0.75, 
         addCoef.col = "black",
         number.digits = 2, 
         number.cex = 0.75, 
         col = colorRampPalette(c("darkred","white","midnightblue"))(100))
dev.off()

## t-SNE (what clusters are identified through t-SNE embedding?)
library(Rtsne)

train = datExpr[rownames(datExpr) %in% rownames(phenos), !colSums(is.na(datExpr))>0]

k.max <- 5
wss <- sapply(1:k.max, 
              function(k){kmeans(train, k, nstart=50,iter.max = 15 )$tot.withinss})

BIC2 <- function(fit){
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(data.frame(AIC = D + 2*m*k,
                    BIC = D + log(n)*m*k))}

bicstat = list()
for(k in 1:5){
  fit = kmeans(train, k, nstart=50,iter.max = 15 )
  bicstat[[k]] = BIC2(fit)
}


bic = unlist(lapply(bicstat, function(x) x$BIC))
aic = unlist(lapply(bicstat, function(x) x$AIC))

bic = data.frame(Val = bic, Class = 'bic', Ncluster = 1:length(bic))
aic = data.frame(Val = aic, Class = 'aic', Ncluster = 1:length(aic))

crt = ldply(list(bic, aic))

ggplot(crt, aes(y = Val, x = Ncluster, col = toupper(Class))) + 
  geom_line() +
  ylab("Value (BIC/AIC)") + 
  xlab("K-means") +
  scale_color_discrete(NULL) +
  theme_bw() + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))


tsne <- Rtsne(train, 
              dims = 2, 
              perplexity=30, 
              verbose=TRUE, 
              Theta = 0.5,
              max_iter = 500)

tPlot = data.frame(tsne$Y)
tPlot$Dx = toupper(phenos$dx)
tPlot$studyID = as.character(phenos$studyID)

ggplot(tPlot, aes(x = X1, y = X2, col = Dx)) + 
  geom_point() + 
  theme_classic() + 
  xlab('Componenet 1') + 
  ylab("Component 2")


# soft-threshold power analyis for WGCNA
powers = c(c(1:10), seq(from = 12, to=20, by=2))

iters = 5
sft1 = list()
for( y in 1:iters){
  random = datExpr[,sample(colnames(datExpr), 500, replace = FALSE)]
  sft1[[y]] = pickSoftThreshold(random, 
                                powerVector = powers, 
                                corFnc="bicor",
                                networkType="signed")
}
sft1_power = lapply(sft1, function(x) x[[1]])
topPower = names(sort(table(unlist(sft1_power)),decreasing = T))
sft1_in = lapply(sft1, function(x) x$fitIndices)
names(sft1_in) = iters
sft1_df = ldply(sft1_in)
sft1_avg = aggregate(sft1_df$SFT.R.sq, by = list(sft1_df$Power), mean)
sft1_se = aggregate(sft1_df$SFT.R.sq, by = list(sft1_df$Power), function(x) sd(x)/sqrt(length(x)))
sft1_avg$CI_LOWER = sft1_avg$x - ( 1.96 * sft1_se$x)
sft1_avg$CI_HIGH = sft1_avg$x + ( 1.96 * sft1_se$x)
names(sft1_avg)[names(sft1_avg) %in% "x"] = "mean"
sft1_avg$label = ifelse(sft1_avg$Group.1 == topPower[[i]], "*", NA)

datGroups = unique(phenos$dx)
datGroups = paste(datGroups, collapse = ", ",sep ="")

r2 = ggplot(sft1_avg, aes(x = factor(Group.1), y = mean)) +
  geom_bar(stat = 'identity', fill = heat.colors(15), colour = 'black') + 
  theme_classic() +
  xlab("Soft threshold power") +
  ylab(expression(paste("Model fit ",italic(R)^2))) +
  geom_errorbar(aes(ymin = CI_LOWER, ymax = CI_HIGH), width=0.2) +
  geom_hline(aes(yintercept = .8),col='navyblue',linetype='dashed') + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))

print(r2)

png(paste("~/Google Drive/mac_storage/TWAS/WGCNA_BD-SZ-CT.png",sep=""),res=300,units="in",height=6,width=6)
print(r2)
dev.off()


pwr = pickSoftThreshold(datExpr[,colSums(is.na(datExpr))==0], 
                  powerVector = powers, 
                  corFnc="bicor",
                  networkType="signed")

pwr = pwr$fitIndices

ggplot(pwr, aes(x = factor(Power), y = SFT.R.sq)) +
  geom_bar(stat = 'identity', fill = heat.colors(15), colour = 'black') + 
  theme_classic() +
  xlab("Soft threshold power") +
  ylab(expression(paste("Model fit ",italic(R)^2))) +
  geom_hline(aes(yintercept = .8),col='navyblue',linetype='dashed') + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))



# 3. set parameters for network algorithm
sft.power = 10; ## SFT 8, R-squared = 0.88013307 (SE = 0.05562321)
deepSplit = 2;
minModuleSize = 30;

# 4.one-step automated network construction
net_a = blockwiseModules(datExpr, 
                                 power = sft.power, 
                                 networkType = "signed",
                                 deepSplit= deepSplit, 
                                 TOMType = "signed",
                                 minModuleSize = minModuleSize, 
                                 minCoreKME = 0.5, 
                                 minCoreKMESize = minModuleSize/3,
                                 minKMEtoStay=0, 
                                 reassignThreshold = 1e-6, 
                                 mergeCutHeight = 0.25, 
                                 detectCutHeight = 0.995,
                                 numericLabels = TRUE, 
                                 pamRespectsDendro = FALSE, 
                                 pamStage = TRUE,
                                 saveTOMs = TRUE, 
                                 verbose = 3,
                                 maxBlockSize = 5000)

net = net_a

# saveRDS(net , file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/WGCNA_network_BD-SZ-CT.pwr7.Rdata")
# saveRDS(net , file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/WGCNA_network_BD-SZ-CT.pwr8.Rdata")
saveRDS(net , file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/freeze_WGCNA_network_BD-SZ-CT.pwr10.Rdata")
# saveRDS(net, file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/WGCNA_network_BD-SZ-CT.pwr9.Rdata")


# 6. extract network meta-data and eigengenes
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
MEs = data.frame(sampleIDs = rownames(datExpr), MEs)

metadata = data.frame(colors = moduleColors, labels = paste("ME",net$colors, sep=""))
metadata = metadata[!duplicated(metadata$colors), ]

mergedColors = labels2colors(net$colors)


png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD-SZ-CT_dendrogram.png", res=300, units="in",height = 5.5, width = 8)

plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()


module = list()
colors = unique(moduleColors)
for(y in 1:length(colors)){
  genesInModule = colnames(datExpr)[which(moduleColors %in% colors[[y]])]
  module[[y]] = data.frame(color = colors[[y]], label = unique(net$colors)[[y]], symbol = genesInModule)
}
module = ldply(module)


fwrite(data.table(module), 
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD-SZ-CT_WGCNA_module-membership.txt", 
       quote = F, row.names = F, sep = "\t")



## ===== STATISTICAL ANALYSIS


# combine MEs and phenos 
phenos = phenos[!is.na(phenos$sex), ]
datAll = merge(phenos, MEs, by='sampleIDs',all=T)
datAll$dx = as.character(datAll$dx)
datAll$dx[grepl("BD|scz", datAll$dx)] = 'psychosis'

fwrite(datAll, file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD-CT-SZ_WGCNA-network.txt", quote = F, row.names = F, sep=  "\t")

y = datAll[,grepl("ME", colnames(datAll))]
lmod = lm(as.matrix(y) ~ age + sex * dx + SV1 + SV2 + SV3 + SV4 + studyID , data=datAll)
sum = summary(lmod)
coefs = lapply(sum, function(x) as.data.frame(x$coefficients))

# add predictor and response columns to coefficient tables
for( i in 1:length(coefs)){
  coefs[[i]]$labels = colnames(y)[[i]]
  coefs[[i]]$predictors = rownames(coefs[[i]])
}
# combine list elements into dataframe
coef_all = plyr::ldply(coefs)

# convert p-values to -log10 scale
coef_all$logP = -log10(coef_all$`Pr(>|t|)`)

# delete unwanted sub-string
coef_all$predictors = gsub("design", "", coef_all$predictors)

# grab effects for diagnostic status and adjust p-values for multiple testing
dx.est = coef_all[coef_all$predictors %in% "dxpsychosis", ]
dx.est = merge(metadata, dx.est, by = "labels")
dx.est = dx.est[!dx.est$labels %in% "ME0",] # discard ME0 --> "grey" module
dx.est = dx.est[order(dx.est$`Pr(>|t|)`, decreasing = F), ]
dx.est$fdr = p.adjust(dx.est$`Pr(>|t|)`, "fdr")
dx.est$bonf = p.adjust(dx.est$`Pr(>|t|)`, "bonferroni")

head(dx.est)

mColorCount = data.frame(table(module$color))
colnames(mColorCount) = c("colors", "nGenes")

dx.est = merge(mColorCount, dx.est, by="colors")
dx.est = dx.est[order(dx.est$`Pr(>|t|)`, decreasing = F), ]

head(dx.est)

fwrite(dx.est, 
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/psychosis_WGCNA_lmresults.txt",
       quote = F, row.names = F, sep = "\t")

dx.est$lab_point = ifelse(dx.est$bonf < .05, as.character(dx.est$labels), NA)

require(ggrepel)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD-SZ-CT.volcano.png",res=300,units="in",height=6,width=6)
ggplot(dx.est, aes( x = Estimate, y = logP)) +
  geom_point(bg = as.character(dx.est$colors),pch=21,col='black', size = 2.5) + 
  theme_classic() +
  geom_hline(yintercept = -log10(.05/nrow(dx.est)), linetype=1,col = 'red') + 
  geom_hline(yintercept = -log10(.05), col = 'orange', linetype=2) + 
  theme(axis.text = element_text(size =15), axis.title=element_text(size=15)) + 
  xlab(expression(paste("Effect size (",beta,")"))) + 
  ylab(expression(paste("-log"[10],"(P-value)"))) +
  geom_text_repel(aes(label = lab_point), size = 5)
dev.off()



dx.est$condition = 'Psychosis vs. Combined CT'
dx.est$LABELS = paste(round(dx.est$Estimate, 4), "\n", "(",round(dx.est$`Std. Error`,3),")", sep = "")
PSYCHOSIS = dx.est


# Psychosis vs. CT (sex-by-diagnosis)

# grab effects for diagnostic status and adjust p-values for multiple testing
dx.est = coef_all[grepl(":", coef_all$predictors), ]
dx.est = merge(metadata, dx.est, by = "labels")
dx.est = dx.est[!dx.est$labels %in% "ME0",] # discard ME0 --> "grey" module
dx.est = dx.est[order(dx.est$`Pr(>|t|)`, decreasing = F), ]
dx.est$fdr = p.adjust(dx.est$`Pr(>|t|)`, "fdr")
dx.est$bonf = p.adjust(dx.est$`Pr(>|t|)`, "bonferroni")

head(dx.est)

mColorCount = data.frame(table(module$color))
colnames(mColorCount) = c("colors", "nGenes")

dx.est = merge(mColorCount, dx.est, by="colors")
dx.est = dx.est[order(dx.est$`Pr(>|t|)`, decreasing = F), ]

head(dx.est)

fwrite(dx.est, 
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/psychosis_WGCNA-sex-by-diagnosis_lmresults.txt",
       quote = F, row.names = F, sep = "\t")

dx.est$lab_point = ifelse(dx.est$bonf < .05, as.character(dx.est$labels), NA)

require(ggrepel)


dx.est$condition = 'Psychosis vs. Combined CT\n(Sex-by-diagnosis)'
dx.est$LABELS = paste(round(dx.est$Estimate, 4), "\n", "(",round(dx.est$`Std. Error`,3),")", sep = "")
PSYSEX = dx.est



# combine MEs and phenos 
datAll = merge(phenos, MEs, by='sampleIDs',all=T)
datAll$dx = factor(datAll$dx)
datAll$dx = relevel(datAll$dx, ref = "control")

y = datAll[,grepl("ME", colnames(datAll))]
lmod = lm(as.matrix(y) ~ age + sex * dx + SV1 + SV2 + SV3 + SV4 + studyID, data=datAll)
sum = summary(lmod)
coefs = lapply(sum, function(x) as.data.frame(x$coefficients))

# add predictor and response columns to coefficient tables
for( i in 1:length(coefs)){
  coefs[[i]]$labels = colnames(y)[[i]]
  coefs[[i]]$predictors = rownames(coefs[[i]])
}
# combine list elements into dataframe
coef_all = plyr::ldply(coefs)

# convert p-values to -log10 scale
coef_all$logP = -log10(coef_all$`Pr(>|t|)`)

# delete unwanted sub-string
coef_all$predictors = gsub("design", "", coef_all$predictors)

# grab effects for diagnostic status and adjust p-values for multiple testing
dx.est = coef_all[coef_all$predictors %in% c("dxBD","dxscz"), ]
dx.est = merge(metadata, dx.est, by = "labels")
dx.est = dx.est[!dx.est$labels %in% "ME0",] # discard ME0 --> "grey" module

dx.est = split(dx.est, dx.est$predictors)

for(i in 1:length(dx.est)){
dx.est[[i]]$fdr = p.adjust(dx.est[[i]]$`Pr(>|t|)`, "fdr")
dx.est[[i]]$bonf = p.adjust(dx.est[[i]]$`Pr(>|t|)`, "bonferroni")
}
dx.est = ldply(dx.est)


mColorCount = data.frame(table(module$color))
colnames(mColorCount) = c("colors", "nGenes")

dx.est = merge(mColorCount, dx.est, by="colors")
dx.est = dx.est[order(dx.est$`Pr(>|t|)`, decreasing = F), ]

head(dx.est)
  
  fwrite(dx.est, 
         file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/controlref-vs-BD-SZ_WGCNA_lmresults.txt",
         quote = F, row.names = F, sep = "\t")
  
  dx.est$lab_point = ifelse(dx.est$bonf < .05, as.character(dx.est$labels), NA)
  
  dx.est$predictors = paste(toupper(gsub("dx", "", dx.est$predictors)), ' vs. Combined CT', sep="")
  
  
  dx.est$condition = dx.est$predictors
  dx.est$LABELS = paste(round(dx.est$Estimate, 4), "\n", "(",round(dx.est$`Std. Error`,3),")", sep = "")
  GROUPS = dx.est
  
  #### 
  
  
  png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/CTref.vs.BD-SZ.volcano.png",res=300,units="in",height=6,width=8)
  ggplot(dx.est, aes( x = Estimate, y = logP)) +
    geom_point(aes(bg = I(colors)), col = 'black', pch=21, size = 2.5) +
    facet_grid(~predictors)  +
    theme_classic() +
    geom_hline(aes(yintercept = -log10(.05/nrow(dx.est))), linetype=1,col = 'red') +
    geom_hline(aes(yintercept = -log10(.05)), col = 'orange', linetype=2) +
    theme(axis.text = element_text(size =15), axis.title=element_text(size=15), strip.text = element_text(size = 15)) + 
    xlab(expression(paste("Effect size (",beta,")"))) + 
    ylab(expression(paste("-log"[10],"(P-value)"))) +
    geom_text_repel(aes(label = lab_point), size = 4, segment.color = 'grey', segment.size = 0.1)
  dev.off()


# grab effects for diagnostic status and adjust p-values for multiple testing
dx.est = coef_all[grepl(":", coef_all$predictors), ]
dx.est = merge(metadata, dx.est, by = "labels")
dx.est = dx.est[!dx.est$labels %in% "ME0",] # discard ME0 --> "grey" module

dx.est = split(dx.est, dx.est$predictors)

for(i in 1:length(dx.est)){
dx.est[[i]]$fdr = p.adjust(dx.est[[i]]$`Pr(>|t|)`, "fdr")
dx.est[[i]]$bonf = p.adjust(dx.est[[i]]$`Pr(>|t|)`, "bonferroni")
}
dx.est = ldply(dx.est)

head(dx.est)

mColorCount = data.frame(table(module$color))
colnames(mColorCount) = c("colors", "nGenes")

dx.est = merge(mColorCount, dx.est, by="colors")
dx.est = dx.est[order(dx.est$`Pr(>|t|)`, decreasing = F), ]


dx.est$predictors = paste(toupper(gsub("sexmale:dx", "", dx.est$predictors)), ' vs. Combined CT', sep="")
dx.est$condition = paste(dx.est$predictors, '\n(Sex-by-diagnosis)')
dx.est$LABELS = paste(round(dx.est$Estimate, 4), "\n", "(",round(dx.est$`Std. Error`,3),")", sep = "")
SEX = dx.est
  
  
  # combine MEs and phenos 
  datAll = merge(phenos, MEs, by='sampleIDs',all=T)
  datAll = datAll[!is.na(datAll$age), ]
  datAll = datAll[datAll$dx %in% c("scz","BD"), ]
  datAll$dx = factor(datAll$dx)
  datAll$dx = relevel(datAll$dx, ref = "BD")
  
  y = datAll[,grepl("ME", colnames(datAll))]
  lmod = lm(as.matrix(y) ~ age + sex*dx + SV1 + SV2 + SV3 + SV4, data=datAll)
  sum = summary(lmod)
  coefs = lapply(sum, function(x) as.data.frame(x$coefficients))
  
  # add predictor and response columns to coefficient tables
  for( i in 1:length(coefs)){
    coefs[[i]]$labels = colnames(y)[[i]]
    coefs[[i]]$predictors = rownames(coefs[[i]])
  }
  # combine list elements into dataframe
  coef_all = plyr::ldply(coefs)
  
  # convert p-values to -log10 scale
  coef_all$logP = -log10(coef_all$`Pr(>|t|)`)
  
  # delete unwanted sub-string
  coef_all$predictors = gsub("design", "", coef_all$predictors)
  
  # grab effects for diagnostic status and adjust p-values for multiple testing
  dx.est = coef_all[coef_all$predictors %in% c("dxBD","dxscz"), ]
  dx.est = merge(metadata, dx.est, by = "labels")
  dx.est = dx.est[!dx.est$labels %in% "ME0",] # discard ME0 --> "grey" module
  
  dx.est = split(dx.est, dx.est$predictors)
  
  for(i in 1:length(dx.est)){
    dx.est[[i]]$fdr = p.adjust(dx.est[[i]]$`Pr(>|t|)`, "fdr")
    dx.est[[i]]$bonf = p.adjust(dx.est[[i]]$`Pr(>|t|)`, "bonferroni")
  }
  dx.est = ldply(dx.est)
  
  head(dx.est)
  
  mColorCount = data.frame(table(module$color))
  colnames(mColorCount) = c("colors", "nGenes")
  
  dx.est = merge(mColorCount, dx.est, by="colors")
  dx.est = dx.est[order(dx.est$`Pr(>|t|)`, decreasing = F), ]
  
  head(dx.est)
  
  fwrite(dx.est, 
         file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/BDvs.SZ_WGCNA_lmresults.txt",
         quote = F, row.names = F, sep = "\t")
  
  dx.est$lab_point = ifelse(dx.est$bonf < .05, as.character(dx.est$labels), NA)
  
  dx.est$predictors = paste("BD vs. ", toupper(gsub("dx", "", dx.est$predictors)),sep="")
  
  
  png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDvs.SZ.volcano.png",res=300,units="in",height=6,width=8)
  
  ggplot(dx.est, aes( x = Estimate, y = logP)) +
    geom_point(aes(bg = I(colors)), col = 'black', pch=21, size = 2.5) +
    facet_grid(~predictors)  +
    theme_classic() +
    geom_hline(aes(yintercept = -log10(.05/nrow(dx.est))), linetype=1,col = 'red') +
    geom_hline(aes(yintercept = -log10(.05)), col = 'orange', linetype=2) +
    theme(axis.text = element_text(size =15), axis.title=element_text(size=15), strip.text = element_text(size = 15)) + 
    xlab(expression(paste("Effect size (",beta,")"))) + 
    ylab(expression(paste("-log"[10],"(P-value)"))) +
    geom_text_repel(aes(label = lab_point), size = 4, segment.color = 'grey', segment.size = 0.1)
  
  dev.off()
  
  

  dx.est$condition = 'BD vs. SZ'
  dx.est$LABELS = paste(round(dx.est$Estimate, 4), "\n", "(",round(dx.est$`Std. Error`,3),")", sep = "")
  CASES = dx.est
  
  
  # grab effects for diagnostic status and adjust p-values for multiple testing
  dx.est = coef_all[grepl(":", coef_all$predictors), ]
  dx.est = merge(metadata, dx.est, by = "labels")
  dx.est = dx.est[!dx.est$labels %in% "ME0",] # discard ME0 --> "grey" module
  
  dx.est = split(dx.est, dx.est$predictors)
  
  for(i in 1:length(dx.est)){
    dx.est[[i]]$fdr = p.adjust(dx.est[[i]]$`Pr(>|t|)`, "fdr")
    dx.est[[i]]$bonf = p.adjust(dx.est[[i]]$`Pr(>|t|)`, "bonferroni")
  }
  dx.est = ldply(dx.est)
  
  head(dx.est)
  
  mColorCount = data.frame(table(module$color))
  colnames(mColorCount) = c("colors", "nGenes")
  
  dx.est = merge(mColorCount, dx.est, by="colors")
  dx.est = dx.est[order(dx.est$`Pr(>|t|)`, decreasing = F), ]
  
  
  dx.est$predictors = paste(toupper(gsub("sexmale:dxscz", "BD", dx.est$predictors)), ' vs. SZ', sep="")
  dx.est$condition = paste(dx.est$predictors, '\n(Sex-by-diagnosis)')
  dx.est$LABELS = paste(round(dx.est$Estimate, 4), "\n", "(",round(dx.est$`Std. Error`,3),")", sep = "")
  SEXCASE = dx.est


# tile plot
SEXDF = ldply(list(PSYSEX, SEXCASE))
# SEXDF = ldply(list(PSYSEX, SEXCASE))
# OTHERDF = ldply(list(PSYCHOSIS, GROUPS, CASES))
OTHERDF = ldply(list(PSYCHOSIS, CASES))

SEXDF$FDRALL = p.adjust(SEXDF$`Pr(>|t|)`, "fdr")
SEXDF$BONFALL = p.adjust(SEXDF$`Pr(>|t|)`, "bonferroni")

OTHERDF$FDRALL = p.adjust(OTHERDF$`Pr(>|t|)`, "fdr")
OTHERDF$BONFALL = p.adjust(OTHERDF$`Pr(>|t|)`, "bonferroni")

graphDF = ldply(list(OTHERDF, SEXDF))
graphDF$condition = factor(graphDF$condition, levels=unique(graphDF$condition))
graphDF$labels = factor(graphDF$labels, levels = paste("ME",1:100,sep=""))
graphDF$LABELS = ifelse(graphDF$FDRALL < .05, paste("*",graphDF$LABELS, sep = ''), graphDF$LABELS)
graphDF$FONT = 1
graphDF$FONT[grepl("[*]", graphDF$LABELS)] = 2
graphDF$condition = gsub("SCZ", "SZ", graphDF$condition)
graphDF$condition = factor(graphDF$condition, levels = rev(unique(as.character(graphDF$condition))))

require(ggthemes)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDSZCT_tile.png",res=300,units="in",height=4,width=10.5)
ggplot(graphDF, aes(x = labels, y = condition, fill = logP)) +
  geom_tile(color = 'black',  size = 0.2, show.legend = TRUE) +
  theme_tufte(base_family = "Helvetica") +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_continuous(expression(paste('-log'[10],'(P-value)')), low = 'royalblue', high = c('green','firebrick3')) +
  theme(axis.ticks = element_line(colour = NA), legend.title = element_text(size=8), axis.text.x = element_text(hjust = 1, angle = 90), axis.text = element_text(size=6)) + 
  theme(legend.position =  'bottom') +
  geom_text(aes(label = LABELS, fontface = FONT), size = 1.75) +
  coord_equal()
dev.off()
  

##= stratify by sex 
 
datAll = merge(phenos, MEs, by='sampleIDs',all=T)
datAll$dx = factor(datAll$dx)
datAll$dx = relevel(datAll$dx, ref = "control")  

sexType = unique(datAll$sex)
sexType = sexType[!is.na(sexType)]
sex.est = list()
for( z in 1:length(sexType)){
  
  strat = datAll[datAll$sex %in% sexType[[z]], ]
  
  y = strat[,grepl("ME", colnames(strat))]
  lmod = lm(as.matrix(y) ~ age + dx + SV1 , data=strat)
  sum = summary(lmod)
  coefs = lapply(sum, function(x) as.data.frame(x$coefficients))
  
  
  # add predictor and response columns to coefficient tables
  for( i in 1:length(coefs)){
    coefs[[i]]$labels = colnames(MEs)[[i]]
    coefs[[i]]$predictors = rownames(coefs[[i]])
  }
  # combine list elements into dataframe
  coef_all = plyr::ldply(coefs)
  
  # convert p-values to -log10 scale
  coef_all$logP = -log10(coef_all$`Pr(>|t|)`)
  
  # delete unwanted sub-string
  coef_all$predictors = gsub("design", "", coef_all$predictors)
  
  # grab effects for diagnostic status and adjust p-values for multiple testing
  dx.est = coef_all[coef_all$predictors %in% c('dxBD','dxscz'), ]
  dx.est = merge(metadata, dx.est, by = "labels")
  dx.est = dx.est[!dx.est$labels %in% "ME0",] # discard ME0 --> "grey" module
  dx.est = dx.est[order(dx.est$`Pr(>|t|)`, decreasing = F), ]
  
  dx.est = split(dx.est, dx.est$predictors)
  
  for(x in 1:length(dx.est)){
    dx.est[[x]]$fdr = p.adjust(dx.est[[x]]$`Pr(>|t|)`, 'fdr')
    dx.est[[x]]$bonf = p.adjust(dx.est[[x]]$`Pr(>|t|)`, 'bonferroni')
  }
  dx.est = ldply(dx.est)
  
  # dx.est$fdr = p.adjust(dx.est$`Pr(>|t|)`, "fdr")
  # dx.est$bonf = p.adjust(dx.est$`Pr(>|t|)`, "bonferroni")
  
  mColorCount = data.frame(table(module$color))
  colnames(mColorCount) = c("colors", "nGenes")
  
  dx.est = merge(mColorCount, dx.est, by="colors")
  dx.est = dx.est[order(dx.est$`Pr(>|t|)`, decreasing = F), ]
  
  head(dx.est)
sex.est[[z]] = dx.est
sex.est[[z]]$Gender = sexType[[z]]
}
sex.est = ldply(sex.est)
  

sex.est$lab_point = ifelse(sex.est$bonf < .05, as.character(sex.est$labels), NA)
sex.est$Gender = stringr::str_to_title(sex.est$Gender)
sex.est$predictors = paste("Combined control vs. ", toupper(gsub("dx", "", sex.est$predictors)),sep="")

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/sex_CTref.vs.BD-SZ.volcano.png",res=300,units="in",height=8,width=9)
ggplot(sex.est, aes( x = Estimate, y = logP)) +
  geom_point(aes(bg = I(colors)),pch=21,col='black',size = 2.5) +
  facet_wrap(~ predictors + Gender, ncol = 2)  +
  theme_classic() +
  geom_hline(aes(yintercept = -log10(.05/nrow(dx.est))), linetype=1,col = 'red') +
  geom_hline(aes(yintercept = -log10(.05)), col = 'orange', linetype=2) +
  theme(axis.text = element_text(size =11), axis.title=element_text(size=15), strip.text = element_text(size = 11)) + 
  xlab(expression(paste("Effect size (",beta,")"))) + 
  ylab(expression(paste("-log"[10],"(P-value)"))) +
  geom_text_repel(aes(label = lab_point), size = 4, segment.color = 'grey', segment.size = 0.1)
dev.off()


## Brain cell signatures enriched in disease-associated modules?

require(ABAEnrichment)

## create input vector with candidate genes 
colors = PSYCHOSIS$colors[PSYCHOSIS$fdr < .05]

fwer_list = list()
for(x in 1:length(colors)){
  
genes = ifelse(module$color == colors[[x]], 1, 0)
names(genes) = module$symbol
genes = genes[genes == 1]

## ----results='hide'------------------------------------------------------
## run enrichment analysis
res=aba_enrich(genes,dataset='5_stages',cutoff_quantiles = 0.9,n_randsets=100)
fwers=res[[1]]
fwers = fwers[order(fwers$mean_FWER,  decreasing = F), ]
fwer_list[[x]] = fwers
}
names(fwer_list) = colors

fwer_df = ldply(fwer_list)

fwer_df = fwer_df[order(fwer_df$mean_FWER, decreasing = F), ]

fwer_df = fwer_df[fwer_df$mean_FWER < .05, ]

sort(table(fwer_df$structure),decreasing = T)


# detailed connectivity plot of genes in significant modules
  
  library(network)
  library(sna)
  library(ggplot2)
  library(GGally)
  
sig_mods = as.character(unique(graphDF$colors[graphDF$FDRALL < .05]))
net.grab = sig_mods
  
  slice.networks = list()
  color.list = list()
  for( y in 1:length(net.grab)){
    grab.genes =  which(moduleColors %in% net.grab[[y]])
    slice.networks[[y]] = datExpr[,grab.genes]
    color.list[[y]] = rep(net.grab[[y]], length(grab.genes))
  }
  slice.networks.df = as.data.frame(do.call(cbind, slice.networks))
  
  adjMat <- abs(cor(slice.networks.df, use = "p")) ## compute network-wide adjacency matrix 
  int_mod <- intramodularConnectivity(adjMat, colors = unlist(color.list)) ## computed gene connectivity from adjacency matrix 
  int_mod = data.frame(gene = rownames(int_mod), int_mod)
  mod.gene.col = data.frame(gene = colnames(datExpr), col = moduleColors)
  int_mod = merge(mod.gene.col, int_mod, by="gene")
  
  hubs = split(int_mod, as.character(int_mod$col))
  hubs = lapply(hubs, function(x) x[order(x$kWithin,decreasing=F), ])
  hubs = lapply(hubs, function(x) x[1,])
  hubs = as.character(do.call(rbind, hubs)$gene)
  
  color.list = which(colnames(datExpr) %in% hubs)
  color.list = moduleColors[color.list]
  
  # network graph of hub genes
  g <- graph.adjacency(as.matrix(adjMat),
                       mode="undirected",
                       weighted=T, diag = F) ## Graph information for plotting
  
  
  length(V(g)) -> vertices
  g=igraph::delete.edges(g, which(E(g)$weight <=.5)) # filter edjges (defined by correlation coefficient)
  
  seq(.1, 1.0, by = .1) -> corr.seq
  
  edges <- list()
  for(x in 1:length(corr.seq)){
    x - 1 -> k
    length(which(E(g)$weight >= corr.seq[[x]])) -> edges[[x]]  }
  names(edges) <- corr.seq
  
  ## Extract proper vertex names
  
  V(g)$colors <- as.character(unlist(color.list))
  
  values = as_data_frame(g, "vertices")
  values$hub_col = NA
  values$hub_col[values$name %in% hubs] = "gold"
  # values$hub_col[values$colors %in% dx.est$colors[dx.est$bonf < .05]] = "gold"
  
  V(g)$pie.color = list(c(as.character(net.grab), "gold"))
  
  val.list = list()
  for(n in 1:nrow(values)){
    cola = list(c(as.character(net.grab), "gold"))[[1]] == values$colors[[n]]
    cola = as.numeric(cola) 
    colb = list(c(as.character(net.grab), "gold"))[[1]] == values$hub_col[[n]]
    colb = as.numeric(colb) 
    colb[is.na(colb)] = 0
    val.list[[n]] = cola + colb
  }
  
  sum.take = unlist(lapply(val.list, sum))
  node.type = rep("circle", length(val.list))
  node.type[which(sum.take > 1)] = "pie"
  
  values$node.type = node.type
  
  V(g)$shape = node.type
  
  V(g)$pie.color[node.type == "circle"] = as.character(unlist(color.list)[node.type == "circle"])
  val.list[node.type == "circle"] = 1
  
  label.colors = rep("black", length(V(g)))
  label.colors[which(V(g)$colors == "black")] = "white"
  
  
  # gDel=igraph::delete.edges(g, which(E(g)$weight <= 0.4))
  gDel = g
  netdf = as_data_frame(minimum.spanning.tree(gDel, weights = E(g)$weights, algorithm = 'prim'), "edges")
  netdf = as_data_frame(gDel)
  net = network(netdf[,1:2], directed = FALSE)
  
  # network plot
  png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood//MST_majorpsychosis.png",res=300,units="in",height=5,width=6)
  ggnet2(net,
         vjust = -0.6,
         label.size = 1.2,
         mode = 'fruchtermanreingold',
         node.color = V(g)$colors,
         size = 1.5, 
         label = TRUE,
         edge.alpha = 0.5)
  dev.off()
  
  
  # gDel = igraph::delete.edges(g, which(E(g)$weight <=.0001))
  gDel = g
  
  hubScore = hub_score(gDel)$vector
  
  hubs = sort(hubScore,decreasing = T)
  
  size = rep(0.15, length(hubScore))
  size[which(names(V(gDel)) %in% names(hubs[1:5]))] = 0.7
  
  cols = rep('lightgrey', length(hubScore))
  cols[which(names(V(gDel)) %in% names(hubs[1:5]))] = 'darkgrey'
  
  vSize = rep(2, length(hubScore))
  vSize[which(names(V(gDel)) %in% names(hubs[1:5]))] = 7
  
  
  
  png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood//TOP_minimumSpanningTree.png",res=300,units="in",height=7,width=10)
  par(mar=c(0,0,0,0)+0.001)
  igraph::plot.igraph(gDel,
                      layout=layout.fruchterman.reingold,
                      # layout=layout.lgl,
                      vertex.size = vSize, 
                      mark.groups = NULL,
                      vertex.shape= "circle",
                      vertex.color = V(g)$colors,
                      vertex.label.cex= size,
                      vertex.label.color = 'black',
                      edge.width=E(g)$weight*2,
                      vertex.label.font = NA,
                      vertex.frame.color= NA,
                      edge.color= "grey",
                      vertex.label.family="Arial")
  dev.off()

  # === next module (ME3)
  
  net.grab = "tan"
  
  slice.networks = list()
  color.list = list()
  for( y in 1:length(net.grab)){
    grab.genes =  which(moduleColors %in% net.grab[[y]])
    slice.networks[[y]] = datExpr[,grab.genes]
    color.list[[y]] = rep(net.grab[[y]], length(grab.genes))
  }
  slice.networks.df = as.data.frame(do.call(cbind, slice.networks))
  
  adjMat <- abs(cor(slice.networks.df, use = "p")) ## compute network-wide adjacency matrix 
  int_mod <- intramodularConnectivity(adjMat, colors = unlist(color.list)) ## computed gene connectivity from adjacency matrix 
  int_mod = data.frame(gene = rownames(int_mod), int_mod)
  mod.gene.col = data.frame(gene = colnames(datExpr), col = moduleColors)
  int_mod = merge(mod.gene.col, int_mod, by="gene")
  
  hubs = split(int_mod, as.character(int_mod$col))
  hubs = lapply(hubs, function(x) x[order(x$kWithin,decreasing=F), ])
  hubs = lapply(hubs, function(x) x[1,])
  hubs = as.character(do.call(rbind, hubs)$gene)
  
  color.list = which(colnames(datExpr) %in% hubs)
  color.list = moduleColors[color.list]
  
  # network graph of hub genes
  g <- graph.adjacency(as.matrix(adjMat),
                       mode="undirected",
                       weighted=T, diag = F) ## Graph information for plotting
  
  
  length(V(g)) -> vertices
  # g=delete.edges(g, which(E(g)$weight <=.2)) # filter edjges (defined by correlation coefficient)
  
  seq(.1, 1.0, by = .1) -> corr.seq
  
  edges <- list()
  for(x in 1:length(corr.seq)){
    x - 1 -> k
    length(which(E(g)$weight >= corr.seq[[x]])) -> edges[[x]]  }
  names(edges) <- corr.seq
  
  ## Extract proper vertex names
  
  V(g)$colors <- as.character(unlist(color.list))
  
  values = as_data_frame(g, "vertices")
  values$hub_col = NA
  values$hub_col[values$name %in% hubs] = "gold"
  # values$hub_col[values$colors %in% dx.est$colors[dx.est$bonf < .05]] = "gold"
  
  V(g)$pie.color = list(c(as.character(net.grab), "gold"))
  
  val.list = list()
  for(n in 1:nrow(values)){
    cola = list(c(as.character(net.grab), "gold"))[[1]] == values$colors[[n]]
    cola = as.numeric(cola) 
    colb = list(c(as.character(net.grab), "gold"))[[1]] == values$hub_col[[n]]
    colb = as.numeric(colb) 
    colb[is.na(colb)] = 0
    val.list[[n]] = cola + colb
  }
  
  sum.take = unlist(lapply(val.list, sum))
  node.type = rep("circle", length(val.list))
  node.type[which(sum.take > 1)] = "pie"
  
  values$node.type = node.type
  
  V(g)$shape = node.type
  
  V(g)$pie.color[node.type == "circle"] = as.character(unlist(color.list)[node.type == "circle"])
  val.list[node.type == "circle"] = 1
  
  label.colors = rep("black", length(V(g)))
  label.colors[which(V(g)$colors == "black")] = "white"
  
  
  # gDel=igraph::delete.edges(g, which(E(g)$weight <= 0.1))
  gDel = g
  netdf = as_data_frame(minimum.spanning.tree(gDel, weights = E(g)$weights, algorithm = 'prim'), "edges")
  # netdf = as_data_frame(gDel)
  net = network(netdf[,1:2], directed = FALSE)
  
  # network plot
  png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood//ME12_network_minimumSpanningTree.png",res=300,units="in",height=6,width=6)
  ggnet2(net,
         vjust = 0,
         label.size = 2,
         mode = 'fruchtermanreingold',
         node.color = "tan",
         size = 2, 
         label = TRUE,
         edge.alpha = 0.5)
  dev.off()

  

  
  hubScore = hub_score(gDel)$vector

  hubs = sort(hubScore, decreasing=T)
  hubs = names(hubs[1:5])
  names(V(g))[names(V(g)) %in% hubs]
  
  size = rep(0.01, length(hubScore))
  size[which(names(V(g)) %in% hubs)] = 1.1
  
  png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood//ME12_network_minimumSpanningTree.png",res=600,units="in",height=6,width=6)
  par(mar=c(0,0,0,0)+.1)
  igraph::plot.igraph(mst(gDel),
                      vertex.size = 2, 
                      mark.groups = NULL,
                      vertex.shape= "circle",
                      vertex.color = 'tan',
                      vertex.label.cex= size,
                      vertex.label.color = label.colors,
                      edge.width=0.5,
                      vertex.label.font = 3,
                      vertex.frame.color='grey',
                      edge.color= "grey",
                      vertex.label.family="Arial")
  dev.off()
  
  
  
  
  
  # === next module (ME4)
  
  net.grab = "blue"
  
  slice.networks = list()
  color.list = list()
  for( y in 1:length(net.grab)){
    grab.genes =  which(moduleColors %in% net.grab[[y]])
    slice.networks[[y]] = datExpr[,grab.genes]
    color.list[[y]] = rep(net.grab[[y]], length(grab.genes))
  }
  slice.networks.df = as.data.frame(do.call(cbind, slice.networks))
  
  adjMat <- abs(cor(slice.networks.df, use = "p")) ## compute network-wide adjacency matrix 
  int_mod <- intramodularConnectivity(adjMat, colors = unlist(color.list)) ## computed gene connectivity from adjacency matrix 
  int_mod = data.frame(gene = rownames(int_mod), int_mod)
  mod.gene.col = data.frame(gene = colnames(datExpr), col = moduleColors)
  int_mod = merge(mod.gene.col, int_mod, by="gene")
  
  hubs = split(int_mod, as.character(int_mod$col))
  hubs = lapply(hubs, function(x) x[order(x$kWithin,decreasing=F), ])
  hubs = lapply(hubs, function(x) x[1,])
  hubs = as.character(do.call(rbind, hubs)$gene)
  
  color.list = which(colnames(datExpr) %in% hubs)
  color.list = moduleColors[color.list]
  
  # network graph of hub genes
  g <- graph.adjacency(as.matrix(adjMat),
                       mode="undirected",
                       weighted=T, diag = F) ## Graph information for plotting
  
  
  length(V(g)) -> vertices
  # g=delete.edges(g, which(E(g)$weight <=.2)) # filter edjges (defined by correlation coefficient)
  
  seq(.1, 1.0, by = .1) -> corr.seq
  
  edges <- list()
  for(x in 1:length(corr.seq)){
    x - 1 -> k
    length(which(E(g)$weight >= corr.seq[[x]])) -> edges[[x]]  }
  names(edges) <- corr.seq
  
  ## Extract proper vertex names
  
  V(g)$colors <- as.character(unlist(color.list))
  
  values = as_data_frame(g, "vertices")
  values$hub_col = NA
  values$hub_col[values$name %in% hubs] = "gold"
  # values$hub_col[values$colors %in% dx.est$colors[dx.est$bonf < .05]] = "gold"
  
  V(g)$pie.color = list(c(as.character(net.grab), "gold"))
  
  val.list = list()
  for(n in 1:nrow(values)){
    cola = list(c(as.character(net.grab), "gold"))[[1]] == values$colors[[n]]
    cola = as.numeric(cola) 
    colb = list(c(as.character(net.grab), "gold"))[[1]] == values$hub_col[[n]]
    colb = as.numeric(colb) 
    colb[is.na(colb)] = 0
    val.list[[n]] = cola + colb
  }
  
  sum.take = unlist(lapply(val.list, sum))
  node.type = rep("circle", length(val.list))
  node.type[which(sum.take > 1)] = "pie"
  
  values$node.type = node.type
  
  V(g)$shape = node.type
  
  V(g)$pie.color[node.type == "circle"] = as.character(unlist(color.list)[node.type == "circle"])
  val.list[node.type == "circle"] = 1
  
  label.colors = rep("black", length(V(g)))
  label.colors[which(V(g)$colors == "black")] = "white"
  
  
  gDel=igraph::delete.edges(g, which(E(g)$weight <= 0.1))
  # gDel = g
  netdf = as_data_frame(minimum.spanning.tree(gDel, weights = E(g)$weights, algorithm = 'prim'), "edges")
  # netdf = as_data_frame(gDel)
  net = network(netdf[,1:2], directed = FALSE)

  
  
  hubScore = hub_score(gDel)$vector
  
  hubs = sort(hubScore, decreasing = T)
  hubs = names(hubs)[1:5]
  names(V(g))[names(V(g)) %in% hubs]
  
  size = rep(0.001, length(hubScore))
  size[which(names(V(g)) %in% hubs)] = 1.1
  
  png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood//ME2_network_minimumSpanningTree.png",res=600,units="in",height=6,width=6)
  par(mar=c(0,0,0,0)+.1)
  igraph::plot.igraph(mst(gDel),
                      vertex.size = 1.5, 
                      mark.groups = NULL,
                      vertex.shape= "circle",
                      vertex.color = 'blue',
                      vertex.label.cex= size,
                      vertex.label.color = label.colors,
                      edge.width=0.5,
                      vertex.label.font = 3,
                      vertex.frame.color='grey',
                      edge.color= "grey",
                      vertex.label.family="Arial")
  dev.off()
  
  
  
  
  # === next module (ME20)
  
  net.grab = "yellow"
  
  slice.networks = list()
  color.list = list()
  for( y in 1:length(net.grab)){
    grab.genes =  which(moduleColors %in% net.grab[[y]])
    slice.networks[[y]] = datExpr[,grab.genes]
    color.list[[y]] = rep(net.grab[[y]], length(grab.genes))
  }
  slice.networks.df = as.data.frame(do.call(cbind, slice.networks))
  
  adjMat <- abs(cor(slice.networks.df, use = "p")) ## compute network-wide adjacency matrix 
  int_mod <- intramodularConnectivity(adjMat, colors = unlist(color.list)) ## computed gene connectivity from adjacency matrix 
  int_mod = data.frame(gene = rownames(int_mod), int_mod)
  mod.gene.col = data.frame(gene = colnames(datExpr), col = moduleColors)
  int_mod = merge(mod.gene.col, int_mod, by="gene")
  
  hubs = split(int_mod, as.character(int_mod$col))
  hubs = lapply(hubs, function(x) x[order(x$kWithin,decreasing=F), ])
  hubs = lapply(hubs, function(x) x[1,])
  hubs = as.character(do.call(rbind, hubs)$gene)
  
  color.list = which(colnames(datExpr) %in% hubs)
  color.list = moduleColors[color.list]
  
  # network graph of hub genes
  g <- graph.adjacency(as.matrix(adjMat),
                       mode="undirected",
                       weighted=T, diag = F) ## Graph information for plotting
  
  
  length(V(g)) -> vertices
  # g=delete.edges(g, which(E(g)$weight <=.2)) # filter edjges (defined by correlation coefficient)
  
  seq(.1, 1.0, by = .1) -> corr.seq
  
  edges <- list()
  for(x in 1:length(corr.seq)){
    x - 1 -> k
    length(which(E(g)$weight >= corr.seq[[x]])) -> edges[[x]]  }
  names(edges) <- corr.seq
  
  ## Extract proper vertex names
  
  V(g)$colors <- as.character(unlist(color.list))
  
  values = as_data_frame(g, "vertices")
  values$hub_col = NA
  values$hub_col[values$name %in% hubs] = "gold"
  # values$hub_col[values$colors %in% dx.est$colors[dx.est$bonf < .05]] = "gold"
  
  V(g)$pie.color = list(c(as.character(net.grab), "gold"))
  
  val.list = list()
  for(n in 1:nrow(values)){
    cola = list(c(as.character(net.grab), "gold"))[[1]] == values$colors[[n]]
    cola = as.numeric(cola) 
    colb = list(c(as.character(net.grab), "gold"))[[1]] == values$hub_col[[n]]
    colb = as.numeric(colb) 
    colb[is.na(colb)] = 0
    val.list[[n]] = cola + colb
  }
  
  sum.take = unlist(lapply(val.list, sum))
  node.type = rep("circle", length(val.list))
  node.type[which(sum.take > 1)] = "pie"
  
  values$node.type = node.type
  
  V(g)$shape = node.type
  
  V(g)$pie.color[node.type == "circle"] = as.character(unlist(color.list)[node.type == "circle"])
  val.list[node.type == "circle"] = 1
  
  label.colors = rep("black", length(V(g)))
  label.colors[which(V(g)$colors == "black")] = "white"
  
  
  gDel=igraph::delete.edges(g, which(E(g)$weight <= 0.1))
  # gDel = g
  netdf = as_data_frame(minimum.spanning.tree(gDel, weights = E(g)$weights, algorithm = 'prim'), "edges")
  # netdf = as_data_frame(gDel)
  net = network(netdf[,1:2], directed = FALSE)
  
  
  
  hubScore = hub_score(gDel)$vector
  
  hubs = sort(hubScore, decreasing = T)
  hubs = names(hubs)[1:5]
  names(V(g))[names(V(g)) %in% hubs]
  
  size = rep(0.001, length(hubScore))
  size[which(names(V(g)) %in% hubs)] = 0.9
  
  png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood//ME4_network_minimumSpanningTree.png",res=600,units="in",height=6,width=6)
  par(mar=c(0,0,0,0)+.1)
  igraph::plot.igraph(mst(gDel),
                      vertex.size = 2, 
                      mark.groups = NULL,
                      vertex.shape= "circle",
                      vertex.color = 'yellow',
                      vertex.label.cex= size,
                      vertex.label.color = label.colors,
                      edge.width=0.5,
                      vertex.label.font = 3,
                      vertex.frame.color='grey',
                      edge.color= "grey",
                      vertex.label.family="Arial")
  dev.off()
  
  
  
  
  
  ## Gene set enrichment analysis for module 
  
  
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  require(org.Hs.eg.db)
  #
  ## Perform gene set enrichment analysis ...
  msig = list.files(path = "~/Google Drive/mac_storage/gene_sets/", pattern = "msig.txt", full.names = T)
  msig = msig[!grepl("c1.all",msig)]
  msig = lapply(msig, function(x) fread(x, h = T))
  msig = rbindlist(msig)
  msig = msig[grepl("KEGG|GO|REACTOME|PANTHER",msig$PATHWAY), ]
  msig = msig[!grepl("_DN|_UP",msig$PATHWAY)]
  
  
  msig$SLIM = NA
  msig$SLIM[grepl("metabol|peptidase|lase",ignore.case = T, msig$PATHWAY)]="Metabolism and enzymatic activity"
  msig$SLIM[grepl("immun|inflammato|lupus|humoral|innnate|defense|leukocyte|neutrophil|platelet|B cell|lymphoblast|monocyte|blood|inflammation|infereon|cytokine|nf kappa",ignore.case = T, msig$PATHWAY)] = "Immunity"
  msig$SLIM[grepl("cell death|apoptosis|apotosi|survival|autophagy|growth|repair|division",ignore.case = T, msig$Set_name)] = "Cell growth, repair, and survival"
  msig$SLIM[grepl("chromatin|histone|epigeno|methyltransferase|acetyltransferase",ignore.case = T, msig$PATHWAY)] = "Chromatin regulation"
  msig$SLIM[grepl("cytoske|part|lysosome|endosom|golgi|plasma|mitochondri|ribos|nucleu|secreto|vesicl|extracellu",ignore.case = T, msig$PATHWAY)] = "(Extra)cellular structures and organelles"
  msig$SLIM[grepl("fungu|bacteri|external|yeast|symbio|stimulus|microbia|parasit|other",ignore.case = T, msig$PATHWAY)] = "Response to external stimuli"
  msig$SLIM[grepl("neurotrans|cephalon|myelin|dopamine|gliogenesis|neurogenesis|astrogenesis|glutama|gaba|calcium|synapse|neuro|glial|astrocyte|axon|brain|behavior|cognitive",ignore.case = T, msig$PATHWAY)] = "Neuronal gene sets"
  msig$SLIM[grepl("adhesion|motility|movement|integrin|cell-cell",ignore.case = T, msig$PATHWAY)] = "Cell adhesion and contact"
  msig$SLIM[is.na(msig$SLIM)] = "Other processes"
  
  
  # pathCounts = table(msig$PATHWAY)
  # pathCounts = pathCounts[pathCounts > 5 & pathCounts < 300]
  
  msig = msig[msig$PATHWAY %in% names(pathCounts), ]
  
  sig_mods = as.character(unique(graphDF$colors[graphDF$FDRALL < 0.05]))
  
  # remove number of characters from long pathway names
  ncharcount = nchar(msig$PATHWAY)
  toomanychar = which(ncharcount > 250)
  
  substr = substr(x = msig$PATHWAY[toomanychar], 1, 55)
  substr = paste(substr, "...",sep="")
  msig$PATHWAY[toomanychar] = substr
  
  
  # number of top ranked gene sets to report per SLIM set
  topValue = 10
  
  hypStats_save = list()
  for( i in 1:length(sig_mods)){
    
    cat("\nRunning hypergeometric test for module:",sig_mods[[i]])
    
    gene_grab = module[module$color %in% sig_mods[[i]], ]
    
    msig = msig[msig$SYMBOL %in% module$symbol]
    msig = msig[grepl("chr|KEGG_|GO_|REACTOME|PANTHER", msig$PATHWAY)]
    msig = msig[!is.na(msig$SYMBOL)]
    
    
    # filter sets
    set_count = table(as.character(msig$PATHWAY))
    
    list = unique(names(set_count))
    
    msig.keep = msig[msig$PATHWAY %in% list]
    
    msig.split = split(msig, msig$PATHWAY)
    
    # calculate hypergeometric input values
    universe = length(unique(module$symbol)) # total # of genes in transcriptome 
    overlaps = lapply(msig.split, function(x) length(intersect(x$SYMBOL, gene_grab$symbol)))
    geneSetSize = lapply(msig.split, nrow)
    geneSetOut = universe - unlist(geneSetSize)
    
    # hypergeomtric p-values
    hypStats = data.frame(
      Module = sig_mods[[i]],
      PATHWAY = names(msig.split),
      Population = universe, 
      Sample_success = unlist(overlaps),
      Population_success = unlist(geneSetSize),
      Population_fails = unlist(geneSetOut),
      Sample_size = length(unique(gene_grab$symbol))
    )
    
    hypStats = hypStats[hypStats$Sample_success >= 2, ]
    
    # enrichment p-value test
    pvals = phyper(hypStats$Sample_success - 1, 
                   hypStats$Population_success, 
                   hypStats$Population_fails, 
                   hypStats$Sample_size, 
                   lower.tail = FALSE, log.p = FALSE)
    
    hypStats$P = pvals
    
    hypStats = hypStats[order(hypStats$P, decreasing = F), ]
    hypStats$FDR = p.adjust(hypStats$P, "fdr")
    hypStats$BONF = p.adjust(hypStats$P, "bonferroni")
    
    hypStats = merge(hypStats, msig[,c("SLIM","PATHWAY")][!duplicated(msig$PATHWAY),], by = "PATHWAY")
    
    SLIMSETS = unique(hypStats$SLIM)
    tmp_list= list()
    for(x in 1:length(SLIMSETS)){
      
      tmp = hypStats[hypStats$SLIM %in% SLIMSETS[[x]], ]
      tmp = tmp[order(tmp$P, decreasing = F), ]
      if(topValue == 0){topValue = nrow(hypStats)}
      tmp = tmp[1:topValue, ]
      tmp_list[[x]] = tmp
      
    }
    hypStats = ldply(tmp_list)
    names(hypStats)[names(hypStats) %in% "PATHWAY"] = "SET_NAME"
    hypStats_save[[i]] = hypStats # save sum stats to list obj 
    
  }
  
  hypStats_save = ldply(hypStats_save)
  
  
  sets = strsplit(as.character(hypStats_save$SET_NAME), "_")
  pfx = lapply(sets, function(x) x[[1]])
  suffix = lapply(sets, function(x) x[-1])
  
  hypStats_save$relabel = NA
  for( i in 1:length(sets)){
    sfx = paste(unlist(suffix[[i]]), collapse = " ", sep = "")
    hypStats_save$relabel[[i]] = paste(pfx[[i]], ":", sfx, sep = "")
  }
  
  
  hypStats_save$p_label = paste("P = ", format(hypStats_save$P,scientific=T,digits=3), sep = "")
  
  labels = data.frame(Label = moduleLabels, Module = moduleColors)
  labels$Label = paste("ME", labels$Label,sep="")
  labels = labels[!duplicated(labels$Module), ]
  
  hypStats_save = hypStats_save[order(hypStats_save$P,decreasing = F), ]
  hypStats_save = hypStats_save[!duplicated(hypStats_save$SET_NAME),]
  
  
  hypStats_save$foldEnrichment = (hypStats_save$Sample_success/hypStats_save$Sample_size)/(hypStats_save$Population_success/hypStats_save$Population)
  
  hypStats_save = merge(hypStats_save, labels, by='Module')
  
  hypStats_save = hypStats_save[order(hypStats_save$SLIM, hypStats_save$P, decreasing = F), ]
  
  # remove categories with no significance
  hypStats_save = hypStats_save[hypStats_save$FDR < .05, ]
  
  threshold = -log10(.05/length(unique(msig$PATHWAY)))
  
  hypStats_save$relabel = gsub("GO:|KEGG:|Reactome:", "", hypStats_save$relabel)
  hypStats_save$relabel = stringr::str_to_title(hypStats_save$relabel)
  
  
  split_set = split(hypStats_save, hypStats_save$Module)
  split_set = lapply(split_set, function(x) x[order(x$P,decreasing = F), ])
  split_set = lapply(split_set, function(x) x[1:5, ])
  top_sets = data.frame(do.call(rbind, split_set))
  top_sets$LOGP = -log10(top_sets$P)
  
  require(ggplot2)
  
  
  ySize = rep(0.5, 35)
  ySet = seq(from = ceiling(5/2), length.out = length(unique(top_sets$Label)))
  ySize[ySet] = 1e-5
  
  png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/GSEA_barplot.png",res=300,units = "in",height=5,width=6.5)
  par(mar=c(3,3,0,1)+1)
  barplot(top_sets$LOGP,
          xlim=c(0,80),
          xlab = expression(paste("Gene set enrichment, -log"[10],"(P-value)")),
          col = top_sets$Module, las = 1, cex.names = 0.5,
          names.arg = top_sets$Label, horiz= T)
  abline(v = -log10(.05), col = 'black', lty = 1, lwd = .5)
  
  col = ifelse(grepl("black|purple|dark", ignore.case = T, top_sets$Module ) == T, "grey", "black")
  text(x = 1, col = col, cex = 0.5, pos = 4, y = seq(from = 0.7, by = 1.2, length.out = nrow(top_sets)), 
       labels = top_sets$relabel)
  dev.off()
  
  
  # semantic similarity between modules for FDRp < 0.05 enriched gene sets
  split_set = split(hypStats_save, hypStats_save$Label)
  split_set = lapply(split_set, function(x) x$relabel)
  
  sets = expand.grid(names(split_set), names(split_set))
  sim = list()
  for(x in 1:nrow(sets)){
    a = sets$Var1[[x]]
    b = sets$Var2[[x]]
  calc = stringdist(method = "lcs", a = paste(split_set[names(split_set) %in% a], collapse = " "), b = paste(split_set[names(split_set) %in% b], collapse = ", "))
  sets$val[[x]] = calc
  }
  sets$std = (sets$val - max(sets$val))/(min(sets$val) - max(sets$val))
  sim = reshape2::acast(sets, Var1 ~ Var2, value.var = 'std')
  
  png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/similarity_text_gsea.png",res=300,units="in",height=5,width=5.5)
  pheatmap(sim, border_color = 'black', lwd = 1.75)
  dev.off()
  
  require(Biostrings)
  # calculate distances between strings
  
  colors = unique(hypStats_save$Module)
  
  require(pheatmap)
  pdf("~/Google Drive/mac_storage/TWAS/bd_mega/MajorPsychosis_similarity_genesets.pdf",width = 12, height =8)
  par(mfrow=c(1,1))
  for(y in 1:length(colors)){
    
  grab = hypStats_save[hypStats_save$Module %in% colors[[y]], ]
  strings = expand.grid(grab$relabel, grab$relabel)
  
  strings = strings[!duplicated(strings), ]
  strings$Var1 = as.character(strings$Var1)
  strings$Var2 = as.character(strings$Var2)
  for(x in 1:nrow(strings)){
    cat("\rCalculating string distance:",100*round(x/nrow(strings),3),"%")
    strings$dist[[x]] = stringDist(c(strings$Var1[[x]], strings$Var2[[x]]))[[1]]
    s1 = strsplit(strings$Var1[[x]], " ")[[1]]
    s2 = strsplit(strings$Var2[[x]]," ")[[1]]
    strings$similarity[[x]] = length(intersect(s1, s2))/(length(union(s1,s2)))
  }
  string_dist = reshape2::acast(strings, Var1 ~ Var2, value.var = 'similarity')
  
  pheatmap(string_dist, main = paste("Module:", grab$Label[[y]]), cluster_rows = T, fontsize_row = 8, fontsize_col = 5,
           cluster_cols = T, show_rownames = TRUE, show_colnames = FALSE)
  
  }
  dev.off()
  
  
  write.csv(hypStats_save,
            file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/wgnca-psychosis-gsea.csv",
            quote= F, row.names= F)
  
  png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/wgcna-psychosis_gsea.png",res=300,units="in",height=6,width=8)
  
  ggplot(hypStats_save, aes(x = relabel, y = -log10(P))) + 
    xlab(NULL) +
    ylab(expression(paste("-log"[10],"(P-value)"))) +
    facet_grid(SLIM ~., scales = 'free_y') + 
    geom_bar(stat = "identity", col = "black",aes(fill = I(hypStats_save$Module))) +
    geom_hline(yintercept = threshold, col = "black", linetype = 2, lwd = 0.5) + 
    theme_classic() + 
    ylim(min = 0, max = max(-log10(hypStats_save$P))*1.75) +
    geom_text(aes(label = p_label, y = max(-log10(P))*1.5), size= 2) +
    theme(axis.text.y = element_text(size = 5), strip.text.y = element_text(size = 8, angle = 0)) +
    coord_flip()
  
  dev.off()
  
  
  
  
  
  
  
  
# find columns that are in both SZ and BD data sets
int = intersect(colnames(szdata)[-c(1:12)], colnames(mydata)[-c(1:11)])
# int = intersect(colnames(szdata), colnames(mydata))


## combine using intersected columns 
datExpr0 = mydata
datExpr1 = szdata

datExprList = list(datExpr0, datExpr1)


GS0 = list()
for( i in 1:length(datExprList)){
  
  cat("\nDifferential expression of models:",i,"\n")
  datExpr = as.data.frame(datExprList[[i]])

  ## Gene significance
  phenos = datList[[i]][,colnames(datList[[i]]) %in% c("dx","studyID", "age","sex"),with=F]
  phenos$sex[phenos$sex %in% ""] = NA
  phenos$studyID = factor(phenos$studyID)
  phenos$dx[grepl("contro|CT", phenos$dx)] = "CT"
  
  phenos$dx = relevel(as.factor(phenos$dx), ref = "CT")
  
  design = model.matrix( ~., data =  phenos)
  
  response = datExpr[,15:ncol(datExpr)]
  response = response[,!colSums(is.na(response))>0]
  
  response = response[rownames(response) %in% rownames(design), ]
  
  GSfit = lm(as.matrix(response) ~ -1 + design)
  
  GScoef = summary(GSfit)
  GScoef = lapply(GScoef,function(x) broom::tidy(x$coefficients))
  names(GScoef) = gsub("Response ", "", names(GScoef))
  GSdf = ldply(GScoef)
  GSdf = GSdf[grepl("dx",GSdf$.rownames),]
  colnames(GSdf) = c("symbol","Term","Beta","SE","T","Pval")
  GSdf = merge(module, GSdf,by="symbol")
  GS0[[i]] = GSdf

}

lapply(dem_stats,head)

GSset = ldply(GS0)

# calculate mean T-value 

geneSig = data.table(GSset)
geneSig = geneSig[,list(Tvalue = mean(T)),by=list(symbol)]

require(org.Hs.eg.db)
require(AnnotationDbi)

genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes = data.frame(genes)
genes$symbol = select(org.Hs.eg.db, keys=as.character(genes$gene_id), keytype="ENTREZID",columns="SYMBOL")$SYMBOL

geneSig = merge(genes, geneSig, by='symbol')

fwrite(geneSig, file = "~/Documents/FLEET/psychosis.txt",
       quote= F, row.names = F, sep = "\t")

# correlation of SCZ v. BD gene significance for top modules

sig.module = GSset[GSset$color %in% "steelblue",]

cast = reshape2::acast(sig.module, symbol ~ Term, value.var="T")
cast = as.data.frame(cast)

ggplot(cast, aes(x = cast[,1],  y= cast[,2])) + 
  geom_point() +
  theme_classic() + 
  xlab("BD gene significance") +
  ylab("SZ gene significance") + 
  geom_smooth(method="lm")









# List-wise approach to WGCNA 
datExprList = list()
for(i in 1:length(datList)){
  
  cat("\n***** Detecting usuable genes and samples on data set:",i,"*****\n")
expr = as.data.frame(datList[[i]])
  
## WGCNA 
datExpr = expr[,colnames(expr) %in% int]
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

keepSamples = expr[,c(1:11)]
keepSamples = keepSamples[which(gsg$goodSamples == TRUE),]

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK


datExprList[[i]] = datExpr

}

lapply(datExprList, dim)

keepGenes = lapply(datExprList, function(x) unique(colnames(x)))
keepGenes = Reduce(function(x,y) intersect(x,y), keepGenes)
length(keepGenes)

net_list = list()
for( i in 1:length(datExprList)){
  
  cat("\n****** Building WGCNA modules for data set:",i,"******\n")
  
  phenos = as.data.frame(datList[[i]])[,c(1:11)]
  datExpr = datExprList[[i]]
  datExpr = datExpr[,colnames(datExpr) %in% keepGenes]
  datExpr = datExpr[,order(colnames(datExpr))]
  

powers = c(c(1:10), seq(from = 12, to=20, by=2))

iters=10
sft1 = list()
for( y in 1:iters){
  random = datExpr[,sample(colnames(datExpr), 500, replace = FALSE)]
  sft1[[y]] = pickSoftThreshold(random, 
                                powerVector = powers, 
                                corFnc="bicor",
                                networkType="signed")
}
sft1_power = lapply(sft1, function(x) x[[1]])
topPower = names(sort(table(unlist(sft1_power)),decreasing = T))
sft1_in = lapply(sft1, function(x) x$fitIndices)
names(sft1_in) = iters
sft1_df = ldply(sft1_in)
sft1_avg = aggregate(sft1_df$SFT.R.sq, by = list(sft1_df$Power), mean)
sft1_se = aggregate(sft1_df$SFT.R.sq, by = list(sft1_df$Power), function(x) sd(x)/sqrt(length(x)))
sft1_avg$CI_LOWER = sft1_avg$x - ( 1.96 * sft1_se$x)
sft1_avg$CI_HIGH = sft1_avg$x + ( 1.96 * sft1_se$x)
names(sft1_avg)[names(sft1_avg) %in% "x"] = "mean"
sft1_avg$label = ifelse(sft1_avg$Group.1 == topPower[[i]], "*", NA)

datGroups = unique(phenos$dx)
datGroups = paste(datGroups, collapse = ", ",sep ="")

r2 = ggplot(sft1_avg, aes(x = factor(Group.1), y = mean)) +
  geom_bar(stat = 'identity', fill = heat.colors(15), colour = 'black') + 
  theme_classic() +
  xlab("Soft threshold power") +
  ylab(expression(paste("Model fit ",italic(R)^2))) +
  geom_errorbar(aes(ymin = CI_LOWER, ymax = CI_HIGH), width=0.2) +
  geom_hline(aes(yintercept = .8),col='navyblue',linetype='dashed') + 
  ggtitle(paste("Groups = ",datGroups,sep="")) +
  geom_text(aes(y = CI_HIGH, label = label),nudge_y = 0.05, size = 10)


png(paste("~/Google Drive/mac_storage/TWAS/WGCNA_network_set-",i,".png",sep=""),res=300,units="in",height=6,width=6)
print(r2)
dev.off()


## WGCNA network algorithm

# 3. set parameters for network algorithm
sft.power = as.integer(topPower[[1]]); ## SFT 8, R-squared = 0.88013307 (SE = 0.05562321)
deepSplit = 2;
minModuleSize = 30;

# 4.one-step automated network construction
net_list[[i]] = blockwiseModules(datExpr, 
                        power = sft.power, 
                        networkType = "signed",
                        deepSplit= deepSplit, 
                        TOMType = "signed", 
                        minModuleSize = minModuleSize, 
                        minCoreKME = 0.5, 
                        minCoreKMESize = minModuleSize/3,
                        minKMEtoStay=0, 
                        reassignThreshold = 1e-6, 
                        mergeCutHeight = 0.25, 
                        detectCutHeight = 0.995,
                        numericLabels = TRUE, 
                        pamRespectsDendro = FALSE, 
                        pamStage = TRUE,
                        saveTOMs = TRUE, 
                        verbose = 3,
                        maxBlockSize = 5000)

}


# relabel modules in netB based on module membership in netA
datSource =  datExprList[[1]]
datSource = datSource[,colnames(datSource) %in% keepGenes]
datRef = datExprList[[2]]
datRef = datRef[,colnames(datRef) %in% keepGenes]

net_match = matchLabels(source = net_list[[1]]$colors, 
                        reference = net_list[[2]]$colors,
                        pThreshold = 5e-2,
                        na.rm = TRUE)


# 6. extract network meta-data and eigengenes


# find columns that are in both SZ and BD data sets
int = intersect(colnames(szdata)[-c(1:12)], colnames(mydata)[-c(1:11)])
# int = intersect(colnames(szdata), colnames(mydata))


## combine using intersected columns 
datExpr0 = mydata
datExpr1 = szdata

datExprList = list(datExpr0, datExpr1)


moduleLabels=list()
moduleColors = list()
MEs = list()
mergedColors=list()
dem_stats = list()
GS0 = list()
for( i in 1:length(net_list)){
  
  cat("\nDifferential expression of models:",i,"\n")
  datExpr = as.data.frame(datExprList[[i]])
  datExpr = datExpr[,colnames(datExpr) %in% keepGenes]
  datExpr = datExpr[,order(colnames(datExpr))]
  
  
  
  ## Gene significance
  phenos = datList[[i]][,colnames(datList[[i]]) %in% c("dx","studyID","age","sex"),with=F]
  phenos$dx[grepl("contro|CT", phenos$dx)] = "CT"
  phenos$dx = relevel(as.factor(phenos$dx), ref = "CT")
  design = model.matrix( ~., phenos)
  response = datExpr
  response = response[rownames(response) %in% rownames(design), ]
  GSfit = lm(as.matrix(datExpr) ~ design)
  GScoef = summary(GSfit)
  GScoef = lapply(GScoef,function(x) broom::tidy(x$coefficients))
  names(GScoef) = colnames(datExpr)
  GSdf = ldply(GScoef)
  GSdf = GSdf[grepl("dx",GSdf$.rownames),]
  colnames(GSdf) = c("GeneSymbol","Term","Beta","SE","T","Pval")
  GSdf$Comparison = paste(unique(phenos$dx),collapse=",",sep="")
  GSdf = merge(geneColor, GSdf,by="GeneSymbol")
  GS0[[i]] = GSdf
  
  # statistical analysis 
  
  phenos = as.data.frame(datList[[i]])[,c(1:12)]
  phenos$dx[grepl("control|CT", phenos$dx)] = "CT"
  phenos$dx = relevel(as.factor(phenos$dx), ref = "CT")
  
  # rownames match
  MEs = MEs[match(rownames(phenos), rownames(MEs)), ]
  datAll = data.frame(phenos, MEs)
  
  pcaCov = pcaCov[match(rownames(datAll), rownames(pcaCov)), ]
  datAll = data.frame(pcaCov, datAll)
  
  lmFitCellMix = lm(as.matrix(datAll[,colnames(datAll) %in% colnames(pcaCov)]) ~ datAll$dx)
  lmCoef = summary(lmFitCellMix)
  lmCoef = lapply(lmCoef, function(x) broom::tidy(x$coefficients))
  names(lmCoef) = colnames(pcaCov)[colnames(pcaCov) %in% colnames(datAll)]
  lmdf = ldply(lmCoef)
  lmCoef = lmdf[grepl("dx",lmdf$.rownames),]
  
  pickCellMixCovs = lmCoef$.id[lmCoef$Pr...t.. < 5e-02] # Keep cellmix PCs that are significantly different between disease groups
  
  # Linear models to test for differetial eigengene expression between groups
  deplist = colnames(MEs)
  

  design = model.matrix( ~ ., datAll[,colnames(datAll) %in% c("sex","dx", "age", "studyID",pickCellMixCovs)])
  y_mat = datAll[,colnames(datAll) %in% deplist]
  y_mat = y_mat[rownames(y_mat) %in% rownames(design), ]
  lmFit = lm(as.matrix(y_mat) ~ -1 + design)
  lmCoef = summary(lmFit)
  lmCoef = lapply(lmCoef, function(x) broom::tidy(x$coefficients))
  names(lmCoef) = colnames(y_mat)
  lmDf = ldply(lmCoef)
  lmDf = lmDf[grepl("designdx", lmDf$.rownames), ]
  lmDf = lmDf[order(lmDf$Pr...t..,decreasing=F), ]
  lmDf = lmDf[!lmDf$.id %in% "ME0",]
  lmDf$FDR = p.adjust(lmDf$Pr...t.., "fdr")
  lmDf$BONF = p.adjust(lmDf$Pr...t..,"bonferroni")
  lmDf$Comparisons = paste(unique(phenos$dx),collapse=",",sep="")
  head(lmDf)
  
  dem_stats[[i]] = lmDf
}

lapply(dem_stats,head)

GSset = ldply(GS0)

statsDf = ldply(dem_stats)
statsDf = statsDf[order(statsDf$.id),]

head(statsDf)

deplist = lapply(dem_stats,function(x) x$.id)
deplist = Reduce(function(x,y)intersect(x,y), deplist)

metaStats = list()
corStat = list()
for( i in 1:length(deplist)){
  cat("\rComparing WGCNA networks:",i)
  # Meta-analyze module
  # tmpstats = statsDf[statsDf$.id %in% deplist[[i]], ]
  # meta = metafor::rma(yi = tmpstats$Estimate, tmpstats$Std..Error, weighted = F)
  # metaB  = meta$b[[1]]
  # metaSE = meta$se
  # metaP = meta$pval
  # metaStats[[i]] = data.frame(.id = deplist[[i]], Estimate = metaB, Std..Error = metaSE, Pr...t.. = metaP, Comparisons = "Psychosis,CT")
  # subset genes from a module, correlation of gene significance
  
  ModuleSub = GSset[GSset$ME %in% deplist[[i]], ]
  frame = reshape2::acast(ModuleSub, GeneSymbol~Comparison, value.var="T")
  
  plot_df = frame
  colnames(plot_df) = toupper(unlist(lapply(strsplit(colnames(plot_df), ","),function(x)x[[1]])))
  plot_df = as.data.frame(plot_df)
  
  
  corStat = data.frame(ME = deplist[[i]], broom::tidy(cor.test(frame[,1], frame[,2])))
  
  g = ggplot(plot_df, aes(x=plot_df[,1], y= plot_df[,2])) + 
    geom_point() +
    theme_classic() +
    xlab(paste("Module ", deplist[[i]], " (",colnames(plot_df)[1]," gene significance)", sep="")) +
    ylab(paste("Module ",deplist[[i]], " (",colnames(plot_df)[2]," gene significance)", sep="")) +
    geom_smooth(method='lm', lwd = 0.5) +
    ggtitle(paste("r = ", round(corStat$estimate,3), ", p = ", format(corStat$p.value,digits=3,scientific=T), sep="" ))
  
  png(paste("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/plots/GScor_",deplist[[i]],".png",sep=""), res=300,units="in",height=6,width=6)
  print(g)
  dev.off()
  
  
  
  metaStats[[i]] = data.frame(.id = deplist[[i]], Estimate = corStat$estimate,  Pr...t.. = corStat$p.value, Comparisons = "BD+SZ,CT")
  
  
}
metaStats_df = ldply(metaStats)
metaStats_df = metaStats_df[order(metaStats_df$Pr...t..,decreasing=F), ]
metaStats_df$FDR = p.adjust(metaStats_df$Pr...t..,"fdr")
metaStats_df$BONF = p.adjust(metaStats_df$Pr...t.., "bonferroni")

head(metaStats_df)

comb = ldply(list(statsDf, metaStats_df))
comb = comb[comb$.id %in% deplist, ]
comb = comb[!is.na(comb$.id) | !comb$.id %in% NA, ]
comb$FDR = p.adjust(comb$Pr...t..,'fdr')
comb$BONF = p.adjust(comb$Pr...t..,"bonferroni")
# comb_split = split(comb, comb$Comparisons)
# for(x in 1:length(comb_split)){
#   comb_split[[x]]$FDR = p.adjust(comb_split[[x]]$Pr...t..,"fdr")
#   comb_split[[x]]$BONF = p.adjust(comb_split[[x]]$Pr...t..,"bonferroni")
# }

# comb$FDR = p.adjust(comb$Pr...t.., "fdr")
# comb$BONF = p.adjust(comb$Pr...t.., "bonferroni")
comb$.id = factor(comb$.id, levels=paste("ME",1:100,sep=""))
comb$Comparisons = gsub("[,]", " v. ", toupper(comb$Comparisons))
comb$Comparisons = factor(comb$Comparisons, levels=c("BD v. CT", "SCZ v. CT", "BD+SZ v. CT"))
comb$Label = ifelse(comb$FDR < .05, "*", NA)

g = ggplot(comb, aes(x = .id, y = -log10(Pr...t..))) +
  geom_bar(stat='identity', col = 'black', width=0.75, position='dodge', aes(fill = Comparisons)) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1, angle = 90), legend.position = "bottom") +
  xlab(NULL) +
  ylab(expression(paste("-log"[10],"(P-value)"))) +
  geom_hline(aes(yintercept = -log10(.05)), col = 'red', linetype=  'dashed')  +
  labs(fill = NULL, col = c("P < 0.05"))  +
  scale_fill_brewer(palette="Paired") + 
  geom_text(aes(label = Label, group = Comparisons),  position = position_dodge(width = 0.8), size = 10)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/plots/WGCNA_combined.png",res=300,units="in",height=6,width=6)
print(g)
dev.off()


# biological annotation of significant modules 










# 
# 
# # table(mdsphenos$dx[rownames(mdsphenos) %in% rownames(datExpr)])
# 
# 
# disableWGCNAThreads() # number of processors reported by system will be used 
# 
# ## one-step automated gene network analysis 
# # 1. find optimal soft-threshold power for network construction
# 
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# 
# iters=10
# sft1 = list()
# for( y in 1:iters){
#   random = datExpr[,sample(colnames(datExpr), 1000, replace = FALSE)]
#   sft1[[y]] = pickSoftThreshold(random, 
#                                 powerVector = powers, 
#                                 corFnc="bicor",
#                                 networkType="signed")
# }
# sft1_power = lapply(sft1, function(x) x[[1]])
# sort(table(unlist(sft1_power)),decreasing = T)
# sft1_in = lapply(sft1, function(x) x$fitIndices)
# names(sft1_in) = iters
# sft1_df = ldply(sft1_in)
# sft1_avg = aggregate(sft1_df$SFT.R.sq, by = list(sft1_df$Power), mean)
# sft1_se = aggregate(sft1_df$SFT.R.sq, by = list(sft1_df$Power), function(x) sd(x)/sqrt(length(x)))
# sft1_avg$CI_LOWER = sft1_avg$x - ( 1.96 * sft1_se$x)
# sft1_avg$CI_HIGH = sft1_avg$x + ( 1.96 * sft1_se$x)
# names(sft1_avg)[names(sft1_avg) %in% "x"] = "mean"
# 
# r2 = ggplot(sft1_avg, aes(x = factor(Group.1), y = mean)) +
#   geom_bar(stat = 'identity', fill = heat.colors(15), colour = 'black') + 
#   theme_classic() +
#   xlab("Soft threshold power") +
#   ylab(expression(paste("Model fit ",italic(R)^2))) +
#   geom_errorbar(aes(ymin = CI_LOWER, ymax = CI_HIGH), width=0.2) +
#   geom_hline(aes(yintercept = .8),col='navyblue',linetype='dashed')
# 
# print(r2)
# 
# png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/SFT_WGCNA_BDvSZvCT-signednetwork_permuted10x.png",res=300,units="in",height = 6, width = 7)
# print(r2)
# dev.off()


# 
# # start expression data at colID
# start_gene = 7
# 
# # estimate circulating leukocyte proportions from gene expression signatures
# require(CellMix)
# require(hgu133a.db)
# 
# exprs = as.data.frame(datExpr[,start_gene:ncol(datExpr)]) # extract normalized gene expression intensities for subjects
# exprs = exprs[,!colSums(is.na(exprs))]
# 
# convert = AnnotationDbi::select(hgu133a.db, keys = colnames(exprs), keytype="SYMBOL", columns="PROBEID")
# 
# minval = min(apply(exprs, 2, function(x) min(x, na.rm = T))) # find minimum value of full matrix
# exprs = exprs + abs(minval) # add |minimum value| to achieve non-negative matrix
# 
# 
# exprs = t(exprs)
# convert = convert[match(rownames(exprs), convert$SYMBOL), ]
# rownames(exprs) = convert$PROBEID
# exprs = exprs[!is.na(rownames(exprs)), ] # remove missing rownames
# exprs = exprs[!duplicated(rownames(exprs)), ]
# 
# exprs = ExpressionSet(exprs) # convert to ExpressionSet object
# 
# 
# res = gedBlood(exprs, verbose = T, normalize = T) # run gedBlood algorithm (non-negative matrix factorization)
# wb.coef = coef(res)
# wb.coef = as.data.frame(t(wb.coef))
# # missing cols
# # miss.col = colSums(wb.coef == 0)/nrow(wb.coef)
# # miss.col = miss.col[miss.col > 0.5]
# # wb.coef = wb.coef[,!colnames(wb.coef) %in% names(miss.col)]
# 
# wb.results = data.frame(sampleIDs = datExpr$sampleIDs, 
#                         dx = datExpr$dx, 
#                         wb.coef)
# 
# 
# ## MDS analysis of leukocytes
# mds = as.data.frame(cmdscale(dist(wb.coef)))
# mds = as.data.frame(prcomp(wb.coef, scale = T, center = T)$x)
# 
# par(mfrow=c(1,1))
# 
# ggplot(mds, aes(x = PC1, y = PC2, col = factor(wb.results$dx))) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_color_discrete('Diagnosis') + 
#   xlab('PC 1') + 
#   ylab('PC 2')
# 
# mdsphenos = data.frame(wb.results, mds)
# mdsphenos$dx = as.character(mdsphenos$dx)
# 
# fwrite(setDT(mdsphenos),
#        file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDvCTvSZ_cellmix.txt", quote = F,
#        row.names = F, sep = "\t")


# regression of PCs derived from predicted leukocyte abundances
groups = as.data.frame(t(combn(unique(as.character(mdsphenos$dx)), 2)))
groups$V1 = as.character(groups$V1)
groups$V2 = as.character(groups$V2)

coefs.out = list()
for( i in 1:nrow(groups)){
  tmp = mdsphenos[mdsphenos$dx %in% c(groups[i,1], groups[i,2]), ]
  y_var = tmp[,colnames(tmp) %in% paste("PC",1:ncol(mds), sep = ""), ]
  design = model.matrix( ~ dx , tmp)
  lmod = lm(as.matrix(y_var) ~ design)
  lmod = summary(lmod)
  lmod = lapply(lmod, function(x) broom::tidy(x$coefficients))
  coefs = ldply(lmod)
  coefs$comparison = paste(groups[i,], collapse = " - ")
  coefs.out[[i]] = coefs
}
coefs.out = ldply(coefs.out)
coefs.out = coefs.out[!grepl("Intercept", coefs.out$.rownames), ]
coefs.out$.id = gsub("Response ", "", coefs.out$.id)

coefs.out$.id = factor(coefs.out$.id, levels=paste("PC",1:nrow(coefs.out)/2, sep = ""))
coefs.out$.rownames = gsub("designdx", "", coefs.out$.rownames)

coefs.out$label = NA
coefs.out$label = ifelse(coefs.out$Pr...t.. < .05, '*', NA)

coefs.out[coefs.out$Pr...t.. < .05, ]

fwrite(setDT(coefs.out),
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDvSZvCT_cellMixSumStats.txt",
       quote = F,
       row.names = F, 
       sep = "\t")

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDvSZvCT_blood-PCA-cellmix.png", res = 300, units = "in", height = 6, width = 9)

ggplot(coefs.out, aes(x = .id, y = -log10(Pr...t..))) + 
  geom_bar(stat='identity', col = 'black') +
  # geom_text(data = coefs.out, aes(x = .id, y = -log10(Pr...t..), label = coefs.out$label)) +
  theme_bw() +
  facet_grid(~gsub("[-]", " v. ", toupper(comparison))) +
  theme(axis.text.x = element_text(hjust = 1, angle = 90)) +
  geom_hline(aes(yintercept = -log10(.05), col ='darkred'), linetype='dashed', lwd = 0.5) +
  xlab("PCs from predicted blood leukocyte abundances") + 
  ylab(expression(paste("-log"[10],"(P-value)"))) +
  guides(col = F)

dev.off()



# 3. set parameters for network algorithm
sft.power = 7; ## SFT 8, R-squared = 0.88013307 (SE = 0.05562321)
deepSplit = 2;
minModuleSize = 30;

# 4.one-step automated network construction
netA = blockwiseModules(datExpr0[,c(6:ncol(datExpr0))], 
                       power = sft.power, 
                       networkType = "signed",
                       deepSplit= deepSplit, 
                       TOMType = "signed", 
                       minModuleSize = minModuleSize, 
                       minCoreKME = 0.5, 
                       minCoreKMESize = minModuleSize/3,
                       minKMEtoStay=0, 
                       reassignThreshold = 1e-6, 
                       mergeCutHeight = 0.25, 
                       detectCutHeight = 0.995,
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE, 
                       pamStage = TRUE,
                       saveTOMs = TRUE, 
                       verbose = 3,
                       maxBlockSize = 5000)

netB = blockwiseModules(datExpr1[,c(6:ncol(datExpr1))], 
                        power = sft.power, 
                        networkType = "signed",
                        deepSplit= deepSplit, 
                        TOMType = "signed", 
                        minModuleSize = minModuleSize, 
                        minCoreKME = 0.5, 
                        minCoreKMESize = minModuleSize/3,
                        minKMEtoStay=0, 
                        reassignThreshold = 1e-6, 
                        mergeCutHeight = 0.25, 
                        detectCutHeight = 0.995,
                        numericLabels = TRUE, 
                        pamRespectsDendro = FALSE, 
                        pamStage = TRUE,
                        saveTOMs = TRUE, 
                        verbose = 3,
                        maxBlockSize = 5000)













# 5. save the network to a .Rdata file for future use
# saveRDS(net, file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDSZCT_WGCNA_signed.pwr8.Rdata")
# saveRDS(net, file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDSZCT_WGCNA_signed.pwr8_sep28-2017.Rdata")


# 6. extract network meta-data and eigengenes
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs

metadata = data.frame(colors = moduleColors, labels = paste("ME",net$colors, sep=""))
metadata = metadata[!duplicated(metadata$colors), ]

mergedColors = labels2colors(net$colors)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDSZCT_dendrogram.png", res=300, units="in",height = 5.5, width = 8)

plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dev.off()



module = list()
colors = unique(moduleColors)
for(y in 1:length(colors)){
  genesInModule = colnames(datExpr)[which(moduleColors %in% colors[[y]])]
  module[[y]] = data.frame(color = colors[[y]], label = unique(net$colors)[[y]], symbol = genesInModule)
}
module = ldply(module)

fwrite(setDT(module), 
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDSZCT_WGCNA_module-membership.txt",
       quote = F, row.names = F, sep = "\t")



# combine MEs and phenos 

datExpr = as.data.frame(rbindlist(list(szdata, mydata))) # combine into one matrix 
datExpr = datExpr[!is.na(datExpr$studyID), ]

datExpr = as.data.frame(datExpr)

datExpr$dx[datExpr$dx %in% c("CT", "control")] = "CT"

phenos = datExpr[rownames(datExpr) %in% rownames(MEs), c(1:6)]

datAll = data.frame(phenos, MEs)


datAll = merge(mdsphenos[,-c(2)],
                    datAll, by = 'sampleIDs')

fwrite(data.table(datAll),
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDvSZvCT_WGCNA_networkModuleEigengenesPhenotypes.txt",
       row.names = F, 
       quote = F, 
       sep = "\t")



# Statistical analysis 
# 1. make a model matrix

datAll = data.table(datAll)

group_set = unique(coefs.out$comparison)

dx.est.out = list()
for( i in 1:length(group_set)){
  
# extract cell type PCs that are significantly different between groups, add as covariates
cell_select = coefs.out[coefs.out$comparison %in% group_set[[i]], ]
cell_select = unique(as.character(cell_select$.id[cell_select$Pr...t.. < .05]))
  
deconv = strsplit(group_set[i], " - ")[[1]]
# set up model matrix
tmp.matrix = datAll[,colnames(datAll) %in% c("dx", "age",   "sex", cell_select),with=F]
tmp.matrix = tmp.matrix[tmp.matrix$dx %in% deconv, ]
tmp.matrix$sex[tmp.matrix$sex %in% "male"] = "M"
tmp.matrix$sex[tmp.matrix$sex %in% "female"] = "F"
table(tmp.matrix$dx)
design = model.matrix( ~ ., tmp.matrix )

# model equation
equation = paste(colnames(design)[-1], sep = "", collapse = " + ")

# 2. subset rows that are present in model matrix (NAs will be removed)
y = MEs
y  = y[rownames(y) %in% rownames(design), ]
design = design[rownames(design) %in% rownames(y), ]
# 3. fit linear model
lmod = lm(as.matrix(y) ~ -1 + design)
# 4. obtain summary for linear model (coefficients, standard errors, t-stats, p-values, r-squared of model, etc.)
sum = summary(lmod)
# 5. extract coefficients table from each model fit
coefs = lapply(sum , function(x) broom::tidy(x$coefficients))
coefs = ldply(coefs)
dx.est = coefs
# fix up table
dx.est$.id = gsub( "Response ", "", dx.est$.id)
names(dx.est)[names(dx.est) %in% '.id'] = 'labels'
dx.est = merge(metadata, dx.est, by = "labels")
dx.est = dx.est[!dx.est$labels %in% "ME0",] # discard ME0 --> "grey" module
dx.est = dx.est[order(dx.est$Pr...t.., decreasing = F), ]
dx.est$comparison = group_set[[i]]
dx.est$model = equation
dx.est.out[[i]] = dx.est

}
dx.df = ldply(dx.est.out)
dx.df = dx.df[grepl("designdx", dx.df$.rownames), ]
table(dx.df$comparison)
dx.df$labels = factor(dx.df$labels, levels = paste("ME", 1:ncol(MEs), sep = ""))


### SCZ + BD --> CT analysis 

  # extract cell type PCs that are significantly different between groups, add as covariates
  cell_select = unique(as.character(coefs.out$.id[coefs.out$Pr...t.. < .05]))
  
  # set up model matrix
  tmp.matrix = datAll[,colnames(datAll) %in% c("dx", "age", "sex", cell_select),with=F]
  tmp.matrix$dx[tmp.matrix$dx %in% c("scz", "BD")] = "psychosis"
  tmp.matrix$sex[tmp.matrix$sex %in% "male"] = "M"
  tmp.matrix$sex[tmp.matrix$sex %in% "female"] = "F"
  table(tmp.matrix$dx)
  design = model.matrix( ~ .,tmp.matrix )
  
  # model equation
  equation = paste(colnames(design)[-1], sep = "", collapse = " + ")
  
  # 2. subset rows that are present in model matrix (NAs will be removed)
  y = MEs
  y  = y[rownames(y) %in% rownames(design), ]
  design = design[rownames(design) %in% rownames(y), ]
  # 3. fit linear model
  lmod = lm(as.matrix(y) ~ design)
  # 4. obtain summary for linear model (coefficients, standard errors, t-stats, p-values, r-squared of model, etc.)
  sum = summary(lmod)
  # 5. extract coefficients table from each model fit
  coefs = lapply(sum , function(x) broom::tidy(x$coefficients))
  coefs = ldply(coefs)
  psy.est = coefs
  psy.est = psy.est[grepl("dx", psy.est$.rownames), ]
  # fix up table
  psy.est$.id = gsub( "Response ", "", psy.est$.id)
  names(psy.est)[names(psy.est) %in% '.id'] = 'labels'
  psy.est = merge(metadata, psy.est, by = "labels")
  psy.est = psy.est[!psy.est$labels %in% "ME0",] # discard ME0 --> "grey" module
  psy.est = psy.est[order(psy.est$Pr...t.., decreasing = F), ]
  psy.est$comparison = paste("psychosis - CT")
  psy.est$model = equation
  
dx.all = ldply(list(dx.df, psy.est))
dx.all$bonf = p.adjust(dx.all$Pr...t.., "bonferroni")
dx.all$fdr = p.adjust(dx.all$Pr...t.., "fdr")

dx.all$comparison[dx.all$comparison %in% "CT - BD"] = "Control v. BD"
dx.all$comparison[dx.all$comparison %in% "psychosis - CT"] = "Control v. SCZ + BD"
dx.all$comparison[dx.all$comparison %in% "scz - BD"] = "SCZ v. BD"
dx.all$comparison[dx.all$comparison %in% "scz - CT"] = "Control v. SCZ"

MEsizes = data.frame(table(module$label))
colnames(MEsizes) = c("labels",'MEfreq')
MEsizes$labels = paste("ME", MEsizes$labels, sep = "")
MEsizes$colFreq = paste(MEsizes$labels, " (", MEsizes$MEfreq, ")", sep = "")

dx.all = merge(MEsizes, dx.all, by= "labels")

dx.all = dx.all[order(dx.all$Pr...t..,decreasing=F),]

head(dx.all)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDvSZvCT_WGCNA_DiffExpr_barplot.png", res = 300, units= "in", height = 8, width = 11.5)

ggplot(dx.all, aes(x = colFreq,  y = -log10(Pr...t..))) + 
  geom_bar(stat='identity', col = 'black', width=0.75, position='dodge', aes(fill = comparison)) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1, angle = 90)) + 
  xlab(NULL) + 
  guides(col = F) +
  ylab(expression(paste("-log"[10],"(P-value)")))+ 
  geom_hline(aes(yintercept = -log10(.05)), col = 'red', linetype=  'dashed')  +
  geom_hline(aes(yintercept = -log10(.05/(nrow(dx.all)))), col = 'green', linetype=  'solid')  +
  labs(fill = 'Group comparisons', col = c("P < 0.05", "Bonferroni P < .05"))  +
  scale_fill_brewer(palette="Paired") 
  
dev.off()



## Perform gene set enrichment analysis ...
msig = list.files(path = "~/Google Drive/mac_storage/gene_sets/", pattern = "msig.txt", full.names = T)
msig = lapply(msig, function(x) fread(x, h = T))
msig = rbindlist(msig)

sig_mods = as.character(dx.all$colors[dx.all$bonf < .05])
dx.all[dx.all$colors %in% sig_mods, ]

hypStats_save = list()
for( i in 1:length(sig_mods)){
  cat("\rRunning hypergeometric test for module:",sig_mods[[i]])
  gene_grab = module[module$color %in% sig_mods[[i]], ]
  
  msig = msig[msig$SYMBOL %in% module$symbol]
  msig = msig[grepl("chr|KEGG_|GO_|REACTOME|PANTHER", msig$PATHWAY)]
  msig = msig[!is.na(msig$SYMBOL)]
  
  
  # filter sets
  set_count = table(as.character(msig$PATHWAY))
  set_count = set_count[set_count >= 10 & set_count <= 200]
  
  list = unique(names(set_count))
  
  msig.keep = msig[msig$PATHWAY %in% list]
  
  msig.split = split(msig, msig$PATHWAY)
  
  # calculate hypergeometric input values
  universe = length(unique(module$symbol)) # total # of genes in transcriptome 
  overlaps = lapply(msig.split, function(x) length(intersect(x$SYMBOL, gene_grab$symbol)))
  geneSetSize = lapply(msig.split, nrow)
  geneSetOut = universe - unlist(geneSetSize)
  
  # hypergeomtric p-values
  hypStats = data.frame(
    Module = sig_mods[[i]],
    Set_name = names(msig.split),
    Population = universe, 
    Sample_success = unlist(overlaps),
    Population_success = unlist(geneSetSize),
    Population_fails = unlist(geneSetOut),
    Sample_size = length(unique(gene_grab$symbol))
  )
  
  # enrichment p-value test
  pvals = phyper(hypStats$Sample_success - 1, hypStats$Population_success, hypStats$Population_fails, hypStats$Sample_size, lower.tail = FALSE, log.p = FALSE)
  
  hypStats$P = pvals
  
  hypStats = hypStats[order(hypStats$P, decreasing = F), ]
  hypStats$FDR = p.adjust(hypStats$P, "fdr")
  hypStats$BONF = p.adjust(hypStats$P, "bonferroni")
  
  hypStats_save[[i]] = hypStats
  
}
hypStats_save = ldply(hypStats_save)

hypStats_save[hypStats_save$BONF < .05, ]

sets = strsplit(as.character(hypStats_save$Set_name), "_")
pfx = lapply(sets, function(x) x[[1]])
suffix = lapply(sets, function(x) x[-1])

hypStats_save$relabel = NA
for( i in 1:length(sets)){
  sfx = paste(unlist(suffix[[i]]), collapse = " ", sep = "")
  hypStats_save$relabel[[i]] = paste(pfx[[i]], ":", sfx, sep = "")
}

hypStats_save$label = paste(hypStats_save$relabel, "\n(", hypStats_save$Sample_success,"/", hypStats_save$Population_success, ")", sep = "")
hypStats_save$p_label = paste("P = ", format(hypStats_save$P,scientific=T,digits=3), sep = "")

threshold = -log10(.05 / nrow(hypStats_save))

split = split(hypStats_save, hypStats_save$Module)
split = lapply(split, function(x) x[1:50, ])
hypStats_bind = ldply(split)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDvSZvCT_WGCNA_hypergeometricGeneset.png", res = 300,units="in",height=7,width=14)

ggplot(hypStats_bind, aes(x = relabel, y = -log10(P))) + 
  xlab(NULL) +
  facet_wrap(~Module) +
  ylab("-log10(Hypergeometric test P-value)") +
  geom_bar(stat = "identity", col = "black", fill = hypStats_bind$Module) +
  geom_hline(yintercept = threshold, col = "red", linetype = "solid", lwd = 0.3) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 4)) +
  geom_text(aes(label = p_label, y = threshold), nudge_y = 3, size=1.5) + 
  coord_flip() 

dev.off()
