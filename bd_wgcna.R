
# load packages
require(data.table)
require(WGCNA)
require(plyr)
require(ggplot2)
require(igraph)
require(massiR)
require(biomaRt)
require(org.Hs.eg.db)
require(CellMix)

# read data 
# data = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_7studiescombined_zscaled.txt", header = T)
data = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/expr_covs/BD_allstudies_standardized.txt",h=T)
# data = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_combat_7studies.txt",h=T)

data= data[!duplicated(data$FACTOR_sampleIDs), ]
data$FACTOR_sex = ifelse(grepl("witt", data$FACTOR_studyID), "M", data$FACTOR_sex)
data$FACTOR_sex[data$FACTOR_sex %in% ""] = NA


round(table(data$FACTOR_studyID, data$FACTOR_sex)/rowSums(table(data$FACTOR_studyID, data$FACTOR_sex))*100, 1)
table(data$FACTOR_studyID, data$FACTOR_dx)

aggregate(data$FACTOR_age, by=list(data$FACTOR_studyID), function(x) round(mean(x),2))
aggregate(data$FACTOR_age, by=list(data$FACTOR_studyID), function(x) sd(x)/sqrt(length(x)))


res_df = fread("~/Google Drive/mac_storage/TWAS/bd_mega/freeze_qcBD_blood_meta_Nmin4_LeaveOneOut.txt",h=T)
trim_df = res_df[res_df$Nstudy >= 4 & !is.na(res_df$ENTREZID)]

mean.age = aggregate(data$FACTOR_age, by=list(data$FACTOR_studyID), function(x) mean(x, na.rm = T))
mean.age$x = round(mean.age$x, 2)
sd.age = aggregate(data$FACTOR_age, by=list(data$FACTOR_studyID), function(x) sd(x, na.rm = T))
sd.age$x = round(sd.age$x/sqrt(table(data$FACTOR_studyID)) + 1.96, 2)
paste(mean.age$x, " + ", sd.age$x, sep = "")


datExpr = data[!is.na(data$FACTOR_sex) | !is.na(data$FACTOR_age), ]
datExpr = datExpr[,!grepl("FACTOR_",colnames(datExpr)),with=F]


datExpr = datExpr[,colnames(datExpr) %in% trim_df$GeneSymbol,with=F]

datExpr = as.data.frame(datExpr)



## 

gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

keepSamples = data[,grepl("FACTOR_",colnames(data)),with=F]
keepSamples = as.data.frame(keepSamples)
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

rownames(keepSamples) = keepSamples$FACTOR_sampleIDs
rownames(datExpr) = keepSamples$FACTOR_sampleIDs

disableWGCNAThreads() # number of processors reported by system will be used 

table(data$FACTOR_sex, data$FACTOR_dx)

# surrogate variable analysis 
## Surrogate variable analysis - default method 
mod = model.matrix(~ FACTOR_dx + FACTOR_sex + FACTOR_age + FACTOR_studyBatch, data=keepSamples) # model with known factors and covariates
mod0 = model.matrix(~1,data=keepSamples[rownames(keepSamples) %in% rownames(mod), ]) # intercept only model
nonMissingExpr = datExpr[,!colSums(is.na(datExpr))>1]
nonMissingExpr = nonMissingExpr[rownames(nonMissingExpr) %in% rownames(mod), ]
exprs = ExpressionSet(as.matrix(t(nonMissingExpr)))

require(sva)
n.sv = num.sv(exprs(exprs), mod, method="leek", B = 20) # number of significant surrogate variables.. Using "be" works well with smaller samples

if(n.sv > 0){
  svobj = sva(exprs(exprs),mod,mod0,n.sv=n.sv)
  svdf = data.frame(NULL)
  svdf = as.data.frame(svobj$sv)
  colnames(svdf) = paste("SV",1:ncol(svdf), sep = "")
  phens = data.frame(mod, svdf)
}

t.test(phens$SV1 ~ phens$FACTOR_dxCT)

exprs = as.data.frame(datExpr) # extract normalized gene expression intensities for subjects
exprs = exprs[,colSums(is.na(exprs))==0]

exprs = exprs

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

res = gedBlood(exprs, verbose = TRUE, normalize = TRUE) # run gedBlood algorithm (non-negative matrix factorization)

wb.coef = coef(res)
wb.coef = as.data.frame(t(wb.coef))

nrow(wb.coef) == nrow(phens)

wb_res = ldply(lapply(wb.coef, function(x) broom::tidy(t.test(x ~ phens$FACTOR_dxCT))))
wb_res = wb_res[order(wb_res$p.value, decreasing = F), ]
head(wb_res)


wb_res  = lm(as.matrix(wb.coef) ~ keepSamples$FACTOR_studyID + keepSamples$FACTOR_age + keepSamples$FACTOR_sex + keepSamples$FACTOR_dx)
wb_res = summary(wb_res)
wb_res = lapply(wb_res, function(x) broom::tidy(x$coefficients))
wb_res = ldply(wb_res)
wb_res = wb_res[grepl("FACTOR_dx", wb_res$.rownames), ]
wb_res = wb_res[order(wb_res$Pr...t.., decreasing = F), ]
wb_res$FDR = p.adjust(wb_res$Pr...t.., "fdr")
wb_res[wb_res$Pr...t.. < .05, ]

# correlation of blood cell concentrations to surrogate variables
cor_res = lapply(wb.coef, function(x) broom::tidy(cor.test(x, phens$SV1)))
cor_res = ldply(cor_res)
cor_res = cor_res[order(cor_res$p.value, decreasing = F), ]
dim(cor_res[cor_res$p.value < .05])
cor_res[cor_res$p.value < .05, ]



## ===== Module preservation of control gene network in cases =====  # 
caseExpr = datExpr[which(data$FACTOR_dx %in% "BD"), ]
controlExpr = datExpr[which(data$FACTOR_dx %in% "CT"), ]

multiExpr = list(caseExpr, controlExpr)
lapply(multiExpr, dim)

for(x in 1:length(multiExpr)){
  
  temp = multiExpr[[x]]
  gsg = goodSamplesGenes(temp, verbose = 3);
  gsg$allOK
  
  keepSamples = data[rownames(data) %in% rownames(temp), ]
  
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(temp)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(temp)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    temp = temp[gsg$goodSamples, gsg$goodGenes]
  }
  
  gsg = goodSamplesGenes(temp, verbose = 3);
  gsg$allOK
  
  multiExpr[[x]] = temp
  
}

intersect_genes = intersect(colnames(multiExpr[[1]]), colnames(multiExpr[[2]]))

multiExpr = lapply(multiExpr, function(x) x[,colnames(x) %in% intersect_genes])

# Create gene co-expression network in controls

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft1 = pickSoftThreshold(multiExpr[[1]],
                         powerVector = powers,
                         corFnc="bicor",
                         networkType="signed")

sft.power = 12;
deepSplit = 2;
minModuleSize = 30;

control_net = blockwiseModules(multiExpr[[1]], 
                       power = sft.power, 
                       networkType = "signed",
                       deepSplit= deepSplit, 
                       TOMType = "signed", impute = FALSE,
                       minModuleSize = minModuleSize, 
                       minCoreKME = 0.5, 
                       minCoreKMESize = minModuleSize/3,
                       minKMEtoStay=0, 
                       reassignThreshold = 1e-6, 
                       mergeCutHeight = 0.25, 
                       detectCutHeight = 0.995,
                       numericLabels = TRUE, 
                       corType = 'bicor',
                       pamRespectsDendro = FALSE, 
                       pamStage = TRUE,
                       saveTOMs = TRUE, 
                       verbose = 3,
                       maxBlockSize = 5000)

case_net = blockwiseModules(multiExpr[[2]], 
                               power = sft.power, 
                               networkType = "signed",
                               deepSplit= deepSplit, 
                               TOMType = "signed", impute = FALSE,
                               minModuleSize = minModuleSize, 
                               minCoreKME = 0.5, 
                               minCoreKMESize = minModuleSize/3,
                               minKMEtoStay=0, 
                               reassignThreshold = 1e-6, 
                               mergeCutHeight = 0.25, 
                               detectCutHeight = 0.995,
                               numericLabels = TRUE, 
                               corType = 'bicor',
                               pamRespectsDendro = FALSE, 
                               pamStage = TRUE,
                               saveTOMs = TRUE, 
                               verbose = 3,
                               maxBlockSize = 5000)

setLabels = c("Control", "Case");
nSets = 2;
# Object that will contain the expression data
controlExpr = multiExpr[[1]]
caseExpr = multiExpr[[2]]
multiExpr = list();
multiExpr[[1]] = list(data = controlExpr);
multiExpr[[2]] = list(data = caseExpr);
names(multiExpr) = setLabels
multiColor = list(Control = labels2colors(control_net$colors), Case = labels2colors(case_net$colors));

mp = modulePreservation(multiData = multiExpr, 
                        loadPermutedStatistics = FALSE,
                        dataIsExpr = TRUE,
                        multiColor = multiColor, 
                        verbose = 3,
                        nPermutations = 100,
                        networkType = 'signed', 
                        corFnc = 'bicor',
                        referenceNetworks = 1)


ref = 1 # Select the human data as reference
test = 2 # Select the chimp data as test
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

print(signif(statsZ[, "Zsummary.pres", drop = FALSE],2));
# Compare preservation to quality:
print(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

ref = 1;
test = 2;
ind = 1;
stats= mp$preservation$observed[[ref]][[test]];
labelsX = rownames(stats)
labelsX[labelsX=="gold"] = "orange"
modColors = labelsX;
plotMods = !(modColors %in% c("grey", "orange"));
moduleSizes = stats[plotMods, 1];
textLabels = match(modColors, standardColors(20))[plotMods];
colorLabels = labelsX[plotMods];

stats = data.frame(Module = rownames(stats), stats, Zsummary = statsZ)
stats
write.csv(stats, file="~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/WGCNA_preservation_BD-CT.csv",quote=F,row.names=F)

Zscore = mp$preservation$Z[[ref]][[test]]$Zsummary.pres[plotMods];
Rank = mp$preservation$observed[[ref]][[test]]$medianRank.pres[plotMods];
plot_df = data.frame(Zscore, Rank, ModuleSize = moduleSizes)
plot_df$colors = colorLabels

ggplot(plot_df, aes(x = ModuleSize, y = Zscore, col = colors)) +
  geom_point(size = 4, pch = 21, bg = as.character(plot_df$colors), col = 'black') +
  theme_classic() +
  xlab("Module size") +
  ylab("Module preservation (z-score)") +
  geom_hline(aes(yintercept = 2), col = 'blue', lwd = 0.5, lty = 2 ) +
  geom_hline(aes(yintercept = 0), col = 'black', lwd=0.35, lty = 2) +
  geom_hline(aes(yintercept = 10), col = 'salmon', lwd=0.75, lty = 2) +
  theme(panel.border = element_rect(size = 1, fill = NA), axis.text=element_text(size = 12, colour = 'black'), axis.title=element_text(size = 14))


res_df = fread("~/Google Drive/mac_storage/TWAS/bd_mega/freeze_qcBD_blood_meta_Nmin4_LeaveOneOut.txt",h=T)
res_df = data.frame(res_df)

# --> compare DEG list with module preservation 

colors = control_net$colors
colors = labels2colors(colors)

moduleLabels = data.frame(Colors = colors, GeneSymbol = colnames(multiExpr$Control$data))

res_df = merge(moduleLabels, res_df, by='GeneSymbol')
res_df = res_df[!is.na(res_df$P), ]

mean_p_val = ddply(res_df, ~Colors, summarize, Mean = mean(-log10(P)), SE = sd(-log10(P))/sqrt(length(P)), N = length(P))
mean_p_val = mean_p_val[order(mean_p_val$Mean, decreasing = T), ]

mean_p_val = mean_p_val[match(plot_df$colors, mean_p_val$Colors), ]

plot_df$meanGeneSig = mean_p_val$Mean

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/WGCNA_preservation_zscore_geneSigSize.png",res=300,units="in",height=5,width=6)
ggplot(plot_df, aes(x = ModuleSize, y = Zscore, col = colors)) +
  geom_point(size = (1+plot_df$meanGeneSig)^2, pch = 21, bg = as.character(plot_df$colors), col = 'black') +
  theme_classic() +
  xlab("Module size") +
  ylab("Module preservation (z-score)") +
  geom_hline(aes(yintercept = 2), col = 'blue', lwd = 0.5, lty = 2 ) +
  geom_hline(aes(yintercept = 0), col = 'black', lwd=0.35, lty = 2) +
  geom_hline(aes(yintercept = 10), col = 'salmon', lwd=0.75, lty = 2) +
  theme(panel.border = element_rect(size = 1, fill = NA), axis.text=element_text(size = 12, colour = 'black'), axis.title=element_text(size = 14))
dev.off()

## one-step automated gene network analysis 
# 1. find optimal soft-threshold power for network construction
powers = c(c(1:10), seq(from = 12, to=20, by=2))

iters=5
sft1 = list()
for( y in 1:iters){
  random = datExpr[,sample(colnames(datExpr), 2000, replace = FALSE)]
  sft1[[y]] = pickSoftThreshold(random,
                                powerVector = powers,
                                corFnc="bicor",
                                networkType="signed")
}
sft1_power = lapply(sft1, function(x) x[[1]])
sort(table(unlist(sft1_power)),decreasing = T)
sft1_in = lapply(sft1, function(x) x$fitIndices)
names(sft1_in) = iters
sft1_df = ldply(sft1_in)
sft1_avg = aggregate(sft1_df$SFT.R.sq, by = list(sft1_df$Power), mean)
sft1_se = aggregate(sft1_df$SFT.R.sq, by = list(sft1_df$Power), function(x) sd(x)/sqrt(length(x)))
sft1_avg$CI_LOWER = sft1_avg$x - ( 1.96 * sft1_se$x)
sft1_avg$CI_HIGH = sft1_avg$x + ( 1.96 * sft1_se$x)
names(sft1_avg)[names(sft1_avg) %in% "x"] = "mean"

sft1 = pickSoftThreshold(datExpr,
                  powerVector = powers,
                  corFnc="bicor",
                  networkType="signed")

r2 = ggplot(sft1_avg, aes(x = factor(Group.1), y = mean)) +
  geom_bar(stat = 'identity', fill = heat.colors(15), colour = 'black') + 
  theme_classic() +
  xlab("Soft threshold power") +
  ylab(expression(paste("Model fit ",italic(R)^2))) +
  geom_errorbar(aes(ymin = CI_LOWER, ymax = CI_HIGH), width=0.2) +
  geom_hline(aes(yintercept = .8),col='navyblue',linetype='dashed') + 
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))

print(r2)

# 2. figure of possible SFT thresholds to pick from
png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/SFT_WGCNA_BD-signednetwork.png",res=300,units="in",height = 5, width=5)
# plot(x = sft1$fitIndices$Power, 
#      y = sft1$fitIndices$SFT.R.sq, 
#      type = "o", 
#      pch = 19, 
#      col = "red",
#      xlab = "Soft threshold power", 
#      ylab ="Model fit (R-squared)")
# abline(h = 0.8, col = "darkblue", lty = 2)
print(r2)
dev.off()

## Combat procedure to remove batch effects (study ID)

modcombat = model.matrix(~1, data=phenos)
edata = datExpr[,colSums(is.na(datExpr))==0]
combat_edata = ComBat(dat= t(edata), batch=phenos$FACTOR_studyBatch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

edata = data.frame(t(combat_edata))


pwr = pickSoftThreshold(edata,
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
sft.power = 12; ## SFT 8, Rsq = 0.80982071 (SE = 0.05211014)
deepSplit = 2;
minModuleSize = 30;


# 4.one-step automated network construction
net = blockwiseModules(datExpr, 
                       power = sft.power, 
                       networkType = "signed",
                       deepSplit= deepSplit, 
                       TOMType = "signed", 
                       impute = FALSE,
                       minModuleSize = minModuleSize, 
                       minCoreKME = 0.5, 
                       minCoreKMESize = minModuleSize/3,
                       minKMEtoStay=0, 
                       reassignThreshold = 1e-6, 
                       mergeCutHeight = 0.25, 
                       detectCutHeight = 0.995,
                       numericLabels = TRUE, 
                       corType = 'bicor',
                       pamRespectsDendro = FALSE, 
                       pamStage = TRUE,
                       saveTOMs = TRUE, 
                       verbose = 3,
                       maxBlockSize = 5000)

# 5. save the network to a .Rdata file for future use
# saveRDS(net, file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_7studiescombined_WGCNA_signed.pwr9.Rdata")
# saveRDS(net, file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_7studiescombined_WGCNA_signed.pwr12.Rdata")
# saveRDS(net, file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_7studiescombined_WGCNA_signed.pwr14.Rdata")
# saveRDS(net, file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/freeze_BD_7studiescombined_WGCNA_signed.pwr12.Rdata")
# saveRDS(net, file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/freeze_BD_7studiescombined_WGCNA_signed.pwr12_noimpute.Rdata")

net = readRDS("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/freeze_BD_7studiescombined_WGCNA_signed.pwr12.Rdata")


# 6. extract network meta-data and eigengenes
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
rownames(MEs) = rownames(datExpr)

message("\rDetected ", ncol(MEs), " modules!")

metadata = data.frame(colors = moduleColors, labels = paste("ME",net$colors, sep=""))
metadata = metadata[!duplicated(metadata$colors), ]

mergedColors = labels2colors(net$colors)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_dendrogram.png", res=300, units="in",height = 5.5, width = 8)

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
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_WGCNA_module-membership.txt", 
       quote = F, row.names = F, sep = "\t")


# combine MEs and phenos 
phenos = keepSamples[rownames(keepSamples) %in% rownames(datExpr), grepl("FACTOR_",colnames(keepSamples))]
order = rownames(MEs)
datAll = data.frame(phenos, MEs[rownames(MEs) %in% rownames(phenos), ])
a = datAll
phens$FACTOR_sampleIDs = phenos$FACTOR_sampleIDs
datAll = merge(phens[,c("SV1","FACTOR_sampleIDs")], datAll, "FACTOR_sampleIDs", all.y=T)
datAll = datAll[match(a$FACTOR_sampleIDs, datAll$FACTOR_sampleIDs), ]

# 
# # 
# # # === CellMix section
# exprs = as.data.frame(datExpr) # extract normalized gene expression intensities for subjects
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
# # wb.coef = wb.coef[,!colSums(wb.coef == 0)/nrow(wb.coef) > 0.99]
# wb.pca = prcomp(wb.coef)
# pcaCov = as.data.frame(wb.pca$x)
# 
# lmFit = lm(as.matrix(pcaCov) ~ phenos$FACTOR_dx + phenos$FACTOR_sex + phenos$FACTOR_age)
# lmCoef = summary(lmFit)
# lmCoef = lapply(lmCoef, function(x) broom::tidy(x$coefficients))
# names(lmCoef) = colnames(pcaCov)
# lmdf = ldply(lmCoef)
# lmdf = lmdf[grepl("dx",lmdf$.rownames), ]
# lmdf = lmdf[order(lmdf$Pr...t..,decreasing=F),]
# lmdf$.id = factor(lmdf$.id , levels=paste("PC",1:nrow(lmdf),sep=""))
# # 
# coefs = lmdf
# 
# pickPC = lmdf$.id[lmdf$Pr...t.. < 5e-02]
# pickPC = pickPC[!is.na(pickPC)]
# 
# ggplot(coefs, aes( x = .id, y = -log10(Pr...t..))) + 
#   geom_bar(stat = 'identity') + 
#   geom_hline(aes(yintercept = -log10(.05)), col = 'red', linetype = 'dashed') +
#   theme_bw() +
#   theme(axis.text.x = element_text(hjust = 1, angle = 90)) + 
#   xlab(NULL) + 
#   ylab(expression(paste("-log"[10],"(P)")))
# 
# # datAll = data.frame(leuk.pc, datAll)
# 
# pcaCov = pcaCov[match(rownames(datAll), rownames(pcaCov)), ]
# 
# # datAll = data.frame(datAll, pcaCov)
# 
# datAll = datAll
datAll$FACTOR_medicated[datAll$FACTOR_medicated %in% ""] = NA
datAll$FACTOR_sex[datAll$FACTOR_sex %in% ""] = NA
datAll$FACTOR_race[datAll$FACTOR_sex %in% ""] = NA
datAll$FACTOR_dx = factor(datAll$FACTOR_dx, levels=c("CT","BD"))

datAll = datAll[!is.na(datAll$FACTOR_sex) & !is.na(datAll$FACTOR_age), ]

table(datAll$FACTOR_sex, datAll$FACTOR_dx)

# fwrite(datAll, file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_WGCNA_MEphenos.txt",
#        quote = F,
#        row.names = F,
#        sep = "\t")


# Statistical analysis 
y = datAll[,grepl("ME", colnames(datAll))]

# lmer - random effect regression
parms = list()
for(n in 1:ncol(y)){
  frm = formula(paste(colnames(y)[[n]], " ~ FACTOR_age + FACTOR_sex + FACTOR_dx + SV1 + (1|FACTOR_studyBatch)"))
  fit = lmer(frm, data = datAll)  
  parms[[n]] = broom::tidy(summary(fit)$coefficients)
}
names(parms) = colnames(y)
parms = ldply(parms)
parms = parms[parms$.rownames %in% "FACTOR_dxBD", ]
parms = parms[order(parms$Pr...t..), ]
colnames(parms)[1:2] = c("labels", "predictors")
parms = merge(metadata, parms, by='labels')
parms$FDRp = p.adjust(parms$Pr...t.., 'fdr')
parms = parms[order(parms$Pr...t.., decreasing = F), ]

# linear model - estimate Betas with normal equation
lmod = lm(as.matrix(y) ~ FACTOR_studyBatch + FACTOR_sex + FACTOR_age + FACTOR_dx + SV1, data = datAll)
# lmod = lm(as.matrix(y)  ~ FACTOR_dx, data = datAll)
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
dx.est = coef_all[coef_all$predictors %in% "FACTOR_dxBD", ]
dx.est = merge(metadata, dx.est, by = "labels")
dx.est = dx.est[!dx.est$labels %in% "ME0",] # discard ME0 --> "grey" module
dx.est = dx.est[order(dx.est$`Pr(>|t|)`, decreasing = F), ]
dx.est$fdr = p.adjust(dx.est$`Pr(>|t|)`, "fdr")
dx.est$bonf = p.adjust(dx.est$`Pr(>|t|)`, "bonferroni")

head(dx.est)

dx.est[dx.est$fdr < .05, ]

fwrite(dx.est, 
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/BD_WGCNA_lmresults.txt",
       quote = F, row.names = F, sep = "\t")

dx.est$lab_point = ifelse(dx.est$bonf < .05, as.character(dx.est$labels), NA)

require(ggrepel)
png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/wgcna_bd_volcano.png",res=300,units="in",height=6,width=6)
ggplot(dx.est, aes(x = Estimate, y = logP)) +
  geom_point(size = 3, pch = 21, col = 'black', aes(bg = I(colors))) +
  geom_text_repel(aes(label = lab_point), size = 4) + 
  geom_hline(aes(yintercept = -log10(.05)), col = 'orange', linetype=2) + 
  geom_hline(aes(yintercept = -log10((.05)/nrow(dx.est))), col = 'red') +
  theme_classic() +
  ylab(expression(paste("-log"[10],"(P-value)"))) + 
  xlab(expression(paste("Effect size (",beta,")")))  +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15), strip.text = element_text(size = 15)) 
dev.off()

require(ggthemes)

dx.est$condition = 'BD vs. CT'
dx.est$LABELS = paste( round(dx.est$Estimate, digits = 3), "\n", "(",round(dx.est$`Std. Error`,3),")", sep = "")


DIAG = dx.est


# grab effects for diagnostic status and adjust p-values for multiple testing
# dx.est = coef_all[grepl(":", coef_all$predictors), ]
# dx.est = merge(metadata, dx.est, by = "labels")
# dx.est = dx.est[!dx.est$labels %in% "ME0",] # discard ME0 --> "grey" module
# dx.est = dx.est[order(dx.est$`Pr(>|t|)`, decreasing = F), ]
# dx.est$fdr = p.adjust(dx.est$`Pr(>|t|)`, "fdr")
# dx.est$bonf = p.adjust(dx.est$`Pr(>|t|)`, "bonferroni")

head(dx.est)

# png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/ME15_sex-by-diagnosis.png",res=300,units="in",height=5,width=5.5)
# ggplot(datAll[!is.na(datAll$FACTOR_sex), ], aes(y = ME15, x = FACTOR_dx, fill = FACTOR_sex)) +
#   geom_boxplot(lwd=0.7, outlier.colour = 'grey') + 
#   theme_classic() +
#   scale_fill_brewer(NULL, palette = 'Set1') +
#   xlab('Diagnosis') + 
#   ylab("ME15 eigengene expression") + 
#   theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12))
# dev.off()

# dx.est$condition = 'Sex-by-diagnosis'
# dx.est$LABELS = paste(round(dx.est$Estimate, 4), "\n", "(",round(dx.est$`Std. Error`,3),")", sep = "")

# tile plot
# graphDF = ldply(list(dx.est, DIAG))
graphDF = DIAG
graphDF$condition = factor(graphDF$condition, levels=unique(graphDF$condition))
graphDF$labels = factor(graphDF$labels, levels = paste("ME",1:100,sep=""))
graphDF$LABELS = ifelse(graphDF$fdr < .05, paste("*",graphDF$LABELS, sep = ''), graphDF$LABELS)
graphDF$FONT = 1
graphDF$FONT[grepl("[*]", graphDF$LABELS)] = 2

require(ggthemes)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BDCT_tile.png",res=300,units="in",height=4.5,width=10)
ggplot(graphDF, aes(x = labels, y = condition, fill = logP)) +
  geom_tile(color = 'black',  size = 0.2, show.legend = TRUE) +
  theme_tufte(base_family = "Helvetica") +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_continuous(expression(paste('-log'[10],'(P-value)')), low = 'green',  high = c("white","salmon")) +
  theme(axis.ticks = element_line(colour = NA), legend.title = element_text(size=8), axis.text.x = element_text(hjust = 1, angle = 90, size = 8), axis.text = element_text(size=10) ) + 
  theme(legend.position =  'bottom') +
  geom_text(aes(label = LABELS, fontface = FONT), size = 1.75) +
  coord_equal()
dev.off()


# === stratify by gender
gender = as.character(unique(datAll$FACTOR_sex))
gender = gender[!is.na(gender)]

gender.list = list()
for(x in 1:length(gender)){
  
# Statistical analysis 
# 1. make a model matrix
dat = datAll[datAll$FACTOR_sex %in% gender[[x]], ]
groupCounts = table(dat$FACTOR_dx)
y = dat[,grepl("ME",colnames(dat))]

lmod = lm(as.matrix(y) ~ FACTOR_age + FACTOR_studyID + FACTOR_dx + SV1, data=dat)
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
gender.est = coef_all[coef_all$predictors %in% "FACTOR_dxBD", ]
gender.est = merge(metadata, gender.est, by = "labels")
gender.est = gender.est[!gender.est$labels %in% "ME0",] # discard ME0 --> "grey" module
gender.est = gender.est[order(gender.est$`Pr(>|t|)`, decreasing = F), ]
gender.est$fdr = p.adjust(gender.est$`Pr(>|t|)`, "fdr")
gender.est$bonf = p.adjust(gender.est$`Pr(>|t|)`, "bonferroni")
gender.est = gender.est[,!colnames(gender.est) %in% ".id"]
gender.est$cases = groupCounts[[1]]
gender.est$controls = groupCounts[[2]]

gender.list[[x]] = gender.est
}

names(gender.list) = gender
gender.df = ldply(gender.list)
names(gender.df)[names(gender.df) %in% '.id'] = 'gender'
gender.df = gender.df[order(gender.df$`Pr(>|t|)`, decreasing = F), ]

gender.df$FDR = p.adjust(gender.df$`Pr(>|t|)`, "fdr")

head(gender.df)

gender.df$lab_point = ifelse(gender.df$fdr < .05, as.character(gender.df$labels), NA)

# = volcano plot : association of module eigengenes with diagnosis, stratified by gender
counts = nrow(gender.df)/2

gender.df$gender = ifelse(gender.df$gender %in% "F", 
                          "Female BD  vs. Female CT", 
                          "Male BD vs. Male CT")


png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/wgcna_BD-gender-stratified.png",res=300,units="in",height=5,width=8)

ggplot(gender.df, aes(x = Estimate, y = logP)) +
  geom_point(size = 3, pch = 21, col = 'black', aes(bg = I(colors))) +
  geom_text_repel(aes(label = lab_point), size = 4) + 
  geom_hline(aes(yintercept = -log10(.05)), col = 'orange', linetype=2) + 
  geom_hline(aes(yintercept = -log10((.05)/counts)), col = 'red') +
  theme_classic() +
  facet_wrap( ~ gender) + 
  ylab(expression(paste("-log"[10],"(P-value)"))) + 
  xlab(expression(paste("Effect size (",beta,")")))  +
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15), strip.text = element_text(size = 15)) +
  scale_x_continuous(breaks = pretty(gender.df$Estimate, n = 3)) 

dev.off()
  


# plot p-values for predictors
png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_WGCNA_ggplot_sumstats.png", res=300,units="in",height = 7, width = 8)
ggplot(coef_all, aes( y = coef_all$logP,x = coef_all$predictors, fill = factor(coef_all$predictors))) + 
  geom_violin() +
  guides(fill = F) + 
  xlab(NULL) +
  theme_bw() +
  ylab(expression(paste("-log"[10],"(P-values)"))) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45)) + 
  geom_boxplot(width = 0.05, outlier.shape = NA, fill = "lightgrey") + 
  geom_hline(yintercept = -log10(.05), col = "red", linetype = 'dashed')
dev.off()

# plot -log10 P-values for each module
barplot(dx.est$Estimate, ylim=c(-0.1,0.1), 
        col = as.character(dx.est$colors), 
        ylab= "Differential expression effect size (+/- SE)")
mtext(side = 1, line = 0.1, "Modules ranked by significance")
# vertical bars
segments(x0 = seq(from = 0.7, by = 1.2, length.out = nrow(dx.est)),
         x1 = seq(from = 0.7, by = 1.2, length.out = nrow(dx.est)),
         y0 = dx.est$Estimate + dx.est$`Std. Error`,
         y1 = dx.est$Estimate - dx.est$`Std. Error`)
# horizontal bottom
segments(x0 = seq(from = 0.7, by = 1.2, length.out = nrow(dx.est)) - .2,
         x1 = seq(from = 0.7, by = 1.2, length.out = nrow(dx.est)) + .2,
         y0 = dx.est$Estimate + dx.est$`Std. Error`,
         y1 = dx.est$Estimate + dx.est$`Std. Error`)
# horizontal top
segments(x0 = seq(from = 0.7, by = 1.2, length.out = nrow(dx.est)) - .2,
         x1 = seq(from = 0.7, by = 1.2, length.out = nrow(dx.est)) + .2,
         y0 = dx.est$Estimate - dx.est$`Std. Error`,
         y1 = dx.est$Estimate - dx.est$`Std. Error`)


# Plot average gene significance by module
# res_df = fread("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_sumstats_7studiescombined.megaanalysis.txt", h = T)
res_df= fread("~/Google Drive/mac_storage/TWAS/bd_mega/freeze_qcBD_blood_meta_NoGeneFilter.txt",h=T)
net.grab = dx.est$colors
prange = list()
for(i in 1:length(net.grab)){
  gene.grab = colnames(datExpr)[which(moduleColors %in% net.grab[[i]])]
  pvals = res_df$P[res_df$GeneSymbol %in% gene.grab]
  prange[[i]] = data.frame(k = length(gene.grab), 
                           MEp = dx.est$`Pr(>|t|)`[dx.est$colors %in% net.grab[[i]]], 
                           p = mean(pvals), 
                           se = 1.96*(sd(-log10(pvals))/sqrt(length(pvals))))
}
prange = ldply(prange)
prange$logP = -log10(prange$p)
# prange$logSE = -log10(prange$se)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_AvgGeneSignificance-WGCNA.png", res = 300, units = "in", height=5.5,width=8)
barplot(prange$logP,
        las = 2,
        names.arg = paste(dx.est$labels, " (", prange$k,")",sep=""),
        font = 1,
        cex.names = 0.7,
        xlab = NA,
        ylab = "Mean gene significance (95% confidence interval)",
        col = as.character(net.grab), ylim = c(0,1.0))
# vertical bars
segments(x0 = seq(from = 0.7, by = 1.2, length.out = nrow(prange)),
         x1 = seq(from = 0.7, by = 1.2, length.out = nrow(prange)),
         y0 = prange$logP + prange$se,
         y1 = prange$logP - prange$se)
# horizontal bottom
segments(x0 = seq(from = 0.7, by = 1.2, length.out = nrow(prange)) - .2,
         x1 = seq(from = 0.7, by = 1.2, length.out = nrow(prange)) + .2,
         y0 = prange$logP + prange$se,
         y1 = prange$logP + prange$se)
# horizontal top
segments(x0 = seq(from = 0.7, by = 1.2, length.out = nrow(prange)) - .2,
         x1 = seq(from = 0.7, by = 1.2, length.out = nrow(prange)) + .2,
         y0 = prange$logP - prange$se,
         y1 = prange$logP - prange$se)
dev.off()


# Volcano plot of module eigengenes

dx.est$lab_point = ifelse(dx.est$bonf < .05, as.character(dx.est$labels), NA)

png("~/Google Drive/mac_storage/TWAS//bd_mega/data/blood/wgcna_bd.png",res=300,units="in",height=6,width=6)
ggplot(dx.est, aes(x = Estimate, y = logP)) + 
  geom_point(size = 3, pch = 21, col = 'lightgrey', bg = dx.est$colors) + 
  theme_classic() +
  xlab(expression(paste("Effect size (",beta,")"))) +
  ylab(expression(paste("-log"[10],"(P-value)"))) +
  theme(axis.title = element_text(size=15), axis.text = element_text(size=15)) + 
  geom_hline(aes(yintercept = -log10(.05/nrow(dx.est))), col = 'red') + 
  geom_hline(aes(yintercept = -log10(.05)), col = 'orange', linetype=2) +
  geom_text_repel(aes(label = lab_point), size = 4)
dev.off()

# plot gene networks with igraph

# net.grab = dx.est$colors[dx.est$`Pr(>|t|)` < .05]
net.grab = dx.est$colors
# net.grab = net.grab[1:2]

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
adj = as.matrix(adjMat[rownames(adjMat) %in% hubs, colnames(adjMat) %in% hubs])

g <- graph.adjacency(adj,
                       mode="undirected",
                       weighted=T, diag = F) ## Graph information for plotting
  
  hub.order = hubs[match(colnames(adjMat)[colnames(adjMat) %in% hubs], hubs)]
  ME.name = which(colnames(datExpr) %in% hub.order)
  ME.name = paste("ME", moduleLabels[ME.name], sep = "")
  ME.name[ME.name %in% dx.est$labels[dx.est$bonf < .05]] = paste(ME.name[ME.name %in% dx.est$labels[dx.est$bonf < .05]], "*", sep = "")
  font.size = rep(3,length(ME.name))
  font.size[which(grepl("[*]", ME.name))] = 4
  
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
  # values$hub_col[values$name %in% hubs] = "gold"
  values$hub_col[values$colors %in% dx.est$colors[dx.est$bonf < .05]] = "gold"
  
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
  

  netdf = as_data_frame(g, "edges")
  net = network::network(netdf[,1:2], directed = FALSE)
  
  
  
  # plot gene networks
  png(paste("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/plots/BD_networkPlot.png", sep = ""), res= 300,units="in",height=5,width=6)
  
  par(mar=c(0,0,0,0)+.1)
  igraph::plot.igraph(g,
            vertex.size = 10, 
            mark.groups = NULL,
            vertex.shape= "circle",
            vertex.label = ME.name,
            vertex.color = as.character(V(g)$colors),
            vertex.label.cex=0.5,
            vertex.label.color = label.colors,
            edge.width= (E(g)$weight - min(E(g)$weight))/(max(E(g)$weight) - min(E(g)$weight))^2,
            vertex.label.font = 4,
            vertex.frame.color= NA,
            layout=layout.fruchterman.reingold,
            edge.color= "lightblue",
            vertex.label.family="Arial")
  
  dev.off()
  

## MDS plot of gene coexpression network

mds = cmdscale(dist(t(MEs)),2)
colrs = as.character(dx.est$colors[match(colnames(MEs), dx.est$labels)])
mds = data.frame(mds, col = colrs)
mds = mds[!rownames(mds) %in% "ME0", ]

nodeSize = rep(1, nrow(mds))
nodeSize[which(mds$col %in% dx.est$colors[dx.est$bonf < .05])] = 3

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_WGCNA-mdsplot.png", res = 300, units = "in", height = 6, width = 6)
plot(mds[,1:2], xlab= "Component 1", 
     ylab= "Component 2", 
     bg = as.character(mds$col), 
     pch = 21, col = "black", cex = nodeSize)
text(x = mds[,1],y = mds[,2], labels = dx.est$labels[match(mds$col, dx.est$colors)], cex = 0.5)
dev.off()

mds$label = as.character(dx.est$labels[match(mds$col, dx.est$colors)])

nodeSize = -log10(dx.est$P[match(mds$col, dx.est$colors)])
nodeSize

require(ggrepel)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_WGCNA-mdsplot-ggplot2.png", res = 300, units = "in", height = 5, width = 5)
ggplot(mds, aes(x = X1, y = X2)) + 
  geom_point(bg = mds$col, col ='black', size = nodeSize*3, pch = 21) + 
  xlab("Component 1") + ylab("Component 2") + 
  guides(col = F) + theme_classic() + 
  geom_text_repel(aes(x = X1, y = X2, label = mds$label), size = 3) + 
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15))
dev.off()




# detailed connectivity plot of genes in significant modules

net.grab = DIAG$colors[DIAG$bonf < .05]
net.grab = as.character(net.grab[1])

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

# plot gene networks
png(paste("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/plots/BD_networkPlot_rank1.png", sep = ""), res= 300,units="in",height=8,width=8)

par(mar=c(0,0,0,0)+.1)
igraph::plot.igraph(g,
                    vertex.size = 0.1, 
                    mark.groups = NULL,
                    # vertex.shape= node.type,
                    vertex.shape= "circle",
                    # vertex.pie = val.list,
                    # vertex.pie.color=V(g)$pie.color,
                    vertex.color = as.character(V(g)$colors),
                    vertex.label.cex= 6/sqrt(length(V(g))),
                    vertex.label.color = label.colors,
                    edge.width= ((E(g)$weight - min(E(g)$weight))/(max(E(g)$weight) - min(E(g)$weight)))^1.5,
                    vertex.label.font = 3,
                    vertex.frame.color='black',
                    layout=layout.fruchterman.reingold,
                    edge.color= "grey",
                    vertex.label.family="Arial")

dev.off()


# gDel = delete_edges(g, edges = which(E(g)$weight < 0.3))
gDel = g

require(GGally)


netdf = as_data_frame(minimum.spanning.tree(gDel, weights = E(gDel)$weights, algorithm = 'prim'), "edges")
# netdf = as_data_frame(gDel)
net = network::network(netdf[,1:2], directed = FALSE)

downRegulation = res_df$GeneSymbol[res_df$Log2FC < 0]
colorAssignment = ifelse(names(V(gDel)) %in% downRegulation, 'red', 'green')

# plot gene networks
png(paste("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/plots/mst_BD_networkPlot_rank1.png", sep = ""), res= 300,units="in",height=5,width=7)

par(mar=c(0,0,0,0)+.3)
ggnet2(net,
       vjust = -0.6,
       label.size = 2,
       edge.color = 'lightblue', 
       mode = 'fruchtermanreingold', 
       edge.size = 0.5,
       color = colorAssignment,
       size = 2, 
       label = TRUE,
       edge.alpha = 0.5)

dev.off()


# detailed connectivity plot of genes in significant modules

net.grab = as.character(gender.df$colors[gender.df$bonf < .05])
net.grab = net.grab[1]

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

# gDel = delete_edges(g, edges = which(E(g)$weight < 0.3))
gDel = g 
require(GGally)

netdf = as_data_frame(minimum.spanning.tree(gDel, weights = E(gDel)$weights, algorithm = 'prim'), "edges")
net = network::network(netdf[,1:2], directed = FALSE)


# plot gene networks
png(paste("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/plots/mst_BD_networkPlot_rank2.png", sep = ""), res= 300,units="in",height=6,width=6)

ggnet2(net,
       vjust = -0.6,
       label.size = 2.1,
       mode = 'fruchtermanreingold',
       node.color = net.grab,
       size = 1.5, 
       label = TRUE,
       edge.alpha = 0.5)

dev.off()


# detailed connectivity plot of genes in significant modules

net.grab = as.character(gender.df$colors[gender.df$bonf < .05])
net.grab = net.grab[2]

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


# gDel = delete_edges(g, edges = which(E(g)$weight < 0.3))
gDel = g 

require(GGally)

netdf = as_data_frame(minimum.spanning.tree(gDel, weights = E(gDel)$weights, algorithm = 'prim'), "edges")
net = network::network(netdf[,1:2], directed = FALSE)


# plot gene networks
png(paste("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/plots/mst_BD_networkPlot_rank3.png", sep = ""), res= 300,units="in",height=6,width=6)

ggnet2(net,
       vjust = -0.6,
       label.size = 2.1,
       mode = 'fruchtermanreingold',
       node.color = as.character(dx.est$colors[3]),
       size = 1.5, 
       label = TRUE,
       edge.alpha = 0.5)

dev.off()


## Gene set enrichment analysis for module 

enrich = userListEnrichment(geneR = module$symbol, labelR = module$color, useBrainRegionMarkers = TRUE)
pvals = head(enrich$pValues)
pvals = pvals[order(pvals$Pvalues, decreasing = F), ]


entrezid = select(org.Hs.eg.db, keys= as.character(module$symbol), keytype="SYMBOL",columns="ENTREZID")
colnames(entrezid) = c("symbol", "entrezid")
entrezid = merge(entrezid, module, by='symbol')

GOenr = GOenrichmentAnalysis(entrezid$color, entrezCodes = entrezid$entrezid, organism = "human", nBestP = 10)
pval = GOenr$enrichmentP
pval = melt(pval)
pval$FDR = p.adjust(pval$value , "fdr")
pval = pval[order(pval$value, decreasing = F), ]
head(pval[pval$Var1 %in% "midnightblue",])



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

metals = msig[grepl("cadmium|mercury|lead_", ignore.case = TRUE, msig$PATHWAY), ]
table(metals$PATHWAY)
ions = metals[grepl("_ION", metals$PATHWAY), ]
ions = ions[!grepl("CELLULAR", ions$PATHWAY), ]

split = split(ions, ions$PATHWAY)
split_genes = lapply(split, function(x) x$SYMBOL)
Reduce(function(x,y) intersect(x,y), split_genes)

df = ldply(split_genes, cbind)
table(df$`1`, df$.id)

msig$SLIM = NA
msig$SLIM[grepl("metabol|peptidase|lase",ignore.case = T, msig$PATHWAY)]="Metabolism and enzymatic activity"
msig$SLIM[grepl("immun|inflammato|lupus|humoral|innnate|defense|leukocyte|neutrophil|platelet|B cell|lymphoblast|monocyte|blood|inflammation|infereon|cytokine|nf kappa",ignore.case = T, msig$PATHWAY)] = "Immunity"
msig$SLIM[grepl("cell death|apoptosis|apotosi|survival|autophagy|growth|repair|division",ignore.case = T, msig$Set_name)] = "Cell growth, repair, and survival"
msig$SLIM[grepl("chromatin|histone|epigeno|methyltransferase|acetyltransferase",ignore.case = T, msig$PATHWAY)] = "Chromatin regulation"
msig$SLIM[grepl("cytoske|part|lysosome|endosom|golgi|plasma|mitochondri|ribos|nucleu|secreto|vesicl|extracellu",ignore.case = T, msig$PATHWAY)] = "(Extra)cellular structures and organelles"
msig$SLIM[grepl("fungu|bacteri|external|yeast|symbio|stimulus|microbia|parasit|other",ignore.case = T, msig$PATHWAY)] = "Response to external stimuli"
msig$SLIM[grepl("neurotrans|cephalon|cortical|myelin|dopamine|gliogenesis|neurogenesis|astrogenesis|glutama|gaba|calcium|synapse|neuro|glial|astrocyte|axon|brain|behavior|cognition",ignore.case = T, msig$PATHWAY)] = "Neuronal gene sets"
msig$SLIM[grepl("adhesion|motility|movement|integrin|cell-cell",ignore.case = T, msig$PATHWAY)] = "Cell adhesion and contact"
msig$SLIM[is.na(msig$SLIM)] = "Other processes"


sig_mods = unique(as.character(graphDF$colors[graphDF$fdr < .05]))



require(piano)
require(AnnotationDbi)
require(GO.db)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)

genes  = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genes = data.frame(genes)
genes$SYMBOL = select(org.Hs.eg.db, keys=as.character(genes$gene_id), keytype="ENTREZID", columns="SYMBOL")$SYMBOL
goid = select(org.Hs.eg.db, keys=genes$SYMBOL, keytype = "SYMBOL", columns="GO")
goid = goid[goid$ONTOLOGY %in% c("BP"),]
goid = goid[,colnames(goid) %in% c("SYMBOL","GO")]
genes = merge(genes, goid, by="SYMBOL")

genes$PATHWAY = select(GO.db , keys=genes$GO, keytype="GOID", columns="TERM")$TERM

msig = genes

# remove number of characters from long pathway names
ncharcount = nchar(msig$PATHWAY)
toomanychar = which(ncharcount > 120)

substr = substr(x = msig$PATHWAY[toomanychar], 1, 120)
substr = paste(substr, "...",sep="")
msig$PATHWAY[toomanychar] = substr


# number of top ranked gene sets to report per SLIM set
topValue = 1e2

# msig = msig[msig$SYMBOL %in% module$symbol, ]
# msig = msig[grepl("chr|KEGG_|GO_|REACTOME|PANTHER", msig$PATHWAY)]
# msig = msig[!is.na(msig$SYMBOL), ]

msig$in_data = ifelse(msig$SYMBOL %in% colnames(datExpr), "in", "out")

count_int = table(msig$PATHWAY, msig$in_data)
perc_overlap = count_int[,1]/rowSums(count_int)
perc_overlap = perc_overlap[perc_overlap > 0.9]

msig = msig[msig$SYMBOL %in% colnames(datExpr), ] # only look at genes in the datExpr object
msig = msig[msig$PATHWAY %in% names(perc_overlap), ]

length(unique(msig$PATHWAY))

hypStats_save = list()
for( i in 1:length(sig_mods)){
  
  cat("\nRunning hypergeometric test for module:",sig_mods[[i]])
  
  gene_grab = module[module$color %in% sig_mods[[i]], ]

  # filter sets
  set_count = table(as.character(msig$PATHWAY))
  
  list = unique(names(set_count))
  
  msig.keep = msig[msig$PATHWAY %in% list, ]
  msig.keep = msig.keep[,colnames(msig.keep) %in% c("SYMBOL","PATHWAY","GO")]
  
  msig.split = split(msig, msig$PATHWAY)
  
  # calculate hypergeometric input values
  universe = length(unique(module$symbol)) # total # of genes in transcriptome 
  overlaps = lapply(msig.split, function(x) length(intersect(unique(x$SYMBOL), unique(gene_grab$symbol))))
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
  
  # enrichment p-value test
  pvals = phyper(hypStats$Sample_success - 1, 
                 hypStats$Population_success, 
                 hypStats$Population_fails, 
                 hypStats$Sample_size, 
                 lower.tail = FALSE, log.p = FALSE)
  
  hypStats$P = pvals
  hypStats = hypStats[hypStats$Sample_success >= 2, ]
  
  hypStats = hypStats[order(hypStats$P, decreasing = F), ]
  hypStats$FDR = p.adjust(hypStats$P, "fdr")
  hypStats$BONF = p.adjust(hypStats$P, "bonferroni")
  
  hypStats_save[[i]] = hypStats # save sum stats to list obj 
  
}

hypStats_save_df = ldply(hypStats_save)
fes = (hypStats_save_df$Sample_success/hypStats_save_df$Sample_size)/(hypStats_save_df$Population_success/hypStats_save_df$Population)

hypStats_save_df$fes = fes
hypStats_save_df = hypStats_save_df[hypStats_save_df$FDR < .05, ]

dim(hypStats_save_df)

metadata

hypStats_save_df$PATHWAY = gsub("[,]", "_", hypStats_save_df$PATHWAY)
names(hypStats_save_df)[names(hypStats_save_df) %in% "Module"] = 'colors'
hypStats_save_df = merge(metadata, hypStats_save_df,by ='colors')

table(as.character(hypStats_save_df$labels))

as = table(hypStats_save_df$PATHWAY)
as[as > 1]

cor(MEs[,colnames(MEs) %in% c("ME15", "ME17", "ME22", "ME23")])
hypStats_save_df[hypStats_save_df$PATHWAY %in% "apoptotic process", ]

hypStats_save_df = merge(msig[,colnames(msig) %in% c("GO","PATHWAY")], hypStats_save_df, by='PATHWAY')
hypStats_save_df = hypStats_save_df[!duplicated(hypStats_save_df), ]

write.csv(hypStats_save_df, 
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/wgcna_bd_gsea_v2.csv",quote = F, row.names = F)

hypStats_save_df = read.csv("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/wgcna_bd_gsea_v2.csv")

split_set = split(hypStats_save_df, as.character(hypStats_save_df$colors))
split_set = lapply(split_set, function(x) return(x[order(x$P,decreasing = F), ]))
split_set = lapply(split_set, function(x) x[1:5, ])
topGS = ldply(split_set)

top_sets = topGS

topGS = topGS[order(topGS$labels), ]

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/WGCNA_gsea_BD_modules.png",res=300,units="in",height=6.5,width=6)
topY = max(-log10(topGS$P)) * 1.5
ggplot(topGS, aes(x = labels, group = PATHWAY, y = -log10(P))) +
  geom_bar(stat='identity', col = 'black', lwd=0.25, position=position_dodge(0.9), fill = as.character(topGS$colors)) +
  theme_classic() +
  guides(fill = FALSE, size = F) +
  ylab(expression(paste("Gene set enrichment, -log"[10],"(P-value)"))) +
  xlab("Modules") +
  ylim(0, topY) +
  theme(axis.text=element_text(size = 12, colour='black'), axis.title=element_text(size=12, colour='black')) +
  geom_text( aes(y = -log10(P)+1, angle = 90, hjust = 0, label = topGS$PATHWAY), colour='black', size = 3.5, position = position_dodge(0.9)) +
  geom_hline(yintercept = -log10(.05), col = 'red', lwd = 0.3)
dev.off()


# == boxplot of expression for significant module eigengenes
sigMods = dx.est$labels[dx.est$fdr < .05]
sigMods = as.character(sigMods)

modBox = datAll[,colnames(datAll) %in% c(sigMods)]
modBox = resid(lm(as.matrix(modBox) ~ datAll$FACTOR_age + datAll$FACTOR_sex + datAll$FACTOR_studyBatch + datAll$SV1))
modBox = data.frame(Dx = datAll$FACTOR_dx, modBox)
modBox = melt(modBox)
modBox$variable = as.character(modBox$variable)

stats_label = dx.est[dx.est$fdr <= 0.05,]


plot_label <- sprintf("~ beta == %0.3f", stats_label$Estimate)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/WGCNA_boxplot.png",res=300,units="in",height=5,width=5)
ggplot(modBox, aes(x = variable, fill = Dx, y = value)) +
  xlab("Modules") +
  theme_classic() +
  ylab("Eigengene expression (residualized)") +
  scale_fill_manual('Diagnostic group', values = c('lightgrey','dodgerblue4')) +
  geom_hline(yintercept=0, colour='darkgrey', lty = 2) +
  geom_boxplot(notch = TRUE, lwd=0.25, outlier.shape = 5, outlier.color = 'darkgrey', notchwidth = 0.7) +
  theme(axis.text.x=element_text(size = 0, colour='black'), strip.background = element_rect(size = 1.2), panel.border = element_rect(size = 1.1, fill = NA), legend.text = element_text(size = 12), legend.position = 'top', axis.text.y = element_text(size = 12, colour='black'), strip.text.x = element_text(size = 12, colour='black'), axis.title=element_text(size = 12, colour='black')) +
  facet_wrap(~variable, scales = 'free_x') +
  # annotate(geom='text', x = 0.7, y = 0.3, label = paste("B = ",round(stats_label$Estimate,3)))+
  annotate(geom='text', x = 0.7, y = 0.3, label = plot_label, parse = T) +
  annotate(geom='text',x=0.7, y=0.25,   label = paste("p = ",format(stats_label$`Pr(>|t|)`, digits=2), sep=""))
dev.off()



### === Semantic similarity among GO sets
require(GOSemSim)

mods = split(hypStats_save_df, as.character(hypStats_save_df$labels))
d <- godata('org.Hs.eg.db', ont="BP", computeIC=T)

mix = expand.grid(1:length(mods),1:length(mods))

sim_index = list()
for(x in 1:nrow(mix)){
  
  set1 = mods[mix[x,1]]
  set2 = mods[mix[x,2]]
  
  sim = mgoSim(set1[[1]]$GO, set2[[1]]$GO, semData = d, measure = "Wang")
  print(sim)
  sim_index[[x]] = data.frame(Set1 = names(set1), Set2 = names(set2), sim = sim)
  
}
sim_df = ldply(sim_index)
sim_df = reshape2::acast(sim_df, Set1 ~ Set2, value.var ='sim')

require(RColorBrewer)

scaleredblue = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/WGCNA_bd_v_ctrl_semanticSimilarity_genesets.png",res=300,units="in",height=5,width=5.5)
corrplot::corrplot(sim_df, col=scaleredblue, addgrid.col = T,  tl.srt=45, diag = T,
                   tl.col = 'black', method='color', addCoef.col = 'black')
dev.off()

split_set = split(hypStats_save_df, as.character(hypStats_save_df$labels))
split_set = lapply(split_set, function(x) x$PATHWAY)

require(stringdist)
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

require(pheatmap)
png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/similarity_text_gsea.png",res=300,units="in",height=5,width=5.5)
# pheatmap(sim, border_color = 'black', lwd = 1.75)
corrplot::corrplot(sim, method='color', number.cex = 0.85, 
                   order = 'hclust', addCoef.col = 'black', tl.srt = 45, tl.col = 'black', tl.cex = 0.6)
dev.off()


# select top 5 pathways per module

split = split(hypStats_save_df,hypStats_save_df$Module)

for(x in 1:length(split)){
  temp = split[[x]]
  temp = temp[order(temp$P, decreasing = F),]
  temp = temp[1:5, ]
  split[[x]] = temp
}

stats = ldply(split)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/BD_WGCNA_hypergeometricGeneset.png", res = 300,units="in",height=7,width=12)

ggplot(hypStats_save_df, aes(x = PATHWAY, y = -log10(P))) + 
  xlab(NULL) +
  ylab(expression(paste("-log"[10],"(P-value)"))) +
  facet_grid(SLIM ~., scales = "free_y", space = "free_y") + 
  geom_bar(stat = "identity", col = "black", fill = hypStats_save$Module) +
  geom_hline(yintercept = threshold, col = "black", linetype = 2, lwd = 0.5) + 
  theme_classic() + 
  theme(axis.text.y = element_text(size = 6), strip.text.y = element_text(size = 10, angle = 0)) +
  geom_text(aes(label = p_label, y = 20), size= 2) +
  coord_flip()

dev.off()



## == Compare significant modules to schizophrenia blood, schizophrenia brain, bipolar disorder brain, bipolar disorder blood

sig_mods = as.character(graphDF$colors[graphDF$fdr <= 0.05])
sig_mods

scz_mega = fread("~/Google Drive/mac_storage/TWAS/scz_mega/mega.blood.full.csv")
scz_blood = scz_mega$dx_all_ref_merge[scz_mega$qvalues <= 0.05]

scz_brain = fread("~/Google Drive/mac_storage/TWAS/scz_mega/mega.brain.full.csv")
scz_brain = scz_brain$dx_ref2[scz_brain$qvalues <= 0.05]

bd_blood = fread("~/Google Drive/mac_storage/TWAS/bd_mega/freeze_qcBD_blood_meta_Nmin4_LeaveOneOut.txt")
bd_blood = bd_blood$GeneSymbol[bd_blood$FDR <= .05]

bd_brain = readxl::read_excel("~/Google Drive/mac_storage/TWAS/bd_mega/bd_mega_brain_bmc/12888_2013_1399_MOESM2_ESM.xlsx")
bd_brain = data.frame(bd_brain)
bd_brain = bd_brain$genes

set_list = list(scz_brain, scz_blood, bd_brain, bd_blood)
names(set_list) = c("SCZ_brain", "SCZ_blood", "BD_brain", "BD_blood")

stats_df = list()
for(x in 1:length(sig_mods)){
  
  gene_in_mod = as.character(module$symbol[module$color %in% sig_mods[[x]]])
  
  # select genes from set list
  stats = list()
  for(y in 1:length(set_list)){
    
  genes_subject = as.character(unique(set_list[[y]]))
  size = length(genes_subject)
  
  # overlap 
  overlap = intersect(gene_in_mod, genes_subject)
  
  # statistical test
  size_overlap = length(overlap)
  
  enrichPval = phyper(size_overlap - 1,
         length(genes_subject),
         ncol(datExpr) - length(genes_subject),
         length(gene_in_mod), lower.tail = FALSE)
  
  # save stats
  stats[[y]] = data.frame(colors = sig_mods[[x]], Tissue = names(set_list)[[y]], Overlap = size_overlap, SizeModule = length(gene_in_mod), SizeSet = size, Pval = enrichPval)
  
  }
  stats_df[[x]] = ldply(stats)
  
}
stats_all = ldply(stats_df)

write.csv(stats_all, file="~/Google Drive/mac_storage/TWAS/bd_mega/GeneListOverlap_SCZ_BD_brain-blood.csv",quote=F)
