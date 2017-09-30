
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
library(Homo.sapiens)
require(org.Hs.eg.db)
require(ggplot2)
require(ggrepel)




# Random effect meta-analysis

## perform meta-analysis of genes (REML model)



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
                                            method = "DL")));
    
    res.std = try(suppressWarnings(metafor::rma(yi = sub_combn$BetaStd,
                                                sei = sub_combn$SEstd, 
                                                verbose = FALSE,
                                                weighted = TRUE, 
                                                method = "DL")))
    
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
res_df = res_df[res_df$Nstudy >= 3, ]

res_df = res_df[order(res_df$P,decreasing=F),]
res_df$FDR = p.adjust(res_df$P, "fdr")
res_df$BONF = p.adjust(res_df$P, "bonferroni")
res_df$GeneSymbol = gsub("[.]", "-", res_df$GeneSymbol)

# Location of genes 

locs <- select(Homo.sapiens, keys=as.character(res_df$GeneSymbol),keytype="SYMBOL", columns=c("ENTREZID", "TXCHROM","TXSTART","TXEND","TXSTRAND"))
locs$WIDTH = abs(locs$TXSTART - locs$TXEND)
locs = locs[order(locs$WIDTH,decreasing = T), ]
locs = locs[!duplicated(locs$SYMBOL), ]
locs$LOC = paste(locs$TXCHROM, ":",locs$TXSTART,"-",locs$TXEND,"(",locs$TXSTRAND,")",sep="")
locs = locs[!is.na(locs$TXCHROM),]
loc_df = locs[,c("SYMBOL","LOC", "ENTREZID")]
colnames(loc_df)[1] = "GeneSymbol"

# Combine results with chromosome map

res_df = merge(loc_df, res_df, by="GeneSymbol",all.y=T)

fwrite(data.table(res_df), 
       sep = "\t",
       file = "~/Google Drive/mac_storage/TWAS/bd_mega/qcBD_blood_meta_Nmin3_LeaveOneOut.txt", quote  = F, row.names = F)

# qq-plot


# association p-values
assoc = data.frame(P = res_df$P, source='DGE Meta-analysis')

observed <- sort(assoc$P)
lobs <- -(log10(as.numeric(observed)))

expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))


assoc = assoc[order(assoc$P, decreasing = F), ]
assoc$lobs = lobs
assoc$lexp = lexp

ci = .95
N = length(assoc$P)
observed = -log10(sort(assoc$P))
expected = -log10(1:N / N)
clower   = -log10(qbeta(ci,     1:N, N - 1:N + 1))
cupper   = -log10(qbeta(1 - ci, 1:N, N - 1:N + 1))
assoc$clower = clower
assoc$cupper = cupper

# heterogeneity p-values
het = data.frame(P = res_df$HetP, source='Heterogeneity')

observed <- sort(het$P)
lobs <- -(log10(as.numeric(observed)))

# qqplot statistics
chisq1 <- qchisq(1-observed, 1)
medianchi = median(chisq1)
lambdaout = medianchi/.454


expected <- c(1:length(observed))
lexp <- -(log10(expected / (length(expected)+1)))


het = het[order(het$P, decreasing = F), ]
het$lobs = lobs
het$lexp = lexp

ci = .95
N = length(het$P)
observed = -log10(sort(het$P))
expected = -log10(1:N / N)
clower   = -log10(qbeta(ci,     1:N, N - 1:N + 1))
cupper   = -log10(qbeta(1 - ci, 1:N, N - 1:N + 1))
het$clower = clower
het$cupper = cupper

set = ldply(list(assoc,het))

mclab = substitute(paste("Median ", chi^2, "=", MC, " ", lambda , "=", LD), list(MC=format(medianchi, digits = 3), LD=format(lambdaout,digits=3)))

qplot = ggplot(assoc, aes(x = lexp, y = lobs)) +
  geom_point(colour="black", fill= 'dodgerblue', shape = 21, size = 2, stroke = 0.3) + 
  geom_abline(slope = 1, intercept = 0, col = 'black', lwd = 0.3) +
  theme_classic() + 
  ylab(expression(paste("Observed -log"[10]," P-value"))) +
  xlab(expression(paste("Expected -log"[10]," P-value"))) +
  scale_fill_discrete(NULL) +
  geom_line(aes(x = lexp, y = clower), colour ='grey' ,lwd = 0.75) +
  geom_line(aes(x = lexp, y = cupper), colour = 'grey', lwd = 0.75) +
  geom_ribbon(aes(x = lexp, ymin = clower, ymax = cupper), fill="grey", alpha="0.2") +
  ggtitle(mclab)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/qcBD_qqplot.png",res=300,units="in",height=6,width=6)
print(qplot)
dev.off()


# Manhattan-style plot


## Manhattan plot

trim = res_df[!is.na(res_df$LOC),]
locs = strsplit(trim$LOC, "[:]")
chr = lapply(locs, function(x) x[[1]])
bp = lapply(locs, function(x) x[[2]])
pos = data.frame(CHR=unlist(chr),BP=unlist(bp))
pos$GeneSymbol = trim$GeneSymbol
pos$bonf = trim$BONF
pos$fdr = trim$FDR
pos$P = trim$P
pos$BP = gsub("[(-)]", "", pos$BP)
pos$BP = gsub("[(+)]", "", pos$BP)
pos$START = unlist(lapply(strsplit(pos$BP, "[-]"), function(x)x[[1]]))
pos$END = unlist(lapply(strsplit(pos$BP, "[-]"), function(x)x[[2]]))
pos$CHR = gsub("chr","",pos$CHR)
pos$CHR[pos$CHR %in% "X"] = 23
pos$CHR[pos$CHR %in% "Y"] = 24

col = data.frame(CHR = 1:24, col = c("forestgreen", "lightgreen"))

sub = merge(pos, col, by = "CHR")
sub$CHR = as.integer(sub$CHR)
sub$START = as.integer(sub$START)
sub$pos = NA
sub$pos = ifelse(sub$CHR  == 1, as.integer(sub$START), NA)

median_pos = list()
chr_grab = unique(col$CHR)
for( i in 1:length(chr_grab)){
  
  if(i == 1){
    prior_max = min(sub$pos[sub$CHR == 1]) - 1
  }
  
  if(i > 1){
    k = i - 1
    prior_max = max(sub$pos[sub$CHR == k])
  }
  
  chr_seq = sub$START[sub$CHR == chr_grab[[i]]]
  true_min = min(chr_seq)
  diff = chr_seq - true_min
  new_chr_seq = diff + prior_max + 1
  
  sub$pos[sub$CHR == i] <- new_chr_seq
  
  median_pos[[i]] =  (min(sub$pos[sub$CHR %in% chr_grab[[i]]]) + max(sub$pos[sub$CHR %in% chr_grab[[i]]]))/2
}

names(median_pos) = 1:24

median_pos = ldply(median_pos)
median_pos$.id[median_pos$.id %in% c(23,24)] = c("X","Y")


sub = sub[order(sub$CHR, sub$pos), ]
sub$SYMBOL = ifelse(sub$fdr < .05, sub$GeneSymbol, NA)



png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/qcBDvCT_ggplot-manhattan.png", res = 300, units = "in", height = 6, width = 11)
# manhattan plot
ggplot(sub, aes(x = sub$pos, y=-log10(sub$P))) + 
  geom_point(size = 0.7, col  = as.character(sub$col)) + 
  ylim(min = 0, max = 1.2*max(-log10(sub$P))) + 
  xlab("Genomic coordinate") + 
  theme_classic() +
  ylab(expression(paste("-log"[10],"(P-value)"))) +
  scale_x_continuous(name="Genomic coordinate", breaks=median_pos$V1, labels=median_pos$.id) +
  geom_hline(yintercept = -log10(.05/nrow(sub)), col = "red", lwd = 0.5, linetype = "dashed") +
  geom_text_repel(aes(label = sub$SYMBOL), fontface = 'italic', size = 3, col = "black")
# end plot
dev.off()



# volcano plot


# 6. volcano plots for meta-analysis and mega-analysis
psize = -log10(res_df$P)
psize = (psize - min(psize))/(max(psize) - min(psize)) + 0.5
col = rep("darkgrey", nrow(res_df))
col[which(abs(res_df$Log2FC) > 0.2)] = "dodgerblue3"
col[which(res_df$FDR < .05)] = "orange"
col[which(res_df$BONF < .05)] = "darkred"

res_df$vLabel = ifelse(abs(res_df$Log2FC) > .2 | res_df$FDR < .05, res_df$GeneSymbol, NA)

png("~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/qcBD_meta-volcano.png",res=300,units="in",height = 6, width = 6)
ggplot(res_df, aes(x = Log2FC, y = -log10(P))) + 
  geom_point(size = psize, col = col) + 
  theme_bw() + 
  geom_hline(yintercept = -log10(.05), col = 'red', linetype = "dashed") + 
  xlab( expression(paste("BD Meta-anaylsis DGE (log"[2]," fold-change)"))) + 
  ylab( expression(paste("DGE significance, -log"[10],"(P-value)"))) + 
  geom_hline(aes(yintercept = -log10(.05/nrow(res_df))), col = 'green', lwd = 0.3) +
  geom_vline(aes(xintercept = 0), col = 'grey',lwd=0.5, linetype='dashed') +
  ggrepel::geom_text_repel(aes(label = res_df$vLabel), size = 2, fontface = 3)
dev.off()
