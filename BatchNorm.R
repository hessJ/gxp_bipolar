

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
require(illuminaHumanv4.db)
require(hgu95av2.db)
require(hgu133plus2.db)
require()
require(hgu133a.db)
require(hgu133b.db)

# list directories containing microarray data 
# Compatible with Affymetrix CEL (gene chip or exon array) and Illumina txt file (ProbeID, Sample, Detection Pval)
dir = list.dirs(path="~/Google Drive/mac_storage/TWAS/bd_mega/data/blood",full.names = T,recursive = F)

# make QC plot folder
dircreate = dir.exists(path = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/QCplots")
if(dircreate == F){ dir.create(path = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/QCplots")}
QCfolder = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/QCplots"

# make normalized expression folder
dircreate = dir.exists(path = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/normalized_data")
if(dircreate == F){ dir.create(path = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/normalized_data")}
NormOut = "~/Google Drive/mac_storage/TWAS/bd_mega/data/blood/normalized_data"


# Serial processing
# Import microarray (RNA expression) data from each study
# 1. Normalize 
# 2. Plot normalized data to new folder
# 3. Write normalized data to new folder
# AFFY: GC-RMA or RMA w/ quantile normalization (log2 transformed)
# ILMN: Background correction based on Detection Pval columns (log2 transformed)

for( i in 1:length(dir)){
  
  cat("\rImporting data from:",dir[[i]])
  
  # List cel files (if present)
  CEL = NULL
  NON = NULL
  
  CEL = list.files(dir[[i]], full.names =  T, pattern = 'CEL')
  
  if(length(CEL) > 1){
    cat("\n   Detected CEL files...")
    
    readCel = read.celfiles(filenames = CEL)
    
    # is this series an exon array? If so, GCRMA will fail
    exonArrayChecker = readCel@annotation[grepl("huex|st",readCel@annotation)]
    
    if(length(exonArrayChecker) < 1){
      # the road to GC-RMA, quantile normalized gene chip data
      Batch = ReadAffy(filenames = CEL)
      RMA = gcrma::gcrma(Batch, normalize = T); # GC-RMA  (gene chip)
      Exprs = exprs(RMA)
      Exprs = data.frame(PROBEID = rownames(Exprs), Exprs)
      } else {
        
        # the road to RMA, quantiled normalized exon array data
      RMA = affy::rma(readCel, normalize = T, target = "probeset");
      
      featureData(RMA) <- getNetAffx(RMA, "probeset")
      annot = pData(featureData(RMA))
      annot = annot[,c("probesetid","transcriptclusterid")]
      
      fwrite(annot, file = paste(NormOut,"/",basename(dir[[i]]),"_PROBEID.txt", sep =""),
             quote = F, row.names = F, sep= "\t")
      
      Exprs = exprs(RMA)
      
      Exprs = data.frame(PROBEID = rownames(Exprs), Exprs)
      
      } # standard RMA (exon array)
    
  }
    
    # Possibly a pre-made expression matrix (common with Illumina data on GEO)
    if(length(CEL) < 1){
      
      NON = list.files(dir[[i]], full.names = T, pattern = "non-normalized")
      NON = NON[!grepl(".gz", NON)]
      if(length(NON) == 0) next
      cat("\n   Detected ILMN matrix...")
      NONread = fread(NON, h=T, sep=  "\t")
      exprnames = colnames(NONread)[!grepl("Probe|ID_REF|Detection", colnames(NONread))]
      
      colnames(NONread)[colnames(NONread) %in% exprnames] = "avgExpr"
      
      fwrite(NONread, 
             file = paste(dir[[i]],"/ILMN_format.txt",sep=""),
             quote = F, row.names = F, sep = "\t")
      
      # read illumina data 
      idata <- read.ilmn(paste(dir[[i]],"/ILMN_format.txt",sep=""), 
                         other.columns="Detection Pval",
                         expr = "avgExpr",
                         probeid = colnames(NONread)[[1]])
    
     
      # background correction using Detection p-value column
      y <- neqc(idata, detection.p = "Detection Pval")
      # expression object
      Exprs = data.frame(y$E)
      colnames(Exprs) = exprnames
      Exprs = data.frame(PROBEID = rownames(Exprs), Exprs)
      
    }
    
    # Boxplot of expression
    if(ncol(Exprs) > 1){
      cat("\nGenerating plot for normalized data")
    png(paste(QCfolder,"/BOXPLOTexprnormed_",basename(dir[[i]]),".png", sep =""), res=300,units="in",height = 6, width = 10)
    boxplot(Exprs[,-1], 
              col = 'tan', 
              las = 2, 
              xaxt = "n",
              outline = F, 
              ylab = expression(paste("log"[2]," expression (quantile normalized)")),
              main = paste("Expression data from:",basename(dir[[i]])),sep="")
    axis(side = 1, 
         labels = paste("S",1:(ncol(Exprs)-1),sep=""), 
         at = 1:(ncol(Exprs)-1), las = 2,cex.axis = 0.6)
      dev.off()
    }
  
  if(ncol(Exprs) > 0){
  # Write normalized expression (probe level) to data file
    cat("\nWriting normalized expression data")
  fwrite(Exprs,
         file = paste(NormOut,"/",basename(dir[[i]]),"_normalizedProbes.txt", sep =""),
         quote = F, row.names = F, sep= "\t")}
  
  if(length(CEL) == 0 & length(NON) == 0) {stop("No expression files detected! Data may be placed in wrong directory.")}
    
}


  
