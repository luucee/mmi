# analisi luciano

setwd("/storage/ocfs2_vol1/PUBLIC/FusionProject")

source("mmi.R")
load("geDataGBM.RData",verbose=T)
load("gseaResGBM160629.RData",verbose=T)

#load("geDataPanglioma.RData",verbose=T)
#load("gseaResPanGlioma160629.RData",verbose=T)

dim(geData)
tf = c("PPARGC1A","ESRRG")
tf %in% rownames(geData)

targets= unique(unlist(netReg[tf]))
targets = c(targets,"IRS1", "IRS2", "IRS4")
targets= targets[targets %in% rownames(geData)]

modulators=unique(read.table("modulators-luciano.txt",stringsAsFactors = F)[,1])
modulators = modulators[modulators %in% rownames(geData)]


good = c("GNB1", "GNB2", "GNB4")


kordering = t(apply(geData[modulators,],1,order,decreasing=T))
colnames(kordering) = colnames(geData)
require(foreach)
require(doParallel)
cl = makePSOCKcluster(10)
registerDoParallel(cl)
out = mmi(geData,tf = tf,target = targets,kordering = kordering,nboot=1000,positiveOnly=F,S=4,sig = 0.01,cl=cl)
stopCluster(cl)
save(out,file="out-GBM.Rdata")
#save(out,file="out-Panglioma.Rdata")

source("mmi.R")
load("out-GBM.Rdata")
#load("out-Panglioma.Rdata")

# tab output
# PVAL - pvalue del test wilcox tra le MI calcolate nei sample 1..k e all
# MODdir - direzione della modulazione (+1 la regolazione avviene in High, -1 la regolazione avviene in Low)
# MI - media della mutua informazione dove avviene la regolazione
# MIall - media della mutua informazione in tutti i sample
# DELTA - rapporto tra MI/Miall
# MOD - gene modulatore
# TF - gene transcription factor
# TRG - target del TF modulato da MOD
stab = summarization(out)
stab$FDR = p.adjust(stab$PVAL,method="fdr")

# in Paglioma i good sono presenti in 3 liste
s = subset(stab,DELTA>4 & MODdir == -1 & TF=="PPARGC1A" & FDR<0.00001)
good %in% s$MOD # PanGlioma  e GBM
s[order(-s$DELTA,s$FDR),]
s = subset(stab,DELTA>4 & MODdir == 1 & TF=="PPARGC1A" & FDR<0.00001)
good %in% s$MOD
s[order(-s$DELTA,s$FDR),]
s = subset(stab,DELTA>4 & MODdir == -1 & TF=="ESRRG" & FDR<0.00001)
good %in% s$MOD # PanGlioma e GBM
s[order(-s$DELTA,s$FDR),]
s = subset(stab,DELTA>4 & MODdir == 1 & TF=="ESRRG" & FDR<0.00001)
good %in% s$MOD # PanGlioma
s[order(-s$DELTA,s$FDR),]

