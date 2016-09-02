# analisi luciano

setwd("/storage/ocfs2_vol1/PUBLIC/FusionProject")

source("mmi2.R")
load("geDataGBM.RData",verbose=T)
load("gseaResGBM160629.RData",verbose=T)

load("geDataPanglioma.RData",verbose=T)
load("gseaResPanGlioma160629.RData",verbose=T)

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
out = mmi(geData,tf = tf,target = targets,kordering = kordering,nboot=1000,positiveOnly=F,S=3,cl=cl)
stopCluster(cl)
save(out,file="out-GBM.Rdata")
#save(out,file="out-Panglioma.Rdata")

source("mmi2.R")
load("out-GBM.Rdata")
#load("out-Panglioma.Rdata")

# tab output
# PVAL - pvalue del test wilcox tra le MI calcolate nei sample 1..k e all
# MI - mutua informazione dove avviene la regolazione
# MIall - mutua informazione in tutti i sample
# MODdir - direzione della modulazione (+1 la regolazione avviene in High, -1 la regolazione avviene in Low)
# DELTA - rapporto tra MI/Miall
# MOD - gene modulatore
# TF - gene transcription factor
# TRG - target del TF modulato da MOD
stab = summarization(out)

stab$PVALmindy = p.adjust(stab$PVALmindy ,method = "fdr")
stab$PVALkn = p.adjust(stab$PVALkn ,method = "fdr")
stab$PVAL1k = p.adjust(stab$PVAL1k ,method = "fdr")

s = subset(stab,PVALmindy<0.01)
good %in% s$MOD
nrow(s)
s[order(-abs(s$DELTAmindy),s$PVALmindy),][1:10,]
s = subset(stab,PVALkn<0.01)
good %in% s$MOD
nrow(s)
s[order(-abs(s$DELTAkn),s$PVALkn),][1:10,]
s = subset(stab,PVAL1k<0.01)
good %in% s$MOD
nrow(s)
s[order(-abs(s$DELTA1k),s$PVAL1k),][1:10,]

require(xlsx)
file.remove("mod-out.xls")
for(t in tf) {
  s = subset(stab,DELTA>4 & TF==t & FDR<0.00001)
  s = s[order(-s$DELTA,s$FDR),]
  write.xlsx(s,file="mod-out.xls",row.names = F,sheetName = t ,append = T)
}

