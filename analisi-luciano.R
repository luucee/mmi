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

fus<-c( "TCGA-CF-A3MF",   "TCGA-CF-A3MG",   "TCGA-CF-A3MH",   "TCGA-DS-A3LQ",   "TCGA-WL-A834",  
        "TCGA-06-6388" ,  "TCGA-06-6390",   "TCGA-12-0820",   "TCGA-12-0826",  
        "TCGA-19-5958"  , "TCGA-27-1835",   "TCGA-27-1835-d", "TCGA-28-2506",   "TCGA-74-6578" , 
        "TCGA-76-4925"   ,"TCGA-76-4925-d", "TCGA-CR-6473",   "TCGA-CV-7100",   "TCGA-FG-7643" , 
        "TCGA-P5-A72U"   ,"TCGA-22-4607",   "TCGA-34-2608",   "TCGA-39-5024",   "TCGA-66-2786"  )

modulators=unique(read.table("modulators-luciano.txt",stringsAsFactors = F)[,1])
modulators = modulators[modulators %in% rownames(geData)]

#load("motifRanges.Rdata")
#motifRanges
#motifRanges=motifRanges[motifRanges$qvalue<=0.05]
#sum(tf %in% motifRanges$TF)
#targets=unique(motifRanges$TARGET[motifRanges$TF %in% tf])

#targets= unique(unlist(netReg[tf]))
targets= setdiff(rownames(geData),c(tf,modulators))
#targets = c(targets,"IRS1", "IRS2", "IRS4")
targets= targets[targets %in% rownames(geData)]

good = c("GNB1", "GNB2", "GNB4")


#kordering = t(apply(geData[modulators,],1,order,decreasing=T))
#colnames(kordering) = colnames(geData)
require(foreach)
require(doParallel)
cl = makePSOCKcluster(10)
registerDoParallel(cl)

# signif dei delta
out = NULL
for(mod in modulators) {
  ptm = proc.time()[3]
  out=rbind(out,mindy2(geData,mod=mod,tf=tf,target = targets,nboot = 10000,nbins=3,h=1,siglev=0.05,method="pearson")) # equiv mindy (nbins=3 h=1)
  print(paste0(mod," took ",proc.time()[3]-ptm," sec."))
}
save(out,file="out-Panglioma.Rdata")
save(out,file="out-GBM.Rdata")
load("out-Panglioma.Rdata")
load("out-GBM.Rdata")
load("out-bins4-GBM.Rdata")
load("out-alltargets-GBM.Rdata")

#out$PVAL = p.adjust(out$PVAL,method = "fdr")
out.sig = subset(out,(DELTA>0.4 | DELTA< -0.4) & PVAL<0.01)
nrow(out.sig)
ntrg = aggregate(out.sig$TRG,by=list(out.sig$MOD,out.sig$TF),function(x) length(unique(x)))
ntrg[,4] = 0
ntrg[,5] = aggregate(out.sig$DELTA,by=list(out.sig$MOD,out.sig$TF),mean)[,3]


# significativita dei num di target per ogni mod-tf
outperm = foreach(b=1:1000, .combine = "rbind") %dopar% {
    outp = mindy2(geData,mod=modulators[1],tf=tf,target = targets,nboot = 80,nbins=5,h=1,perm=T,siglev=0.05,method="pearson")
    outp$B=b
    outp
}
save(outperm,file="outperm-Panglioma.Rdata")
save(outperm,file="outperm-GBM.Rdata")
load("outperm-Panglioma.Rdata")
load("outperm-GBM.Rdata")
load("outperm-bins4-GBM.Rdata")

for(m in unique(as.character(outperm$MOD))) {
  outperm$PVAL[outperm$MOD==m] = p.adjust(outperm$PVAL[outperm$MOD==m],method = "fdr")
}
outperm.sig = subset(outperm,PVAL<0.1)

ntrg.perm = aggregate(outperm.sig$TRG,by=list(outperm.sig$MOD,outperm.sig$TF,outperm.sig$B),function(x) length(unique(x)))
for(i in 1:nrow(ntrg)) {
  mod=as.character(ntrg[i,1])
  tf=as.character(ntrg[i,2])
  ntrg[i,4] = sum(ntrg[i,3] < ntrg.perm[ntrg.perm$MOD==mod & ntrg.perm$TF==tf,4])/1000
}

colnames(ntrg) = c("MOD","TF","N.TRG","PVAL")

require(xlsx)
write.xlsx(ntrg,file="modlist.xls",row.names = F,sheetName = "Panglioma" ,append = T)

mexp = cut.outliers(geData)

source("mmi2.R")
mytf="PPARGC1A" #"ESRRG"
for (mod in as.character(ntrg$MOD[ntrg$TF==mytf])) {
  cat(mod,"\n")
  f=out.sig[out.sig$MOD==mod & out.sig$TF==mytf,]
  trg = unique(as.character(f$TRG))
  cat(length(trg),"\n")
  if (length(trg)>=30) {
    pdf(paste0(mytf,"-",mod,"-high.pdf"))
    plot.mod(mexp,mod,mytf,target = trg,
             nettarget = trg,fus=fus,high=T)
    dev.off()
    pdf(paste0(mytf,"-",mod,"-low.pdf"))
    plot.mod(mexp,mod,mytf,target = trg,
             nettarget = trg,fus=fus,high=F)
    dev.off()
  }
  
}












mdelta=format(mean(t$DELTA),digits = 2)
targets = as.character(t$TRG)
ntargets = unique(unlist(netReg[[tf]]))






#layout(matrix(1:2,nrow=1,ncol=2,byrow = T))


#out = mmi(geData,tf = tf,target = targets,kordering = kordering,nboot=1000,positiveOnly=F,S=3,cl=cl)

out = mindy.perm(geData,tflist = tf,target = targets,kordering = kordering,nboot=1000,S=3,cl=cl)

otab= otab[order(otab[,4]),]
write.xlsx(otab,file="mindy-GBM.xls",row.names = F,sheetName = t ,append = F)

load("geDataGBM.RData")

load("otarget-GBM.Rdata")
s = subset(otarget$OUT,PVAL<0.01)
nrow(s)

modl=c("ATP1A1","DSG2","NEXN","LRRFIP1","PLOD3","TBC1D4","CDC42BPA","PIN4","RANBP10","MLPH","PRKCDBP","DRG1")
mod="LDHB"; mod ="TBC1D4"
tf="PPARGC1A"
mm = mindy2(geData,mod=mod,tf=tf,target = targets,nboot = 30,nbins=3,h=1) # equiv mindy

#heatmap.2(geData[targets,],col=redgreen(100),density.info="none",trace="none",scale="none")


out = mindy(geData,tf = tf,target = targets,kordering = kordering,nboot=1000,S=3,cl=cl)
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
stab = mindy.sum(out)

stab$PVAL = p.adjust(stab$PVAL ,method = "fdr")
stab$PVALkn = p.adjust(stab$PVALkn ,method = "fdr")
stab$PVAL1k = p.adjust(stab$PVAL1k ,method = "fdr")

s = subset(stab,PVAL<0.01)
good %in% s$MOD
nrow(s)
s=s[order(-abs(s$DELTAmindy),s$PVALmindy),]
cumsum(s$MOD %in% good)

s = subset(stab,PVALkn<0.01)
s=s[order(-s$DELTAkn,s$PVALkn),]
nrow(s)
cumsum(s$MOD %in% good)

s = subset(stab,PVAL1k<0.01)
s=s[order(-s$DELTA1k,s$PVAL1k),]
nrow(s)
cumsum(s$MOD %in% good)

require(xlsx)
file.remove("mod-out.xls")
for(t in tf) {
  s = subset(stab,DELTA>4 & TF==t & FDR<0.00001)
  s = s[order(-s$DELTA,s$FDR),]
  write.xlsx(s,file="mod-out.xls",row.names = F,sheetName = t ,append = T)
}

