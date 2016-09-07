source("mmi2.R")

# Dataset lavoro Gambardella et al. BMC Bioinformatics (2015) 16:279
#require(R.matlab)
#d=readMat("ATLAS_dataset.mat",verbose = T)
#mexp=d$expression
#allgenes=unlist(d$map.platform)
#modulators=unlist(d$map.K)
#tf=c("CDX2","E2F1","ELK1","ETS1","GATA1","GATA2","MYC","SMAD3","SMAD4","STAT1","STAT3","STAT6","TCF4","TP53")
#rownames(mexp) = allgenes
#save(tf,modulators,allgenes,mexp,file="ATLAS_dataset.Rdata")
load(file="ATLAS_dataset.Rdata",verbose=T)

# Recupero goldstandard dalle tabelle del paper
t = read.table("tf-modulators-targets.txt",sep="\t",stringsAsFactors = F)
t$V2 = gsub("[ |\\n|\\r]","",t$V2)
t$V3 = gsub("[ |\\n|\\r]","",t$V3)
trg=strsplit(t$V2,",")
names(trg) = t$V1
modlist=strsplit(t$V3,",")
names(modlist) = t$V1

#t=read.table("pathways-genes.txt",sep="\t",header=F,stringsAsFactors = F)
#pwg=strsplit(t$V2,", ")
#names(pwg) = t$V1
#t=read.table("tf-pathways.txt",sep="\t",header=F,stringsAsFactors = F)
#tfpw=lapply(strsplit(t$V2,", "),function(x) gsub(" \\(.*\\)","",x))
#names(tfpw) = t$V1
#tfmodulators = tfpw
#tfmodulators = lapply(tfmodulators, function(x) unique(unlist(pwg[x])))

# rimozione probe con poca variabilita'
#iqr <- apply(mexp, 1, IQR, na.rm = TRUE)
#var.cutoff = quantile(iqr, 0.1, na.rm=TRUE) # percentuale che ne butta fuori
#selected = iqr > var.cutoff
#M <- mexp[selected, ]
#M <- mexp
#dim(M)

# prova con tf dell'oracolo
tf = names(trg)
modulators = unique(unlist(modlist))
targets = unique(unlist(trg))
targets=targets[targets %in% rownames(mexp)]
modulators=modulators[modulators %in% rownames(mexp)]
for (tt in names(modlist)) { # per l'oracolo solo quelli che esistono nei dati
  modlist[[tt]] = modlist[[tt]][modlist[[tt]] %in% modulators]
}
oracle=NULL
for(tt in names(modlist)) {
  oracle=rbind(oracle,data.frame(TF=tt,MOD=modlist[[tt]]))
}
oracle=paste(oracle$MOD,oracle$TF)

source("mmi2.R")
# signif dei delta
out = NULL
for(mod in modulators) {
  ptm = proc.time()[3]
  out=rbind(out,mindy2(mexp,mod=mod,tf=tf,target = targets,nboot = 1000,nbins=3,h=1,siglev=0.01,method="pearson")) # equiv mindy (nbins=3 h=1)
  print(paste0(mod," took ",proc.time()[3]-ptm," sec."))
}
save(out,file="out-ATLAS.Rdata")
load("out-ATLAS.Rdata")
out=subset(out,PVAL<0.01)
out$PVAL=p.adjust(out$PVAL,method = "fdr")
require(ROCR)
layout(matrix(1:3,nrow=1,ncol=3))
for (tt in tf) {
  x = out[out$TF==tt,]
  f=paste(x$MOD,x$TF)
  pred = prediction(abs(x$DELTA),f %in% oracle)
  perf = performance(pred,measure="prec",x.measure="rec")
  plot(perf,ylim=c(0,1),xlim=c(0,0.1),main=tt)
}



perf = performance(pred,measure="tpr",x.measure="fpr")
plot(perf)
perf = performance(pred,measure="auc")

out.sig = subset(out,(DELTA>0.3 | DELTA< -0.3) & PVAL<0.05)
nrow(out.sig)
ntrg = aggregate(out.sig$TRG,by=list(out.sig$MOD,out.sig$TF),function(x) length(unique(x)))
ntrg[,4] = 0
ntrg[,5] = aggregate(out.sig$DELTA,by=list(out.sig$MOD,out.sig$TF),mean)[,3]


# significativita dei num di target per ogni mod-tf
outperm = foreach(b=1:100, .combine = "rbind") %dopar% {
  outp = mindy2(M,mod=modulators[1],tf=tf,target = targets,nboot = 1000,nbins=3,h=1,perm=T,siglev="none",method="pearson")
  outp$B=b
  outp
}
save(outperm,file="outperm-ATLAS.Rdata")
load("outperm-ATLAS.Rdata")

outperm.sig = subset(outperm,PVAL<0.1)

ntrg.perm = aggregate(outperm.sig$TRG,by=list(outperm.sig$MOD,outperm.sig$TF,outperm.sig$B),function(x) length(unique(x)))
for(i in 1:nrow(ntrg)) {
  mod=as.character(ntrg[i,1])
  tf=as.character(ntrg[i,2])
  ntrg[i,4] = sum(ntrg[i,3] < ntrg.perm[ntrg.perm$MOD==mod & ntrg.perm$TF==tf,4])/1000
}

colnames(ntrg) = c("MOD","TF","N.TRG","PVAL")



# costruzione oracolo lista inversa TF -> modulator
tfpred =list()
for(i in names(out)) {
  for (j in names(out[[i]])) {
    if (nrow(out[[i]][[j]])>0) {
      # considero come score del modulatore il minimo del pvalue tra i target regolati dal tf
      pv=min(out[[i]][[j]][,"FDR"])
      names(pv) = i
      tfpred[[j]] = c(tfpred[[j]],pv)
    }
  }
}

require(ROCR)
for(i in names(tfpred)) {
  pred = prediction(1-tfpred[[i]],names(tfpred[[i]]) %in% tfmodulators[[i]])
  perf = performance(pred,measure="prec",x.measure="rec")
  plot(perf)
}
source("mmi.R")
load("outl.Rdata")

s = summarization(out)
s$FDR = p.adjust(s$PVAL,method="fdr")
s = subset(s,FDR<0.01)
s[order(-s$DELTA,s$FDR),][1:30,]

g="TP53"
# STAT3, MYC, TP53
m="EGFR"; m2 = "BTK"; m3="ABL1"
m %in% mod[[g]]
m2 %in% mod[[g]]
m3 %in% mod[[g]]

mm=m2
plot(M[mm,kordering[mm,]])

w=out[[mm]]$DELTA1k[g,,]
max(w)
w[out[[mm]]$PVALkn[g,,]>0.00000001]=0
w[out[[mm]]$DELTA1k[g,,]<2]=0
w=w[rownames(w) != mm,]
w = (out[[mm]]$DELTAkn[g,,]>2 & out[[mm]]$PVALkn[g,,]<0.01) | (out[[mm]]$DELTA1k[g,,]>2 & out[[mm]]$PVAL1k[g,,]<0.01)
rcol=rep("red",nrow(w))
rcol[rownames(w) %in% trg[[g]]]="blue"
require(gplots)
heatmap.2(w[rownames(w) %in% trg[[g]],]*1,dendrogram="row",Colv=F,#RowSideColors=rcol,
          density.info="none",trace="none",scale="none")


h=c("RPL22","TP53","BTK")
heatmap(M[rownames(M) %in% h,])



