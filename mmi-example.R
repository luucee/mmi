source("mmi.R")

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
mod=strsplit(t$V3,",")
names(mod) = t$V1

#t=read.table("pathways-genes.txt",sep="\t",header=F,stringsAsFactors = F)
#pwg=strsplit(t$V2,", ")
#names(pwg) = t$V1
#t=read.table("tf-pathways.txt",sep="\t",header=F,stringsAsFactors = F)
#tfpw=lapply(strsplit(t$V2,", "),function(x) gsub(" \\(.*\\)","",x))
#names(tfpw) = t$V1
#tfmodulators = tfpw
#tfmodulators = lapply(tfmodulators, function(x) unique(unlist(pwg[x])))

# rimozione probe con poca variabilita'
iqr <- apply(mexp, 1, IQR, na.rm = TRUE)
var.cutoff = quantile(iqr, 0.1, na.rm=TRUE) # percentuale che ne butta fuori
selected = iqr > var.cutoff
M <- mexp[selected, ]
dim(M)
colnames(M) = paste0("S",1:ncol(M))

# prova con tf dell'oracolo
tf = names(mod)
modulators = c("ABL1","BTK","EGFR") 
#unique(unlist(mod))
targets = unique(unlist(trg))
targets=targets[targets %in% rownames(M)]
modulators=modulators[modulators %in% rownames(M)]

# Prendo solo i target che hanno un motif di qualche TF sul promotore
#load("motifRanges.Rdata")
#motifRanges
#motifRanges=motifRanges[motifRanges$qvalue<=0.05]

#sum(tf %in% motifRanges$TF)
#trg=unique(motifRanges$TARGET[motifRanges$TF %in% tf])
#trg=trg[trg %in% rownames(M)]
#length(trg)
#tf = tf[tf %in% motifRanges$TF]
#target = trg

#modulators=sample(modulators,3)
kordering = t(apply(M[modulators,],1,order,decreasing=T))
colnames(kordering) = colnames(M)
library(foreach)
library(doParallel)
cl = makePSOCKcluster(10)
registerDoParallel(cl)
out = mmi(M,tf = tf,target = targets,kordering = kordering,nboot=1000,positiveOnly=F,S=4,sig = 0.01,cl=cl)
stopCluster(cl)

# test prediction performance

miotf="E2F1" #"ETS1"
miomod=modulators[modulators %in% tfmodulators[[miotf]]]
notmod=setdiff(modulators,miomod)

w=out[[miomod[2]]][[miotf]]
w=w[w[,"DIRECTION"]==1,]
w=w[order(w[,"DELTA"]*(1-w[,"FDR"]),decreasing=T),]
w[1:10,]
w=out[[notmod[2]]][[miotf]]
w=w[w[,"DIRECTION"]==1,]
w=w[order(w[,"DELTA"]*(1-w[,"FDR"]),decreasing=T),]
w[1:10,]

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

load("out.Rdata")

g="MYC"
m="EGFR"; m2 = "BTK"; m3="ABL1"
m %in% mod[[g]]
m2 %in% mod[[g]]
m3 %in% mod[[g]]

mm=m
plot(M[mm,kordering[mm,]])

w=out[[mm]]$DELTA1k[g,,]
w[out[[mm]]$PVAL1k[g,,]>0.00000001]=0
w[out[[mm]]$DELTA1k[g,,]<3]=0
w=w[rownames(w) != mm,]
rcol=rep("red",nrow(w))
rcol[rownames(w) %in% trg[[g]]]="blue"
require(gplots)
heatmap.2(w,dendrogram="row",Colv=F,RowSideColors=rcol,density.info="none",trace="none",scale="none")






