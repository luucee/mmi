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
t=read.table("pathways-genes.txt",sep="\t",header=F,stringsAsFactors = F)
pwg=strsplit(t$V2,", ")
names(pwg) = t$V1
t=read.table("tf-pathways.txt",sep="\t",header=F,stringsAsFactors = F)
tfpw=lapply(strsplit(t$V2,", "),function(x) gsub(" \\(.*\\)","",x))
names(tfpw) = t$V1
tfmodulators = tfpw
tfmodulators = lapply(tfmodulators, function(x) unique(unlist(pwg[x])))

# rimozione probe con poca variabilita'
iqr <- apply(mexp, 1, IQR, na.rm = TRUE)
var.cutoff = quantile(iqr, 0.1, na.rm=TRUE) # percentuale che ne butta fuori
selected = iqr > var.cutoff
M <- mexp[selected, ]
dim(M)
colnames(M) = paste0("S",1:ncol(M))

# verifica che i geni del goldstandar siano presenti nel dataset
unique(unlist(tfmodulators)) %in% rownames(M)
unique(unlist(tfmodulators)) %in% modulators
names(tfmodulators) %in% rownames(M)

modulators = modulators[modulators %in% rownames(M)]
tf = tf[tf %in% rownames(M)]

target = setdiff(rownames(M),c(tf,modulators))

# Prendo solo i target che hanno un motif di qualche TF sul promotore
load("motifRanges.Rdata")
motifRanges
motifRanges=motifRanges[motifRanges$qvalue<=0.05]

sum(tf %in% motifRanges$TF)
trg=unique(motifRanges$TARGET[motifRanges$TF %in% tf])
trg=trg[trg %in% rownames(M)]
length(trg)
tf = tf[tf %in% motifRanges$TF]
target = trg

#modulators=sample(modulators,3)
kordering = t(apply(M[modulators,],1,order,decreasing=T))
colnames(kordering) = colnames(M)
library(foreach)
library(doParallel)
cl<-makeCluster(ncore)
registerDoParallel(cl)
out = mmi(M,tf = tf,target = target,kordering = kordering,nboot=100,positiveOnly=F,S=200,sig = 0.01,ncore=50,cl=cl)
stopCluster(cl)

# test prediction performance

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


