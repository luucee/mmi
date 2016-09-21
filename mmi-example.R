
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
mexp=mexp[rownames(mexp)!="NA",]


# recupero gold standard paper cindy (chiesto all'autore F. Giorgi)
library(org.Hs.eg.db)
library(annotate)
d=read.table("federico-giorgi/unionset.txt",sep="\t",header=F)
g1=getSYMBOL(as.character(d[,1]),data="org.Hs.eg")
g2=getSYMBOL(as.character(d[,2]),data="org.Hs.eg")

require(xlsx)
ts3.mod = as.character(read.xlsx("federico-giorgi/Table_S3.xlsx",sheetName = "Modulators")[,1])
ts3.tf = as.character(read.xlsx("federico-giorgi/Table_S3.xlsx",sheetName = "TFs")[,1])
length(ts3.mod)

tf1 = which(g1 %in% ts3.tf)
tf2 = which(g2 %in% ts3.tf)

mod = c(g2[tf1],g1[tf2])
tf = c(g1[tf1],g2[tf2])

mod.tf = data.frame(MOD=mod,TF=tf,stringsAsFactors = F)
nrow(mod.tf)
mod.tf = subset(mod.tf,MOD %in% ts3.mod & TF %in% ts3.tf)
nrow(mod.tf)
length(unique(paste(mod.tf$MOD,mod.tf$TF)))

mod.tf = subset(mod.tf,MOD %in% rownames(mexp) & TF %in% rownames(mexp))
nrow(mod.tf)
oracle=paste(mod.tf$MOD,mod.tf$TF)

# Recupero goldstandard dalle tabelle del paper Di Bernardo (dubbi)
#t = read.table("tf-modulators-targets.txt",sep="\t",stringsAsFactors = F)
#t$V2 = gsub("[ |\\n|\\r]","",t$V2)
#t$V3 = gsub("[ |\\n|\\r]","",t$V3)
#trg=strsplit(t$V2,",")
#names(trg) = t$V1
#modlist=strsplit(t$V3,",")
#names(modlist) = t$V1


#t=read.table("pathways-genes.txt",sep="\t",header=F,stringsAsFactors = F)
#pwg=strsplit(t$V2,", ")
#names(pwg) = t$V1
#t=read.table("tf-pathways.txt",sep="\t",header=F,stringsAsFactors = F)
#tfpw=lapply(strsplit(t$V2,", "),function(x) gsub(" \\(.*\\)","",x))
#names(tfpw) = t$V1
#tfmodulators = tfpw
#tfmodulators = lapply(tfmodulators, function(x) unique(unlist(pwg[x])))

source("mmi2.R")

# rimozione relazione mod-tf dirette
m = aracne2(mexp,mod.tf$MOD,mod.tf$TF)

# rimozione probe con poca variabilita'
#iqr <- apply(mexp, 1, IQR, na.rm = TRUE)
#var.cutoff = quantile(iqr, 0.1, na.rm=TRUE) # percentuale che ne butta fuori
#selected = iqr > var.cutoff
#M <- mexp[selected, ]
#M <- mexp
#dim(M)

# prova con tf dell'oracolo
tf = unique(mod.tf$TF)
modulators = unique(mod.tf$MOD)
targets = setdiff(rownames(mexp),c(tf,mod))
#a=t(scale(t(mexp)))

# signif dei delta
out = NULL
for(mod in modulators) {
  ptm = proc.time()[3]
  tf1 = names(m[mod,m[mod,]>0])
  out=rbind(out,mindy2(mexp,mod=mod,tf=tf1,target = rownames(mexp),nboot = 1000,nbins=3,h=1,siglev="none",method="spearman")) # equiv mindy (nbins=3 h=1)
  print(paste0(mod," took ",proc.time()[3]-ptm," sec."))
}
#save(out,file="out-ATLAS.Rdata")
#save(out,file="out-ATLAS-5bins.Rdata")
#save(out,file="out-ATLAS-spearman.Rdata")
load("out-ATLAS.Rdata")
load("out-ATLAS-5bins.Rdata")
load("out-ATLAS-spearman.Rdata")
out$PVAL = p.adjust(out$PVAL,method="bonferroni")

out.sig=subset(out,abs(DELTA) > 0.4)
ntrg = aggregate(out.sig$TRG,by=list(out.sig$MOD,out.sig$TF),function(x) length(unique(x)))
ntrg[,4] = aggregate(out.sig$PVAL,by=list(out.sig$MOD,out.sig$TF),function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower=F))[,3]
ntrg[,4] = aggregate(out.sig$DELTA,by=list(out.sig$MOD,out.sig$TF),max)[,3]
ntrg[,5] = paste(ntrg[,1],ntrg[,2]) %in% oracle

require(ROCR)
tt = c("TP53")
for (tt in tf) {
  #ntrg.tf = ntrg[ntrg[,2]==tt,]
  ntrg.tf = ntrg
  pred = prediction(log(ntrg.tf[,4]),ntrg.tf[,5])
  pred = prediction(log(ntrg.tf[,3]),ntrg.tf[,5])
  perf = performance(pred,measure="prec",x.measure="rec")
  #perf = performance(pred,measure="tpr",x.measure="fpr")
  pdf(paste0("PR2-",tt,".pdf"))
  plot(perf,main=tt,lwd=5,col=2)
  grid()
  dev.off()
}

mexp = cut.outliers(mexp)

mytf="TP53"
for(mod in as.character(ntrg[ntrg[,2]==mytf,1])) {
  f=out.sig[out.sig$MOD==mod & out.sig$TF==mytf,]
  trg = unique(as.character(f$TRG))
  if (length(trg)<50 | !paste(mod,tt) %in% oracle) {
    cat(length(trg),"\n")
    next
  }
  fhigh=paste0(mytf,"-",mod,"-high")
  flow=paste0(mytf,"-",mod,"-low")
  con=file(description = paste0("outfigs2/LOW-",mytf,"-",mod,".tex"), open = "w")
  writeLines('\\documentclass[margin=5mm]{standalone}',con)
  writeLines('\\usepackage{graphics}',con)
  writeLines('\\usepackage{helvet}',con)
  writeLines('\\begin{document}',con)
  writeLines('\\begin{tabular}{cc}',con)
  writeLines(paste0('\\multicolumn{2}{c}{{{\\fontfamily{phv}\\selectfont {\\huge ',mod,' }}}}\\\\'),con)
  writeLines(paste0('\\includegraphics{',flow,'} &'),con)
  writeLines(paste0('\\includegraphics{',fhigh,'} \\\\'),con)
  writeLines('\\end{tabular}',con)
  writeLines('\\end{document}',con)
  close(con)
  
  pdf(paste0("outfigs2/",fhigh,".pdf"))
  plot.mod(mexp,mod,mytf,target = trg,
           nettarget = trg,fus="",high=T)
  dev.off()
  pdf(paste0("outfigs2/",flow,".pdf"))
  plot.mod(mexp,mod,mytf,target = trg,
           nettarget = trg,fus="",high=F)
  dev.off()
}


x=out.sig[out.sig$TF==mytf,]
out.sig = subset(out,DELTA>0.6)
f1=unique(paste(out.sig$MOD,out.sig$TF))

out.sig=out.sig[order(out.sig$PVAL),]
f=c()
for(i in 1:nrow(out.sig)) {
  l=paste(out.sig[i,1],out.sig[i,2])
  if (!l %in% f) {
    f = c(f,l)
  }
}
P=cumsum(f %in% oracle)/(1:length(f))
R=cumsum(f %in% oracle)/(length(oracle))
plot(R,P,type="l")
sum(f1 %in% oracle)/length(f1)
sum(intersect(f1,f) %in% oracle)
length(f1)
ntrg = aggregate(out.sig$TRG,by=list(out.sig$MOD,out.sig$TF),function(x) length(unique(x)))

boxplot(ntrg[ntrg[,4],5],ntrg[!ntrg[,4],5])


pr = c()
rc = c()
for (d in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7)) {
  f=paste(out$MOD,out$TF)
  
  [abs(out$DELTA)>0]
  sum(f %in% oracle)/length(oracle)
  sum(f %in% oracle)/length(f)
  length(oracle)/(length(tf)*length(modulators))
  pr=c(pr,P)
  rc=c(rc,R)
}
plot(rc,pr,type="l")

f=f[!duplicated(f)]
pred = prediction(length(f):1,f %in% oracle)
perf = performance(pred,measure="prec",x.measure="rec")
perf = performance(pred,measure="tpr",x.measure="fpr")
plot(perf)

#out$PVAL=p.adjust(out$PVAL,method = "fdr")
out.sig = subset(out,abs(DELTA)>0.5)
nrow(out.sig)
ntrg = aggregate(out.sig$TRG,by=list(out.sig$MOD,out.sig$TF),function(x) length(unique(x)))
#ntrg[,4] = aggregate(out.sig$PVAL,by=list(out.sig$MOD,out.sig$TF),mean)[,3]
ntrg[,4] = abs(aggregate(out.sig$DELTA,by=list(out.sig$MOD,out.sig$TF),max)[,3])

require(ROCR)
layout(matrix(1:3,nrow=1,ncol=3))
for (tt in tf) {
  x = grep(tt,f,value = T)
  pred = prediction(length(x):1,x %in% oracle)
  perf = performance(pred,measure="prec",x.measure="rec")
  plot(perf,ylim=c(0,1),xlim=c(0,1),main=tt)
}

mexps = cut.outliers(mexp)
mytf="TP53"
mymod="ATM" #modlist[[mytf]][3]
x=mindy2(mexp,mod=mymod,tf=names(modlist),target = rownames(mexp),nboot = 1000,nbins=3,h=1,siglev=0.01,method="pearson")

mytf="MYC"
tt = trg[[mytf]]
tt = as.character(out[out$MOD==mymod & out$TF==mytf,"TRG"])
plot.mod(mexps,mod=mymod,tf = mytf,target = tt,nettarget = tt,high=T)

perf = performance(pred,measure="tpr",x.measure="fpr")
plot(perf)
perf = performance(pred,measure="auc")

out.sig = subset(out,(DELTA>0.3 | DELTA< -0.3) & PVAL<0.05)
nrow(out.sig)
ntrg = aggregate(out.sig$TRG,by=list(out.sig$MOD,out.sig$TF),function(x) length(unique(x)))
ntrg[,4] = 0
ntrg[,5] = aggregate(out.sig$DELTA,by=list(out.sig$MOD,out.sig$TF),mean)[,3]


# significativita dei num di target per ogni mod-tf
outperm = NULL
for(b in 1:1000)  {
  ptm = proc.time()[3]
  outp=mindy2(mexp,mod=modulators[1],tf=tf,target = targets,nboot = 1000,nbins=3,h=1,perm=T,siglev=0.01,method="pearson")
  outp$B=b
  outperm = rbind(outperm,outp)
  print(paste0(b," boot took ",proc.time()[3]-ptm," sec."))
}
#save(outperm,file="outperm-ATLAS.Rdata")
load("outperm-ATLAS.Rdata")

outperm$PVAL = p.adjust(outperm$PVAL,method="fdr")

outperm.sig=subset(outperm,PVAL<0.01)
nrow(outperm.sig)
ntrg.perm = aggregate(outperm.sig$TRG,by=list(outperm.sig$MOD,outperm.sig$TF,outperm.sig$B),function(x) length(unique(x)))
plot(density(ntrg.perm$x))
ntrg[,4]=0
for(i in 1:nrow(ntrg)) {
  mod=as.character(ntrg[i,1])
  tfi=as.character(ntrg[i,2])
  over=sum(ntrg[i,3] < ntrg.perm[ntrg.perm$MOD==mod & ntrg.perm$TF==tfi,4])
  cat(i,mod,tfi,over,"\n")
  ntrg[i,4] = over
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



