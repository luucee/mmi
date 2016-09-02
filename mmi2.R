cohens_d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1
  md  <- abs(mean(x) - mean(y))
  csd <- lx * var(x) + ly * var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)     
  cd  <- md/csd        
  return(cd)
}

mmi = function(mexp,tflist,target,kordering,alltarget=TRUE,positiveOnly=F,ignore = 0.15,S=5,nboot=100,bfrac=0.8,sig=0.05, verbose=T,cl=NULL) {
  # mexp - matrice di espressione (geni sulle righe, samples sulle colonne)
  # tf - lista dei TF (deve essere un sottoinsieme di rownames(mexp))
  # target - lista dei target (deve essere un sottoinsieme di rownames(mexp)\tf)
  # kordering - matrice degli ordinamenti dei samples dei geni ritenuti modulatori
  #             (rownames(kordering) deve essere un sottinsieme di rownames(mexp))
  #             (colnames(kordering) deve essere identico a colnames(mexp)))
  #             Es.
  #             kordering fatto sulla metilazione dei target (mmet matrice di metilazione con matched samples)
  #             kordering = t(apply(mmet[target,],1,order,decreasing=F))
  #             in questo caso ha senso solo alltarget=FALSE e positiveOnly=T
  #
  #             kordering fatto sulla espressione di cofattori che modulano la regolazione TF->target
  #             kordering = t(apply(mexp[cofactors,],1,order,decreasing=F))
  #             in questo caso ha senso alltarget=TRUE e positiveOnly=F
  # alltarget - se uguale a TRUE (default) la mutua inf è calcolata su tutti i target, se FALSE è calcolata solo sul modulatore corrente 
  #             (utile nel caso di kordering fatto sulla metilazione dei target)
  # positiveOnly - se uguale a TRUE (default) effettua solo il test (1..k) vs all, altrimenti considera anche l'inverso (k..n) vs all
  # ignore - percentuale dei sample esterni da ignorare quando si calcola la MI
  # ncore - numero di core per sfruttare il calcolo parallelo
  # S - numero di samples da raggruppare per la ricerca di k (T=1 il test è fatto per ogni sample)
  # nboot - numero di bootstrap per generare le distribuzioni
  # bfrac - frazione dei sample per il bootstrap
  # sig - livello di significativita per filtrare i risultati
  
  require(parmigene)
  # precondition controls
  if (!is.matrix(mexp)) 
    stop("mexp must be a matrix")
  if ( length(intersect(tflist,rownames(mexp))) != length(tflist) ) 
    stop("tf must be included in rownames(mexp)")
  if ( length(intersect(target,rownames(mexp))) != length(target) ) 
    stop("target must be included in rownames(mexp)")
  if ( length(intersect(rownames(kordering),rownames(mexp))) != nrow(kordering) ) 
    stop("rownames(kordering) must be included in rownames(mexp)")
  if ( length(intersect(colnames(kordering),colnames(mexp))) != ncol(mexp) ) 
    stop("colnames(kordering) must be equal to colnames(mexp)")
  
  

  # ignoro il controllo di un range iniziale e finale con pochi sample
  # mi aspetto che il k risieda piu internamente
  #inizio = round(ncol(kordering)*ignore)
  #range= inizio:(ncol(kordering)-inizio)
  #range = range[seq(1,length(range),S)] # prendo ogni S sample per ridurre il numero di k da controllare
  sbin = ncol(kordering) %/% S
  range = seq(0,ncol(kordering),sbin)
  range = range[-1]
  if ((sbin*S)<=ncol(kordering)) {
    range = range[-length(range)]
  }
  
  if(verbose) {
    print(paste0("Exp Matrix: ",paste0(dim(mexp),collapse="x")))
    print(paste0("N° targets: ",length(target)))
    print(paste0("N° modulators: ",nrow(kordering)))
    print(paste0("N° TF: ",length(tflist)))
    print(paste0("N° intervals: ",length(range)))
    print(paste0("Bin size: ",sbin))
    print(paste0("N° boots: ",nboot))
  }
  
  # rimozione mod-tf dipendenti
  mi.mod = knnmi.cross(mexp[rownames(kordering),],mexp[tflist,])
  count = mi.mod*0
  for (bi in 1:nboot) {
    mi = knnmi.cross(mexp[rownames(kordering),],mexp[tflist,sample(1:ncol(mexp))])
    count = count + mi>mi.mod
  }
  pval.mod = count/nboot
  pval.mod[1:length(pval.mod)]=p.adjust(pval.mod[1:length(pval.mod)],method = "fdr")
  mi.mod[pval.mod > 0.01]=0
  
  # struttura dati ritornata
  # lista dei geni modulatori (cofattori o target). Per ogni modulatore una matrice 
  tfmod = list()
  
  if(!verbose) pb = txtProgressBar(min=1,max=nrow(kordering),style=3) 
  ri = 1 # contatore per la progress bar
  for (x in rownames(kordering)) {
    ptm = proc.time()[3]
    
    # considero solo i tf non dipendenti da x
    tf = names(mi.mod[x,tflist]==0)
    
    # per ogni modulatore candidato x-esimo
    # matrice 3-dimensionale delle MI TF x target x range
    mi1k = array(0,dim=c(length(tf),length(target),length(range)),dimnames=list(tf,target,range)) 
    mikn = mi1k 
    for (k in range) {
      ksample = 1:k
      mi1k[,,as.character(k)] = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,ksample]])
      if(!positiveOnly) {
        ksample = (ncol(kordering)-k):ncol(kordering)
        mikn[,,as.character(k)] = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,ksample]])
      }
    }
    
    miall.boot = array(0,dim=c(length(tf),length(target),length(range),nboot),dimnames=list(tf,target,range)) 
    mi1k.perm = miall.boot
    mikn.perm = miall.boot
    tmp = foreach (k = range) %:%
      foreach(bi = 1:nboot) %dopar% {
        require(parmigene)
        # all
        ksample = sample(1:ncol(kordering),k)
        tmp.miall = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,ksample]])

        ksample = 1:k
        tmp.mi1k.perm = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,sample(ksample)]])
        
        ksample = (ncol(kordering)-k):ncol(kordering)
        tmp.mikn.perm = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,sample(ksample)]])
        list(BootAll = tmp.miall,Perm1k = tmp.mi1k.perm,Permkn=tmp.mikn.perm)
        
      }
    for (bi in 1:nboot) {
      for (k in 1:length(range)) {
        miall.boot[,,as.character(range[k]),bi] = tmp[[k]][[bi]]$BootAll
        mi1k.perm[,,as.character(range[k]),bi] = tmp[[k]][[bi]]$Perm1k
        mikn.perm[,,as.character(range[k]),bi] = tmp[[k]][[bi]]$Permkn
      }
    }
    miall = apply(miall.boot,c(1,2,3),median)
    midelta1k = mi1k - miall    
    mideltakn = mikn - miall

    mipval1k = apply(midelta1k[1:length(midelta1k)] < (mi1k.perm - miall[1:length(miall)]),c(1,2,3),sum)/nboot
    mipvalkn = apply(mideltakn[1:length(mideltakn)] < (mikn.perm - miall[1:length(miall)]),c(1,2,3),sum)/nboot

    deltamindy = mi1k[,,as.character(range[1])] - mikn[,,as.character(range[1])]
    deltamindy.perm = mi1k.perm[,,as.character(range[1]),] - mikn.perm[,,as.character(range[1]),]

    pvalmindy = apply(deltamindy[1:length(deltamindy)] < deltamindy.perm,c(1,2),sum)/nboot
    
    if(verbose) {
      print(paste0(x," took ",proc.time()[3]-ptm," sec. "))
    } else {
      setTxtProgressBar(pb, ri)
      ri = ri + 1
    }
    
    tfmod[[x]] = list(DELTA1k=midelta1k,PVAL1k=mipval1k,DELTAkn=mideltakn,PVALkn=mipvalkn,
                      DELTAmindy=deltamindy,PVALmindy=pvalmindy)
  }
  return(tfmod)
}



summarization = function(mmiout) {
  # per ogni modulatore candidato x-esimo 
  # selezione del miglior k per ogni coppia TF->target
  # memorizzo solo le coppie TF-target con pval sig
  retval1 = c()
  retval2 = c()
  for (x in names(mmiout)) {
    cat(x,"\n")
    b1k = apply(mmiout[[x]]$PVAL1k,c(1,2),which.min)
    bkn = apply(mmiout[[x]]$PVALkn,c(1,2),which.min)
    for (i in rownames(b1k)) {
      for (j in colnames(b1k)) {
        pval1k = mmiout[[x]]$PVAL1k[i,j,b1k[i,j]]
        pvalkn = mmiout[[x]]$PVALkn[i,j,bkn[i,j]]
        retval1 = rbind(retval1,c(PVALkn=pvalkn,DELTAkn=mmiout[[x]]$DELTAkn[i,j,bkn[i,j]],Kkn=bkn[i,j],
                                  PVAL1k=pval1k,DELTA1k=mmiout[[x]]$DELTA1k[i,j,b1k[i,j]],K1k=b1k[i,j],
                                  PVALmindy = mmiout[[x]]$PVALmindy[i,j],DELTAmindy = mmiout[[x]]$DELTAmindy[i,j]))
        retval2 = rbind(retval2,c(MOD=x,TF=i,TRG=j))
      }
    }
  }
  return(data.frame(retval1,retval2,stringsAsFactors = F))
}


svm.mod = function(mexp,tf,target,kordering,alltarget=TRUE,positiveOnly=F,ignore = 0.15,S=5,nboot=100,bfrac=0.8,sig=0.05, verbose=T,cl=NULL) {
  # mexp - matrice di espressione (geni sulle righe, samples sulle colonne)
  # tf - lista dei TF (deve essere un sottoinsieme di rownames(mexp))
  # target - lista dei target (deve essere un sottoinsieme di rownames(mexp)\tf)
  # kordering - matrice degli ordinamenti dei samples dei geni ritenuti modulatori
  #             (rownames(kordering) deve essere un sottinsieme di rownames(mexp))
  #             (colnames(kordering) deve essere identico a colnames(mexp)))
  #             Es.
  #             kordering fatto sulla metilazione dei target (mmet matrice di metilazione con matched samples)
  #             kordering = t(apply(mmet[target,],1,order,decreasing=F))
  #             in questo caso ha senso solo alltarget=FALSE e positiveOnly=T
  #
  #             kordering fatto sulla espressione di cofattori che modulano la regolazione TF->target
  #             kordering = t(apply(mexp[cofactors,],1,order,decreasing=F))
  #             in questo caso ha senso alltarget=TRUE e positiveOnly=F
  # alltarget - se uguale a TRUE (default) la mutua inf è calcolata su tutti i target, se FALSE è calcolata solo sul modulatore corrente 
  #             (utile nel caso di kordering fatto sulla metilazione dei target)
  # positiveOnly - se uguale a TRUE (default) effettua solo il test (1..k) vs all, altrimenti considera anche l'inverso (k..n) vs all
  # ignore - percentuale dei sample esterni da ignorare quando si calcola la MI
  # ncore - numero di core per sfruttare il calcolo parallelo
  # S - numero di samples da raggruppare per la ricerca di k (T=1 il test è fatto per ogni sample)
  # nboot - numero di bootstrap per generare le distribuzioni
  # bfrac - frazione dei sample per il bootstrap
  # sig - livello di significativita per filtrare i risultati
  require(kernlab)
  
  for (x in rownames(kordering)) {
    # per ogni modulatore candidato x-esimo
    # generazione del training sample
    kordering[x,]
  tpos = NULL
  tneg = NULL
  for(i in 1:100) {
    trgi = sample(1:nrow(mexp),1)
    ptrg = sort(sample(mexp[trgi,],length(A)))
    ptrg[order(mexp[tfi,A])] = ptrg
    mvec = sort(mmet[trgi,])
    evec = sort(mexp[trgi,])
    tpos=rbind(tpos,c(ptrg,
                      sample(evec[1:length(B)]),
                      sample(mvec[1:length(A)]),
                      sample(mvec[(length(A)+1):length(mvec)])))
    
    
    #trgi = sample(1:nrow(mexp),1)
    ntrg=c(sample(mexp[trgi,A]),sample(mmet[trgi,]))
    tneg=rbind(tneg,c(sample(mexp[trgi,A]),
                      sample(mexp[trgi,],length(B)),
                      sample(mmet[trgi,])))
  }
  rownames(tpos)=NULL
  rownames(tneg)=NULL
  }
  
}



