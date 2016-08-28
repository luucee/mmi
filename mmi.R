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

mmi = function(mexp,tf,target,kordering,alltarget=TRUE,positiveOnly=F,ignore = 0.15,S=5,nboot=100,bfrac=0.8,sig=0.05, verbose=T,cl=NULL) {
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
  if ( length(intersect(tf,rownames(mexp))) != length(tf) ) 
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
    print(paste0("N° TF: ",length(tf)))
    print(paste0("N° intervals: ",length(range)))
    print(paste0("Bin size: ",sbin))
    print(paste0("N° boots: ",nboot))
  }
  
  # struttura dati ritornata
  # lista dei geni modulatori (cofattori o target). Per ogni modulatore una matrice 
  tfmod = list()
  
  if(!verbose) pb = txtProgressBar(min=1,max=nrow(kordering),style=3) 
  ri = 1 # contatore per la progress bar
  for (x in rownames(kordering)) {
    # per ogni modulatore candidato x-esimo
    # matrice 4-dimensionale TF x target x range  x nboot
    mi1k = array(0,dim=c(length(tf),length(target),length(range),nboot),dimnames=list(tf,target,range)) 
    if (!positiveOnly) {
      mikn = mi1k 
      mikn.null = mi1k
    }
    miall = mi1k
    mi1k.null = mi1k
    
    ptm = proc.time()[3]
    tmp = foreach (k = range) %:%
      foreach(bi = 1:nboot) %dopar% {
        require(parmigene)
        retval = list()
        ksample = sample(1:k,sbin*bfrac)
        tmp.mi1k = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,ksample]])
        tmp.mi1k.nl = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,sample(ksample)]])
        #tmp.mi1k[tmp.mi1k<0] = 0
        
        tmp.mikn=NULL
        tmp.mikn.nl=NULL
        if(!positiveOnly) {
          ksample = sample((ncol(kordering)-k):ncol(kordering),sbin*bfrac)
          tmp.mikn = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,ksample]])
          tmp.mikn.nl = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,sample(ksample)]])
          #tmp.mikn[tmp.mikn<0] = 0
        }
        
        # all
        ksample = sample(1:ncol(kordering),sbin*bfrac)
        tmp.miall = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,ksample]])
        #tmp.miall[tmp.miall<0] = 0
        
        retval = list(MI1k=tmp.mi1k,MIkn=tmp.mikn,MIall=tmp.miall,MI1knull=tmp.mi1k.nl,MIknnull=tmp.mikn.nl)
        retval
      }
    
    for (bi in 1:nboot) {
      for (k in 1:length(range)) {
        mi1k[,,as.character(range[k]),bi] = tmp[[k]][[bi]]$MI1k
        miall[,,as.character(range[k]),bi] = tmp[[k]][[bi]]$MIall
        mi1k.null[,,as.character(range[k]),bi] = tmp[[k]][[bi]]$MI1knull
        if(!positiveOnly) {
          mikn[,,as.character(range[k]),bi] = tmp[[k]][[bi]]$MIkn
          mikn.null[,,as.character(range[k]),bi] = tmp[[k]][[bi]]$MIknnull
        }
      }
    }
    

    #azzero con la soglia
    misoglia.1k=apply(mi1k.null,c(1,2,3),max)
    if(!positiveOnly) misoglia.kn=apply(mikn.null,c(1,2,3),max)
    
    # calcolo il cohen's d effect size tra i due bootstrap
    # per ogni TF-target-k
    #midelta1k = (apply(mi1k,c(1,2,3),mean)-apply(miall,c(1,2,3),mean))/sqrt((apply(mi1k,c(1,2,3),var)+apply(miall,c(1,2,3),var))/2)
    media.all = apply(miall,c(1,2,3),median)
    media.1k = apply(mi1k,c(1,2,3),median)
    media.kn = apply(mikn,c(1,2,3),median)
    midelta1k = media.1k/media.all
    midelta1k[media.1k <= misoglia.1k] = 0
    mipval1k = midelta1k*0+1
    mipvalkn=NULL
    mideltakn=NULL
    if (!positiveOnly) {
      #mideltakn = (apply(mikn,c(1,2,3),mean)-apply(miall,c(1,2,3),mean))/sqrt((apply(mikn,c(1,2,3),var)+apply(miall,c(1,2,3),var))/2)
      mideltakn = media.kn/media.all
      mideltakn[media.kn <= misoglia.kn] = 0
      mipvalkn = mideltakn*0+1
    }
    
    # effettuo il test tra i due sample boot
    # per ogni TF-target al sample k-esimo sia direct (positiveOnly==T) che inverse
    for(i in tf) {
      mipval1k[i,,] = foreach(t=target, .combine='rbind') %:% 
        foreach(j=as.character(range), .combine='c') %dopar% {
          if (i==t) {
            1
          } else {
            #t.test(mi1k[i,t,j,],miall[i,t,j,],alternative="greater")$p.value
            wilcox.test(mi1k[i,t,j,],miall[i,t,j,],alternative="greater")$p.value
          }
        }
      if (!positiveOnly) {
        mipvalkn[i,,] = foreach(t=target, .combine='rbind') %:% 
          foreach(j=as.character(range), .combine='c') %dopar% {
            if (i==t) {
              1
            } else {
              #t.test(mikn[i,t,j,],miall[i,t,j,],alternative="greater")$p.value
              wilcox.test(mikn[i,t,j,],miall[i,t,j,],alternative="greater")$p.value
            }
          }
      }
    }

    # correzione dei pvalue
    #mipval1k[1:length(mipval1k)] = p.adjust(mipval1k[1:length(mipval1k)],method="fdr")
    #if (!positiveOnly) mipvalkn[1:length(mipvalkn)] = p.adjust(mipvalkn[1:length(mipvalkn)],method="fdr")
    
    if(verbose) {
      print(paste0(x," took ",proc.time()[3]-ptm," sec. "))
    } else {
      setTxtProgressBar(pb, ri)
      ri = ri + 1
    }
    
    tfmod[[x]] = list(TARGET=target,TF=tf,DELTA1k=midelta1k,PVAL1k=mipval1k,DELTAkn=mideltakn,PVALkn=mipvalkn,
                      MIALL=media.all,MI1k=media.1k,MIkn=media.kn)
  }
  return(tfmod)
}

delta.clustering = function(tfmod,cutt=0.4) {
  # per ogni modulatore candidato x-esimo 
  # e per ogni TF i-esimo effettuo il clustering dei target
  # rispetto al log2 del rapporto
  retval = list()
  for (x in names(tfmod)) {
    tmp = foreach(i=tfmod[[x]]$TF) %dopar% {
      hc1k=hclust(tfmod[[x]]$DELTA1k[i,,]*(1-tfmod[[x]]$PVAL1k[i,,]),method = "ward.D2")
      ct1k=cutree(hc1k,h = cutt)
      ctkn=NULL
      if(!is.null(tfmod[[x]]$DELTAkn)) {
        hckn=hclust(tfmod[[x]]$DELTAkn[i,,]*(1-tfmod[[x]]$PVALkn[i,,]),method = "ward.D2")
        ctkn=cutree(hckn,h = cutt)
      }
      list(CT1k=ct1k,CTkn=CTkn)
    }
    names(tmp) = tf
    retval[[x]] = tmp
  }
  return(retval)
}

summarization = function(mmiout,siglev=0.01) {
  # per ogni modulatore candidato x-esimo 
  # selezione del miglior k per ogni coppia TF->target
  # memorizzo solo le coppie TF-target con pval sig
  retval = c()
  for (x in names(mmiout)) {
    b1k = apply(mmiout[[x]]$DELTA1k,c(1,2),which.max)
    bkn = apply(mmiout[[x]]$DELTAkn,c(1,2),which.max)
    for (i in rownames(b1k)) {
      for (j in colnames(b1k)) {
        pval = mmiout[[x]]$PVAL1k[i,j,b1k[i,j]]
        if(pval<siglev) {
          retval = rbind(retval,data.frame(MOD=x,TF=i,TRG=j,PVAL=pval,DELTA=mmiout[[x]]$DELTA1k[i,j,b1k[i,j]],MODdir= +1,
                                           MI=mmiout[[x]]$MI1k[i,j,b1k[i,j]],MIall=mmiout[[x]]$MIall[i,j,b1k[i,j]]))
        }
        pval = mmiout[[x]]$PVALkn[i,j,bkn[i,j]]
        if(pval<siglev) {
          retval = rbind(retval,data.frame(MOD=x,TF=i,TRG=j,PVAL=pval,DELTA=mmiout[[x]]$DELTA1k[i,j,b1k[i,j]],MODdir= -1,
                                           MI=mmiout[[x]]$MIkn[i,j,bkn[i,j]],MIall=mmiout[[x]]$MIall[i,j,bkn[i,j]]))
        }
      }
    }
  }
  return(retval)
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



