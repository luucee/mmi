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
  range = seq(1,ncol(kordering),sbin)
  range = range[-1]
  range = range[-length(range)]
  
  if(verbose) {
    print(paste0("Exp Matrix: ",paste0(dim(mexp),collapse="x")))
    print(paste0("N° targets: ",length(target)))
    print(paste0("N° modulators: ",nrow(kordering)))
    print(paste0("N° TF: ",length(tf)))
    print(paste0("Testing N° samples: ",length(range)))
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
    if (!positiveOnly) mikn = mi1k 
    miall = mi1k
    minull = mi1k
    ptm = proc.time()[3]
    tmp = foreach (k = range) %:%
      foreach(bi = 1:nboot) %dopar% {
        require(parmigene)
        retval = list()
        ksample = sample(1:k,sbin*bfrac)
        tmp.mi1k = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,ksample]])
        #tmp.mi1k[tmp.mi1k<0] = 0
        
        tmp.mikn=NULL
        if(!positiveOnly) {
          ksample = sample((ncol(kordering)-k):ncol(kordering),sbin*bfrac)
          tmp.mikn = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,ksample]])
          #tmp.mikn[tmp.mikn<0] = 0
        }
        
        # all
        ksample = sample(1:ncol(kordering),sbin*bfrac)
        tmp.miall = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,ksample]])
        #tmp.miall[tmp.miall<0] = 0

        # null
        ksample = sample(1:ncol(kordering),sbin*bfrac)
        tmp.minull = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,sample(ksample)]])
        
        retval = list(MI1k=tmp.mi1k,MIkn=tmp.mikn,MIall=tmp.miall,MInull=tmp.minull)
        retval
      }
    
    for (bi in 1:nboot) {
      for (k in 1:length(range)) {
        mi1k[,,as.character(range[k]),bi] = tmp[[k]][[bi]]$MI1k
        miall[,,as.character(range[k]),bi] = tmp[[k]][[bi]]$MIall
        minull[,,as.character(range[k]),bi] = tmp[[k]][[bi]]$MInull
        if(!positiveOnly) {
          mikn[,,as.character(range[k]),bi] = tmp[[k]][[bi]]$MIkn
        }
      }
    }
    

    #azzero con la soglia
    misoglia=apply(minull,c(1,2,3),max)
    
    # calcolo il cohen's d effect size tra i due bootstrap
    # per ogni TF-target-k
    #midelta1k = (apply(mi1k,c(1,2,3),mean)-apply(miall,c(1,2,3),mean))/sqrt((apply(mi1k,c(1,2,3),var)+apply(miall,c(1,2,3),var))/2)
    midelta1k = (abs(apply(mi1k,c(1,2,3),sum))/abs(apply(miall,c(1,2,3),sum)))
    midelta1k[apply(mikn,c(1,2,3),sum)<misoglia] = 0
    mipval1k = midelta1k*0+1
    mipvalkn=NULL
    mideltakn=NULL
    if (!positiveOnly) {
      #mideltakn = (apply(mikn,c(1,2,3),mean)-apply(miall,c(1,2,3),mean))/sqrt((apply(mikn,c(1,2,3),var)+apply(miall,c(1,2,3),var))/2)
      mideltakn=(abs(apply(mikn,c(1,2,3),sum))/abs(apply(miall,c(1,2,3),sum)))
      mideltakn[apply(mikn,c(1,2,3),sum)<misoglia] = 0
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
    mipval1k[1:length(mipval1k)] = p.adjust(mipval1k[1:length(mipval1k)],method="fdr")
    if (!positiveOnly) mipvalkn[1:length(mipvalkn)] = p.adjust(mipvalkn[1:length(mipvalkn)],method="fdr")
    
    if(verbose) {
      print(paste0(x," took ",proc.time()[3]-ptm," sec. "))
    } else {
      setTxtProgressBar(pb, ri)
      ri = ri + 1
    }
    
    tfmod[[x]] = list(TARGET=target,TF=tf,DELTA1k=midelta1k,PVAL1k=mipval1k,DELTAkn=mideltakn,PVALkn=mipvalkn)
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

summarization = function(tfmod,deltacl) {
  # per ogni modulatore candidato x-esimo 
  # selezione del miglior k per ogni coppia TF->target
  # memorizzo solo le coppie TF-target con pval sig
  retvalSum = list()
  for (x in names(tfmod)) {
    tmp = foreach(i=tfmod[[x]]$TF) %:%
      foreach(j=tfmod[[x]]$TARGET, .combine='rbind') %dopar% {
        b1k = which.min(tfmod[[x]]$PVAL1k[i,j,])
        pv1k = tfmod[[x]]$PVAL1k[i,j,b1k]
        retval=c()
        if (!positiveOnly) {
          bkn = which.min(tfmod[[x]]$PVALkn[i,j,])
          pvkn = tfmod[[x]]$PVALkn[i,j,bkn]
          if (pvkn <= sig) {
            retval = rbind(c(FDR=pvkn,CL=deltacl[[x]]$CTkn,
                             DELTA=tfmod[[x]]$DELTAkn[i,j,bkn],K=as.numeric(names(bkn)),DIRECTION= -1))
            rownames(retval)=j
          }
        }
        if (pv1k <=sig) {
          retval = rbind(retval,c(FDR=pv1k,CL=deltacl[[x]]$CT1k,
                                  DELTA=tfmod[[x]]$DELTA1k[i,j,b1k],K=as.numeric(names(b1k)),DIRECTION= +1))
          rownames(retval)=rep(j,nrow(retval))
        }
        retval
      }
    retvalSum[[x]] = tmp
  }
  return(retvalSum)
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



