
mmi = function(mexp,tf,target,kordering,alltarget=TRUE,positiveOnly=F,ignore = 0.15,ncore=4,S=5,nboot=100,bfrac=0.8,sig=0.05, verbose=T) {
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
  library(foreach)
  library(doParallel)
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
  
  
  cl<-makeCluster(ncore)
  registerDoParallel(cl)  

  # ignoro il controllo di un range iniziale e finale con pochi sample
  # mi aspetto che il k risieda piu internamente
  inizio = round(ncol(kordering)*ignore)
  range= inizio:(ncol(kordering)-inizio)
  range = range[seq(1,length(range),S)] # prendo ogni S sample per ridurre il numero di k da controllare
  
  if(verbose) {
    print(paste0("Exp Matrix: ",paste0(dim(mexp),collapse="x")))
    print(paste0("N° targets: ",length(target)))
    print(paste0("N° modulators: ",nrow(kordering)))
    print(paste0("N° TF: ",length(tf)))
    print(paste0("Testing N° samples: ",length(range)))
  }
  
  # struttura dati ritornata
  # lista dei geni modulatori (cofattori o target). Per ogni modulatore una matrice 
  tfmod = list()
  
  pb = txtProgressBar(min=1,max=nrow(kordering)*nboot*length(range),style=3) 
  ri = 1 # contatore per la progress bar
  for (x in rownames(kordering)) {
    # per ogni modulatore candidato x-esimo
    # matrice 4-dimensionale TF x target x range  x nboot
    mi1k = array(0,dim=c(length(tf),length(target),length(range),nboot),dimnames=list(tf,target,range)) 
    if (!positiveOnly) mikn = mi1k 
    miall = mi1k
    for (bi in 1:nboot) {
      
      for(k in range) {
        ptm = proc.time()[3]
        ksample = sample(1:k,k*bfrac)
        mi1k[,,as.character(k),bi] = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,ksample]])
        mi1k[,,as.character(k),bi][mi1k[,,as.character(k),bi]<0] = 0

        if(!positiveOnly) {
          ksample = sample((ncol(kordering)-k):ncol(kordering),k*bfrac)
          mikn[,,as.character(k),bi] = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,ksample]])
          mikn[,,as.character(k),bi][mikn[,,as.character(k),bi]<0] = 0
        }
        
        # null
        ksample = sample(1:ncol(kordering),k*bfrac)
        miall[,,as.character(k),bi] = knnmi.cross(mexp[tf,kordering[x,ksample]],mexp[target,kordering[x,ksample]])
        miall[,,as.character(k),bi][miall[,,as.character(k),bi]<0] = 0
        if (verbose) print(proc.time()[3] - ptm)
        setTxtProgressBar(pb, ri)
        ri = ri + 1
      }

    }
    
    # calcolo il rapporto tra le medie delle MI tra i due sample boot
    # per ogni TF-target-k
    midelta1k=log2(apply(mi1k,c(1,2,3),sum)/apply(miall,c(1,2,3),sum))
    mipval1k = midelta1k*0+1
    if (!positiveOnly) {
      mideltakn=log2(apply(mikn,c(1,2,3),sum)/apply(miall,c(1,2,3),sum))
      mipvalkn = mideltakn*0+1
    }
    
    # effettuo il test tra i due sample boot
    # per ogni TF-target-k sia direct che inverse
    for(i in tf) {
      mipval1k[i,,] = foreach(t=target, .combine='rbind') %:% 
        foreach(j=as.character(range), .combine='c') %dopar% {
          wilcox.test(mi1k[i,t,j,],miall[i,t,j,],alternative="greater")$p.value
        }
      if (!positiveOnly) {
        mipvalkn[i,,] = foreach(t=target, .combine='rbind') %:% 
          foreach(j=as.character(range), .combine='c') %dopar% {
            wilcox.test(mikn[i,t,j,],miall[i,t,j,],alternative="greater")$p.value
          }
      }
    }

    # correzione dei pvalue
    mipval1k[1:length(mipval1k)] = p.adjust(mipval1k[1:length(mipval1k)],method="fdr")
    if (!positiveOnly) mipvalkn[1:length(mipvalkn)] = p.adjust(mipvalkn[1:length(mipvalkn)],method="fdr")
    
    
    # per ogni modulatore candidato x-esimo 
    # selezione del miglior k per ogni coppia TF->target
    # memorizzo solo le coppie TF-target con pval sig
    tmp = foreach(i=tf) %:% 
      foreach(j=target, .combine='rbind') %do% {
        b1k = which.min(mipval1k[i,j,])
        pv1k = mipval1k[i,j,b1k]
        retval=c()
        if (!positiveOnly) {
          bkn = which.min(mipvalkn[i,j,])
          pvkn = mipvalkn[i,j,bkn]
          if (pvkn <= sig) {
            retval = rbind(c(FDR=pvkn,DELTA=mideltakn[i,j,bkn],K=as.numeric(names(bkn)),DIRECTION= -1))
            rownames(retval)=j
          }
        }
        if (pv1k <=sig) {
          retval = rbind(retval,c(FDR=pv1k,DELTA=midelta1k[i,j,b1k],K=as.numeric(names(b1k)),DIRECTION= +1))
          rownames(retval)=rep(j,nrow(retval))
        }
        retval
      }
    names(tmp) = tf
    tfmod[[x]] = tmp
  }
  stopCluster(cl)
  close(pb)
  return(tfmod)
}

