From_freq_to_cell<-function(SNV_list,FREEC_list=NULL,Sample_names,Genotype_provided=FALSE,save_plot=TRUE,
                            contamination,ncores = 4, restrict.to.AB = FALSE,output_directory=NULL,
                            force.single.copy = FALSE){
  if(save_plot){
    if(is.null(output_directory)){
      dir.create(path = paste(Sample_names[1]), showWarnings = FALSE)
    }
    else{
      dir.create(path=output_directory,showWarnings = FALSE)
    }
  }
  if(Genotype_provided){
    Schrod_out<-Patient_schrodinger_cellularities(SNV_list = SNV_list,Genotype_provided =TRUE,
                                                  contamination = contamination, restrict.to.AB = restrict.to.AB,
                                                  force.single.copy = force.single.copy)
  }
  else{
    Schrod_out<-Patient_schrodinger_cellularities(SNV_list = SNV_list, FREEC_list = FREEC_list,
                                                  contamination = contamination, restrict.to.AB = restrict.to.AB,
                                                  force.single.copy = force.single.copy)
  }
  
  if(save_plot){
    plot_cell_from_Return_out(Schrod_out,Sample_names,output_directory)
  }
  
  return(Schrod_out)
}

Patient_schrodinger_cellularities<-function(SNV_list,FREEC_list=NULL,Genotype_provided=FALSE,
                                            contamination, restrict.to.AB = FALSE,
                                            force.single.copy = FALSE){
  result<-list()
  count<-0
  id<-1
  chr<-SNV_list[[1]][,"Chr"]
  chr_ante<-0
  if(!Genotype_provided){
    ChrCol <- which(colnames(FREEC_list[[1]]) == "Chromosome")
  }
  ### Compute all possible cellularities and join them
  for(i in 1:nrow(SNV_list[[1]])){ ##Exploring all possibilities for each mutation
    Cell<-list()
    test<-TRUE
    #print(SNV_list[[1]][i,'Chr'])
    for(k in 1:length(SNV_list)){
      if(test){ ## Do not look at position if it is invalid in another sample previously explored
        if(!is.null(SNV_list[[k]]$subclone.genotype)){
          if(is.na(SNV_list[[k]]$subclone.genotype)){
            subclone.geno<-NULL
          }
        }
        else{
          subclone.geno<-SNV_list[[k]][i,"subclone.genotype"]
        }
        if(Genotype_provided){
          Cell[[k]]<-cbind(CellularitiesFromFreq(Genotype= as.character(SNV_list[[k]][i,'Genotype']),
                                                 Alt = SNV_list[[k]][i,"Alt"],
                                                 Depth = SNV_list[[k]][i,"Depth"],
                                                 subclone.genotype = subclone.geno,
                                                 subclone.cell = SNV_list[[k]][i,"subclone.cell"],
                                                 chr = SNV_list[[k]][i,'Chr'],
                                                 position = SNV_list[[k]][i,'Start'],
                                                 contamination = contamination[k],
                                                 restrict.to.AB = restrict.to.AB,
                                                 force.single.copy = force.single.copy),
                           id)
        }
        else{
          if(chr[i]!=chr_ante){
            #subset on working chr
            #print(ChrCol)
            #print(FREEC_list[[1]][,1] == chr[i])
            CHR_FREEC<-lapply(FREEC_list,
                              function(z){
                                z[as.character(z[,ChrCol])==as.character(chr[i]),]
                                
                              } )
            chr_ante<-chr[i]
          }
          Cell[[k]]<-cbind(CellularitiesFromFreq(Freec_ratio = CHR_FREEC[[k]],
                                                 Alt = SNV_list[[k]][i,"Alt"],Depth = SNV_list[[k]][i,"Depth"],
                                                 subclone.genotype = subclone.geno,
                                                 subclone.cell = SNV_list[[k]][i,"subclone.cell"],
                                                 chr = SNV_list[[k]][i,'Chr'], position = SNV_list[[k]][i,'Start'],
                                                 contamination = contamination[k],
                                                 restrict.to.AB = restrict.to.AB,
                                                 force.single.copy = force.single.copy),
                           id)
        }
        if(sum(is.na(Cell[[k]]))>0){
          test<-F
          count<-count+1
        }
      }
    }
    if(test){ ## Checking that genotype is available
      L<-list()
      for(r in 1:length(Cell)){
        L[[r]]<-1:(nrow(Cell[[r]]))
      }
      U<-expand.grid(L) ##Table with all Cell row combinations
      for(k in 1:(ncol(U))){
        if(id==1){
          result[[k]]<-Cell[[k]][U[,k],]
        }
        else{
          result[[k]]<-rbind(result[[k]],Cell[[k]][U[,k],])
        }
      }
      id<-id+1
    }
  }
  
  if(count>0){
    warning(paste(count,'mutations excluded due to missing genotype or normalization issues'))
  }
  return(result)
}


FullEM<-function(Schrod, nclust, prior_center, prior_weight=NULL, 
                 contamination, epsilon=5*10**(-3),
                 optim = "default"
){
  if(length(prior_weight!=nclust)){
    prior_weight<-rep(1/nclust,times = nclust)
  }
  E_out<-EM.algo(Schrod = Schrod, nclust = nclust,
                 prior_center = prior_center, prior_weight = prior_weight, 
                 contamination = contamination, epsilon = epsilon,
                 optim = optim)
  if(is.list(E_out)){
    F_out<-filter_on_fik(Schrod = Schrod,fik = E_out$fik)
  }
  return(list(EM.output = E_out, filtered.data=F_out))
}


CellularitiesFromFreq<-function(chr, position,Alt,Depth,
                                Freec_ratio=NULL, Genotype=NULL,subclone.genotype=NULL,
                                subclone.cell=NULL,contamination, restrict.to.AB = FALSE,
                                force.single.copy = FALSE){
  ##For 1 mutation
  if(!is.null(Freec_ratio)){
    if(grepl(pattern = "chr",x = chr,ignore.case = T)){
      FChr<-sapply(X = 'chr',FUN = paste, Freec_ratio[,'Chromosome'],sep='')
    }
    else{
      FChr<-Freec_ratio[,'Chromosome']
    }
    Genotype<-as.character(tail(Freec_ratio[FChr==chr & Freec_ratio[,'Start']<position,'Genotype'],1))
  }
  if(length(Genotype)==0){
    message(paste("Genotype is not available at position",chr,position))
    return(NA)
  }
  else if(is.na(Genotype)){
    message(paste("Genotype is not available at position",chr,position))
    return(NA)
  }
  else if(Genotype==-1){
    message(paste("Genotype is not available at position",chr,position))
    return(NA)
  }
  else if(restrict.to.AB & (Genotype!='A' & Genotype!='AB')){
    message(paste("Position",chr,position,"is not in an A or AB site. Genotype:",Genotype))
    return(NA)
  }
  else if (is.null(subclone.genotype) | is.null(subclone.cell)){
    As<-strcount(x = Genotype, pattern = 'A',split = '')
    Ns<-nchar(Genotype)
    if(force.single.copy){
      ### Only one possibility per mutation
      result<-data.frame(Chr = chr,
                         Start =  position, 
                         Cellularity = as.numeric(
                           {Alt/Depth}*{Ns + contamination/(1-contamination)*2}
                         ),
                         Genotype = Genotype,
                         Alt = Alt,
                         Depth = Depth,
                         NC = 1,
                         NCh = Ns)
      
    }
    else{
      
      ### Vectorized version:
      result<-data.frame(Chr = chr,
                         Start = position,
                         Cellularity = as.numeric(
                           {Alt/Depth}*{Ns + contamination/(1-contamination)*2}/(1:As)
                         ),
                         Genotype = Genotype,
                         Alt = Alt,
                         Depth = Depth,
                         NC = 1:As,
                         NCh = Ns)
    }
  }
  else{ ## Two possibilities: belong to clone or subclone
    result<-data.frame()
    As<-strcount(x = Genotype, pattern = 'A',split = '')
    Ns<-nchar(Genotype)
    for(i in 1:As){
      Cellularity<-as.numeric(frequency/100*{Ns+contamination/(1-contamination)*2}/i)
      spare<-data.frame(chr,position,Cellularity, Genotype,Alt,Depth,i,Ns)
      colnames(spare)<-c('Chr','Start','Cellularity','Genotype',"Alt","Depth","NC","NCh")
      result<-rbind(result,spare)
    }
    A.sub<-strcount(x=subclone.genotype,pattern = "A",split = '')
    N.sub<-nchar(subclone.genotype)
    for(j in A.sub){
      ### Keep only possibilities that have a cellularity lower than the subclonal cellularity
      
      Cellularity<-as.numeric(frequency/100*{N.sub+contamination/(1-contamination)*2}/j)
      if(Cellularity<=subclone.cell){
        spare<-data.frame(chr,position,Cellularity, subclone.genotype,Alt,Depth,j,N.sub)
        colnames(spare)<-c('Chr','Start','Cellularity','Genotype',"Alt","Depth","NC","NCh")
        result<-rbind(result,spare)
      }
    }
  }
  
  result<-data.frame(result)
  colnames(result)<-c('Chr','Start','Cellularity','Genotype',"Alt","Depth","NC","NCh")
  return(result)
}

strcount <- function(x, pattern='', split=''){
  
  unlist(lapply(strsplit(x, split),function(z) na.omit(length(grep(pattern, z)))))
}

create_priors<-function(nclust,nsample,prior=NULL){
  result<-list()
  if(is.null(prior)){
    for(i in 1:nsample){
      result[[i]]<-c(runif(n = nclust-1,min = 0,max = 1),1)
    }
    return(result)
  }
  else if(length(prior[[1]])<nclust){## Need to complete the list
    if(sum(list_prod(prior)==1)>0){ ## there is an ancestral clone in the priors given
      for(i in 1:nsample){
        result[[i]]<-c(prior[[i]],runif(n = nclust-length(prior[[i]])))
      }
      return(result)
    }
    else{##need to add ancestral clone
      for(i in 1:nsample){
        result[[i]]<-c(prior[[i]],runif(n = nclust-1-length(prior[[i]])),1)
      }
      return(result)
    }
  }
  else if(length(prior[[1]])==nclust){
    return(prior)
  }
  else{## need to remove elements
    lp<-list_prod(prior)
    if(sum(lp>0.95**nsample)>0){ ## there is an ancestral clone in the priors given
      w<-which.max(lp>0.95**nsample)
      for(i in 1:nsample){
        result[[i]]<-c(sample(x = prior[[i]],size = nclust-1,replace = F),prior[[i]][w])   
      }
      return(result)
    }
    else{
      for(i in 1:nsample){
        result[[i]]<-c(sample(x = prior[[i]],size = nclust-1,replace = F),1)
      }
      return(result)
    }
  }
}

add.to.list<-function(...){
  c(as.list(...))
}

Compute.adj.fact<-function(Schrod){ ##Factor used to compute the probability of the binomial distribution
  n<-length(Schrod)
  adj.factor<-matrix(ncol = n,nrow=nrow(Schrod[[1]]))
  for(i in 1:n){
    adj.factor[,i]<-Schrod[[i]]$NC/Schrod[[i]]$NCh
  }
  return(adj.factor)
}


EM.algo<-function(Schrod, nclust=NULL,
                  prior_center=NULL,prior_weight=NULL,
                  contamination, epsilon=10**(-2),
                  optim = "default"
){
  if(is.null(prior_weight)){
    prior_weight<-rep(1/nclust,times = nclust)
    cur.weight<-rep(1/nclust,times = nclust)
  }
  else{
    cur.weight<-prior_weight
  }
  if(is.null(prior_center)){
    prior_center<-c(runif(n = (nclust-1)*length(Schrod),min = 0,max = 1),rep(1,times = length(Schrod)))
  }
  else{
    cur.center<-prior_center
  }
  prior_center<-unlist(cur.center)
  cur.val<-NULL
  eval<-1
  adj.factor<-Compute.adj.fact(Schrod = Schrod)
  if(grepl(pattern = optim,x = "compound",ignore.case = TRUE)){
    if(is.matrix(adj.factor) && ncol(adj.factor)>1){
      unicity_test<-TRUE
      for(i in 1:ncol(adj.factor)){
        if(length(unique(adj.factor[,i]))>1){
          unicity_test<-FALSE
        }
      }
      if(unicity_test){
        optim<-"exact"
      }
      else{
        optim<-"default"
      }
    }
    else{
      if(length(unique(adj.factor))==1){
        message("EM evaluation...")
        optim<-"exact"
      }
      else{
        message("default use...")
        optim<-"default"
      }
    }
  }
  if(optim!="DEoptim"){
    iters<-0
    while(eval>epsilon){
      if(optim == "exact"){
        iters<-iters+1 ### exact can be stuck with meta stable values
        ### Contradictory with convergence of EM...
      }
      tik<-e.step(Schrod = Schrod,centers = cur.center,weights = cur.weight,
                  adj.factor = adj.factor)
      m<-m.step(fik = tik,Schrod = Schrod,previous.weights = cur.weight,
                previous.centers =cur.center,
                adj.factor=adj.factor,optim = optim)
      
      if(!is.list(m)){
        test<-create_priors(nclust = 2,nsample = 2)
        eval_1<-max(abs(prior_center-unlist(test)))
        break      
      }
      else{
        n.weights<-unlist(m$weights)
        n.centers<-list()
        n.val<-m$val
        
        for(i in 1:length(cur.center)){
          n.centers[[i]]<-m$centers[((i-1)*length(cur.center[[1]])+1):((i)*length(cur.center[[1]]))]
        }
        eval<-max(abs(c(n.weights,unlist(n.centers))-c(cur.weight,unlist(cur.center))))
        cur.weight<-n.weights
        #prior_center<-c(prior_center,unlist(n.centers))
        cur.center<-n.centers
        ### Add fik*log(weights) if EM not direct optimization
        
      }
      PI<-matrix(nrow = nrow(tik),ncol= ncol(tik))
      for(i in 1:length(cur.weight)){
        PI[,i]<-tik[,i]*log(cur.weight[i])
      }
      PI[PI==0]<-0
      cur.val<-n.val - sum(PI)
    }
    fik<-e.step(Schrod = Schrod,
                centers = cur.center,
                weights = cur.weight,
                adj.factor = adj.factor)
    
    return(list(fik=fik,weights=cur.weight,centers=cur.center,val=cur.val))
  }
  else{
    ### DIRECT EVALUATION WITH DEoptim
    m<-m.step(fik = NULL,Schrod = Schrod,previous.weights = rep(1,times = length(prior_weight)),
              previous.centers =unlist(prior_center),
              adj.factor=adj.factor,optim = "DEoptim")
    fik<-e.step(Schrod = Schrod,
                centers = m$centers,
                adj.factor = adj.factor,
                rep(1,times = length(prior_weight))
    )
    cur.val<-sum(fik * eval.fik.m(Schrod= Schrod,
                                  centers = m$centers,
                                  weights = cur.weight,
                                  adj.factor = adj.factor,
                                  log = TRUE))
    return(list(fik=fik,weights=cur.weight,centers=m$centers,
                val=cur.val,initialpop = m$itialpop))
  }
}

eval.fik<-function(Schrod,centers,weights,keep.all.poss=TRUE,adj.factor,log = FALSE){
  al<-list()
  if(is.list(centers)){
    centers<-unlist(centers)
  }
  idx<-0
  if(log){
    al<-matrix(data = 0,nrow=nrow(Schrod[[1]]),ncol=length(weights))
  }
  else{
    al<-matrix(data = 1,nrow=nrow(Schrod[[1]]),ncol=length(weights))
  }
  for(i in 1:length(Schrod)){ ## i is a sample
    Alt<-Schrod[[i]]$Alt
    Depth<-Schrod[[i]]$Depth
    adj<-adj.factor[,i]
    
    
    for(k in 1:length(weights)){ ## k is a clone
      idx<-idx+1
      pro<-centers[idx]*adj
      test<-pro <=1 & pro >=0
      #pro_0<-pro
      #pro_0[pro>1 | pro<0]<-0
      if(log){
        al[test,k]<-al[test,k]+ifelse(test = Alt[test]==0,
                                      yes = Depth[test]*log(1 - pro[test]),
                                      no =  dbinom(x =Alt[test] ,size = Depth[test],prob = pro[test],log = TRUE)
        )
        al[!test,k]<-log(.Machine$double.xmin)
        
      }
      else{
        al[test,k]<-al[test,k]*dbinom(x =Alt[test] ,size = Depth[test],prob = pro[test],log = FALSE)
        al[!test,k]<-sqrt(.Machine$double.xmin)
      }
    }
  }
  al
}

grzero<-function(fik,adj.factor,Alt,Depth){
  #adj.factor has samples in cols
  #fik has clusters in cols
  if(is.matrix(Alt) && ncol(Alt)>1){
    centers<-numeric(length = ncol(fik)*ncol(Alt))
    index<-0
    for(s in 1:ncol(Alt)){
      for(k in 1:ncol(fik)){
        if(sum(fik[,k])==0){ ### The cluster has 0 probability
          fik[,k]<-.Machine$double.eps
        }
        index<-index+1
        centers[index]<-sum(fik[,k]*Alt[,s])/{adj.factor[1,s]*sum(fik[,k]*Depth[,s])} 
      }
    }
  }
  else{
    centers<-numeric(length = ncol(fik))
    for(k in 1:ncol(fik)){
      centers[k]<-sum(fik[,k]*Alt)/{adj.factor[1]*sum(fik[,k]*Depth)}
      
    }
  }
  centers
}

filter_on_fik<-function(Schrod,fik){
  tmp<-unique(Schrod[[1]]$id)
  keep<-numeric(length = length(tmp))
  
  for(i in 1:length(tmp)){
    u<-Schrod[[1]]$id==tmp[i]
    if(sum(u)>1){
      spare<-fik[u,]
      M<-max(spare)
      if(sum(spare==M)==1){
        l<-which(apply(X = spare,MARGIN = 1,FUN = function(z) sum(z== M)>0))
      }
      else{
        l<-which(apply(X = spare,MARGIN = 1,FUN = function(z) sum(z== M)>0))
        if(length(l)>1){
          l<-l[which.max(apply(X = spare[l,],MARGIN = 1,FUN = sum))]
        }
      }
      keep[i]<-which(u)[l]
    }
    else{
      keep[i]<-which(u)
    }
  }
  result<-Schrod
  for(l in 1:length(Schrod)){
    result[[l]]<-result[[l]][keep,]
  }
  return(result)
}

hard.clustering<-function(EM_out){
  clust<-apply(X = EM_out$fik,MARGIN = 1,FUN = function(z) {
    if(sum(z==max(z))>1){ ### Look for the multiple clones, and attribute with probability proportional to the weight
      if(max(z)>0){
        pos<-which(z==max(z))
        prob<-EM_out$weights[pos]/(sum(EM_out$weights[pos]))
        #sample(x = pos, size = 1, prob = prob))
        return(pos[which.max(prob)])
      }
      else{ ### all possibilities have 0 probability, so choose one randomly
        return(sample(1:length(z),size = 1))
      }
    }
    else{ # only one clone has maximal probability
      return(which.max(z))
    }
  })
  return(clust)
}


BIC_criterion<-function(EM_out_list,model.selection){
  ### Criterion should be minimized
  # Here we assimilate EM.output$val to -ln(L) where L is the likelihood of the model
  # BIC is written -2*ln(L)+k*ln(k)
  # Generalized BIC is written -2*ln(L)+q * k*ln(k)
  # AIC is written 2*k - 2*ln(L)
  
  if(is.numeric(model.selection)){
    ### Modified BIC to relax or add constraints on model selection
    # if q > 1 adding explicative variables should explain observed values better => control overfitting
    # if q = 1 BIC
    # if q < 1 adding explicative variables is less costly
    
    
    Bic<-numeric()
    if(length(EM_out_list)==0){
      return(0)
    }
    Mut_num<-nrow(EM_out_list[[1]]$EM.output$fik)
    for(i in 1:length(EM_out_list)){
      k<-length(unlist(EM_out_list[[i]]$EM.output$centers))
      Bic[i]<-2*EM_out_list[[i]]$EM.output$val+model.selection * k *log(Mut_num)
    }
    W<-which.min(Bic)
    L<-0
    return(Bic)
  }
  else if(model.selection == "BIC"){
    Bic<-numeric()
    if(length(EM_out_list)==0){
      return(0)
    }
    Mut_num<-nrow(EM_out_list[[1]]$EM.output$fik)
    
    for(i in 1:length(EM_out_list)){
      k<-length(unlist(EM_out_list[[i]]$EM.output$centers))
      Bic[i]<-2*EM_out_list[[i]]$EM.output$val+k*log(Mut_num)
    }
    return(Bic)
  }
  else if(model.selection == "AIC"){
    Aic<-numeric()
    if(length(EM_out_list)==0){
      return(0)
    }
    for(i in 1:length(EM_out_list)){
      k<-length(unlist(EM_out_list[[i]]$EM.output$centers))
      Aic[i]<-2*EM_out_list[[i]]$EM.output$val+2*k
    }
    return(Aic)
    
  }
  
}

m.step<-function(fik,Schrod,previous.weights,
                 previous.centers,contamination,adj.factor,
                 optim ="default"
){
  if(!is.null(fik)){
    weights<-apply(X = fik,MARGIN = 2,FUN = function(z) sum(z)/length(unique(Schrod[[1]]$id)))
  }
  else{
    weights<-rep(1/length(previous.weights),times = length(previous.weights))
  }
  # weights<-weights/sum(weights) # Overkill
  cur.cent<-list()
  
  if(optim == "default" | optim == "optimx" | optim == "exact"){
    Alt<-matrix(nrow = nrow(Schrod[[1]]),ncol = length(Schrod))
    Depth<-matrix(nrow = nrow(Schrod[[1]]),ncol = length(Schrod))
    
    for(i in 1:length(Schrod)){
      Alt[,i]<-Schrod[[i]]$Alt
      Depth[,i]<-Schrod[[i]]$Depth
    }
  }
  
  ### Function for maximization step
  fnx<-compiler::cmpfun(function(x) {
    r<--fik*eval.fik.m(Schrod = Schrod,centers = x,adj.factor = adj.factor,
                       weights = weights,
                       log = TRUE)
    
    r[fik==0]<-0
    sum(r,
        na.rm = TRUE
    )},
    options = list(optimize = 3)
  )
  
  ### Exact function with recomputation of fik
  efnx<-compiler::cmpfun(function(x) {
    fik<-e.step(Schrod = Schrod, centers = x, weights = weights,adj.factor = adj.factor)
    r<--fik*eval.fik.m(Schrod = Schrod,centers = x,adj.factor = adj.factor,
                       weights = weights,
                       log = TRUE)
    
    
    r[fik==0]<-0
    PI<-matrix(nrow = nrow(fik),ncol = ncol(fik))
    
    for(i in 1:length(weights)){
      PI[,i]<-fik[,i]*log(weights[i])
    }
    PI[fik==0]<-0
    sum(r-PI,
        na.rm = TRUE
    )},
    options = list(optimize = 3)
  )
  
  if(optim == "default"){
    spare<-tryCatch(optim(par = unlist(previous.centers),
                          fn = fnx ,
                          gr= function(x){grbase(fik = fik,adj.factor = adj.factor,centers = x,Alt = Alt,Depth=Depth)},
                          method = "L-BFGS-B",
                          lower = rep(.Machine$double.eps,times = length(unlist(previous.centers))),
                          upper=rep(1,length(unlist(previous.centers)))),
                    #### IF FAILS DUE TO INFINITE VALUE:
                    ####################################
                    error = function(e){
                      message("Gradient failed")
                      optim(par = unlist(previous.centers),
                            fn = fnx ,
                            method = "L-BFGS-B",
                            lower = rep(.Machine$double.eps,times = length(unlist(previous.centers))),
                            upper=rep(1,length(unlist(previous.centers)))
                      )
                    }
    )
    if(!is.list(spare)){
      return(NA)
    }
    return(list(weights=weights,centers=spare$par,val=spare$val))
  }
  else if(optim =="optimx"){
    
    spare<-optimx::optimx(par = unlist(previous.centers),
                          fn = fnx,
                          method = "L-BFGS-B",
                          lower = rep(.Machine$double.eps,times = length(unlist(previous.centers))),
                          upper=rep(1,length(unlist(previous.centers))))
    
    return(list(weights=weights,centers=spare[1:length(unlist(previous.centers))],val=spare$value))
  }
  else if(optim =="DEoptim"){
    spare<-suppressWarnings(DEoptim::DEoptim(fn = efnx,
                                             lower = rep(0,times = length(unlist(previous.centers))),
                                             upper=rep(1,length(unlist(previous.centers))),
                                             control = DEoptim::DEoptim.control(
                                               NP = min(10*length(unlist(previous.centers)),40),
                                               strategy= 3,
                                               itermax = 200,
                                               initialpop = NULL,
                                               CR = 0.9
                                             )
    )
    )
    
    return(list(weights = weights, centers = spare$optim$bestmem,val = spare$optim$bestval, 
                initialpop = spare$member$pop,itermax = 200)
    )
  }
  else if(optim == "exact"){
    new.centers<-grzero(fik,adj.factor,Alt,Depth)
    val<-fnx(new.centers)
    return(list(weights = weights,centers = new.centers, val = val))
  }
  return(list(weights=weights,centers=spare[1:length(unlist(previous.centers))],val=spare$value))
}

Cluster_plot_from_cell<-function(Cell,Sample_names,simulated,save_plot=TRUE,
                                 contamination, clone_priors,prior_weight,nclone_range,Initializations,preclustering=TRUE,
                                 epsilon=5*(10**(-3)),ncores = 2,output_directory=NULL,
                                 model.selection = "BIC",optim = "default", keep.all.models = FALSE){
  preclustering_success<-FALSE
  preclustering_FLASH<-FALSE
  ##### PRECLUSTERING
  ###################
  if(!is.null(preclustering)){
    if(preclustering=="kmedoid"){
      for(i in 1:length(Cell)){
        if(i==1){
          select<-Cell[[1]][,"Genotype"]=='A'| Cell[[1]][,"Genotype"]=='AB'
          Spare<-Cell[[1]][,'Cellularity']
        }
        else{
          select<- select & (Cell[[i]][,"Genotype"]=='A'| Cell[[i]][,"Genotype"]=='AB')
          Spare<-cbind(Spare,Cell[[i]][,'Cellularity'])
        }
      }
      if(length(Cell)==1){
        Spare<-Spare[select]
        if(length(Spare)==0){
          warning("No A and AB sites to do preclustering")
          p_clone<-clone_priors
          p_weight<-prior_weight
        }
      }
      else{
        Spare<-Spare[select,] ##restriction to A AB sites
      }
      if(is.null(dim(Spare))){
        warning("No A and AB sites to do preclustering")
        p_clone<-clone_priors
        p_weight<-prior_weight
      }
      else{
        if(nrow(Spare)<=max(nclone_range)){
          warning("Too few mutations to cluster. Will use priors / random initial conditions")
          p_clone<-clone_priors
          p_weight<-prior_weight
        }
        else{
          Spare[Spare>1]<-1
          kmeans<-fpc::pamk(Spare,krange = nclone_range,usepam = F)
          p_clone<-list()
          p_weight<-numeric()
          create_prior_weight<- nrow(Spare)>=50
          for(j in 1:(ncol(kmeans$pamobject$medoids))){
            p_clone[[j]]<-as.numeric(kmeans$pamobject$medoids[,j])
            if(create_prior_weight){
              p_weight[j]<-sum(kmeans$pamobject$clustering==j)/length(kmeans$pamobject$clustering)
            }
          }
          if(!create_prior_weight){
            p_weight<-prior_weight
          }
          preclustering_success<-T
        }
      }
    }
    else{ ### "FLASH" preclustering
      preclustering_FLASH<-TRUE
      result<-EM_clustering(Schrod = Cell,contamination = contamination,
                            epsilon = epsilon,
                            Initializations = Initializations,
                            nclone_range = nclone_range,
                            ncores = ncores,
                            model.selection = model.selection,
                            optim = optim, 
                            keep.all.models = keep.all.models,
                            FLASH = TRUE)
    }
  }
  else{
    p_clone<-clone_priors
    p_weight<-prior_weight
  }
  #################
  ### Clustering with EM
  #################
  if(!preclustering_FLASH){
    if(preclustering_success){
      result<-EM_clustering(Schrod = Cell,contamination = contamination,epsilon = epsilon,
                            prior_weight = p_weight,clone_priors = p_clone,Initializations = Initializations,
                            nclone_range = nclone_range,ncores = ncores,
                            model.selection = model.selection,optim = optim, 
                            keep.all.models = keep.all.models)
    }
    else{
      result<-EM_clustering(Schrod = Cell,contamination = contamination,epsilon = epsilon,
                            prior_weight = p_weight,clone_priors = p_clone,Initializations = Initializations,
                            nclone_range = nclone_range,ncores = ncores,
                            model.selection = model.selection,optim = optim,keep.all.models = keep.all.models)
    }
  }
  
  #### Plots possibilities
  
  if(save_plot && length(Cell)>1){
    plot_cell_from_Return_out(result$filtered.data,
                              Sample_names = Sample_names,
                              output_dir=output_directory)
  }
  return(result)
}

e.step<-function(Schrod,centers,weights,adj.factor){
  f<-eval.fik(Schrod = Schrod,centers = centers,weights = weights,
              adj.factor = adj.factor)
  for(k in 1:length(weights)){ ##k corresponds to a clone
    f[,k]<-f[,k]*weights[k]
  }
  ### Normalize fik by mutations
  f_0<-matrix(0,nrow = nrow(f), ncol = ncol(f))
  Id<-Schrod[[1]]$id
  for(m in unique(Id)){
    test<-Id==m
    tot<-sum(f[test,])
    if(tot == 0){
      f_0[test,]<-1/(sum(test)*ncol(f))  
    }
    else{
      f_0[test,]<-f[test,]/tot
    }
  }
  
  f_0
}

eval.fik.m<-function(Schrod,centers,weights,adj.factor,log = TRUE){
  spare<-eval.fik(Schrod = Schrod,
                  centers=centers,
                  weights = weights,
                  adj.factor= adj.factor,
                  log = log)
  test<-is.infinite(spare)
  if(sum(test)){
    spare[test]<-log(.Machine$double.xmin)
  }
  spare
}

EM_clustering<-function(Schrod,contamination,prior_weight=NULL, clone_priors=NULL, Initializations=1,
                        nclone_range=2:5, epsilon=0.01,ncores = 2,
                        model.selection = "BIC",optim = "default",keep.all.models = FALSE,
                        FLASH = FALSE){
  list_out_EM<-list()
  if(FLASH){
    tree<-Cellular_preclustering(Schrod)$tree
  }
  if(ncores >1){
    cl <- parallel::makeCluster( ncores )
    doParallel::registerDoParallel(cl)
    
    list_out_EM<-foreach::foreach(i=paste(rep(nclone_range,each = Initializations),c("",rep("_jit",times = Initializations-1))),
                                  ### jitter around priors if more than 1
                                  .export = c("parallelEM","FullEM","EM.algo","create_priors",
                                              "add.to.list","e.step","m.step","list_prod",
                                              "Compute.adj.fact","eval.fik","eval.fik.m",
                                              "filter_on_fik","Create_prior_cutTree","grbase")) %dopar% {
                                                
                                                if(FLASH){
                                                  if(grepl(pattern= "_",x= i)){
                                                    i<-as.numeric(unlist(strsplit(x = i,split = "_"))[1])
                                                    jitter <- TRUE
                                                  }
                                                  else{
                                                    i<-as.numeric(i)
                                                    jitter <- FALSE
                                                  }
                                                  priors<-Create_prior_cutTree(tree,Schrod,i,jitter)
                                                  return(parallelEM(Schrod = Schrod,nclust = i,epsilon = epsilon,
                                                                    contamination = contamination,prior_center = priors$centers,
                                                                    prior_weight = priors$weights,Initializations = 1,
                                                                    optim = optim,keep.all.models = keep.all.models                                                  )
                                                  )
                                                }
                                                else{
                                                  i<-as.numeric(unlist(strsplit(x = i,split = "_"))[1])
                                                  return(parallelEM(Schrod = Schrod,nclust = i,epsilon = epsilon,
                                                                    contamination = contamination,prior_center = clone_priors,
                                                                    prior_weight = prior_weight,Initializations = 1,
                                                                    optim = optim,keep.all.models = keep.all.models                                                  )
                                                  )
                                                }
                                              }
    #doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
  }
  else{
    index<-0
    for(i in 1:length(nclone_range)){
      for(init in 1:Initializations){
        if(FLASH){
          if(init == 1){
            priors<-Create_prior_cutTree(tree,Schrod,nclone_range[i],jitter = FALSE)
          }
          else{
            priors<-Create_prior_cutTree(tree,Schrod,nclone_range[i],jitter = TRUE)
            
          }
          index<-index+1
          list_out_EM[[index]]<-parallelEM(Schrod = Schrod,nclust = nclone_range[i],
                                           epsilon = epsilon,
                                           contamination = contamination,
                                           prior_center = priors$centers,
                                           prior_weight = priors$weights,
                                           Initializations = 1,
                                           optim = optim,
                                           keep.all.models = keep.all.models)
        }
        else{
          index<-index+1
          list_out_EM[[index]]<-parallelEM(Schrod = Schrod,nclust = nclone_range[i],epsilon = epsilon,
                                           contamination = contamination,prior_center = clone_priors,
                                           prior_weight = prior_weight,Initializations = Initializations,
                                           optim = optim,
                                           keep.all.models = keep.all.models)
        }
      }
    }
  }
  if(!keep.all.models){
    ### 
    # Criterion
    #
    result<-list_out_EM[[which.min(BIC_criterion(EM_out_list = list_out_EM, model.selection = model.selection))]]
    
    
    result$cluster<-hard.clustering(EM_out = result$EM.output)
    
    
    return(result)
  }
  else{
    Crit<-BIC_criterion(EM_out_list = list_out_EM, model.selection = model.selection)
    
    for(i in 1:length(list_out_EM)){
      list_out_EM[[i]]$cluster<-hard.clustering(EM_out = list_out_EM[[i]]$EM.output)
      list_out_EM[[i]]$Crit<-Crit[i]
    }
    
    return(list_out_EM)
  }
}

parallelEM<-function(Schrod,nclust,epsilon,contamination,
                     prior_center=NULL,prior_weight=NULL,
                     Initializations=1,
                     optim = "default",
                     keep.all.models = FALSE
){
  result<-list()
  for(i in 1:Initializations){
    result[[i]]<-FullEM(Schrod = Schrod,nclust = nclust,
                        prior_weight = prior_weight,
                        contamination = contamination,epsilon = epsilon,
                        prior_center = create_priors(nclust = nclust,
                                                     nsample = length(Schrod),
                                                     prior = prior_center),
                        optim = optim
    )
  } 
  if(keep.all.models){
    if(Initializations>1){
      return(result)
    }
    else{
      return(result[[1]])
    }
  }
  else{
    M<-result[[1]]$EM.output$val
    Mindex<-1
    if(length(result)>1){
      for(i in 2:length(result)){
        if(result[[i]]$EM.output$val<M){
          M<-result[[i]]$EM.output$val
          Mindex<-i
        }
      }
    }
    return(result[[Mindex]])
  }
}


Cellular_preclustering<-function(Schrod_cells){
  for(i in 1:length(Schrod_cells)){
    Schrod_cells[[i]]$Norm_Alt<-round(
      Schrod_cells[[i]]$Alt * Schrod_cells[[i]]$NCh / (2 * Schrod_cells[[i]]$NC)
    )
  }
  
  
  DistMat<-ProbDistMatrix(Schrod_cells)
  dissimMatrix<-as.dist(1-DistMat)
  tree<-hclust(d = dissimMatrix,method = "ward.D2")
  #cut_tree<-cutree(tree = tree,k = majority)
  result<-list(similarityMatrix = DistMat,
               distance = dissimMatrix,
               tree = tree
  )
  result 
}


ProbDistMatrix<-function(Schrod_cells){
  result<-matrix(0,nrow = nrow(Schrod_cells[[1]]),ncol = nrow(Schrod_cells[[1]]))
  for(l in 1:length(Schrod_cells)){
    ### z-score should be >0
    result<-result+zscore(Depth = Schrod_cells[[l]]$Depth,
                          Alt = Schrod_cells[[l]]$Norm_Alt
    )
  }
  2*pnorm(-result**(1/2))
  
}

zscore<-function(Depth,Alt){
  n<-length(Depth)
  result<-matrix(ncol = n,nrow = n)
  for(i in 1:n){
    p<-(Alt+Alt[i])/(Depth + Depth[i])
    w0<-which(p==0)
    w1<-which(p>1)
    p1<-Alt[i]/Depth[i]
    p2<-Alt/Depth
    result[,i]<-((p1-p2)**2)/(p*(1-p)*(1/Depth[i]+1/Depth))
    result[w0,i]<-0
    result[w1,i]<-+Inf
    result[is.na(result[,i]),i]<-+Inf
  }
  result
}


Create_prior_cutTree<-function(tree,Schrod_cells,NClus,jitter = FALSE){
  for(i in 1:length(Schrod_cells)){
    Schrod_cells[[i]]$Norm_Alt<-round(
      Schrod_cells[[i]]$Alt * Schrod_cells[[i]]$NCh / (2 * Schrod_cells[[i]]$NC)
    )
  }
  clustering<-cutree(tree,k = NClus)
  
  weights<-table(clustering)/length(clustering)
  if(jitter){
    weights<-jitter(weights)
    weights[weights<0]<-0
    weights<-weights/sum(weights)
  }
  centers<-list()
  for(i in 1:length(Schrod_cells)){
    if(jitter){
      centers[[i]]<-jitter(2 * unlist(sapply(as.numeric(names(weights)),
                                             function(cl){
                                               sum(Schrod_cells[[i]]$Norm_Alt[clustering==cl])/
                                                 sum(Schrod_cells[[i]]$Depth[clustering==cl])
                                             }
      )
      )
      )
    }
    else{
      centers[[i]]<-2 * unlist(sapply(as.numeric(names(weights)),
                                      function(cl){
                                        sum(Schrod_cells[[i]]$Norm_Alt[clustering==cl])/
                                          sum(Schrod_cells[[i]]$Depth[clustering==cl])
                                      }
      )
      )
    }
    centers[[i]][centers[[i]]>1]<-1
    centers[[i]][centers[[i]]<=0]<-.Machine$double.eps
    
  }
  list(weights = weights,centers = centers)
}


Tidy_output<-function(r, Genotype_provided,SNV_list){
  Work_input<-SNV_list[[1]]
  Work_output<-r$filtered.data[[1]]
  rownames(Work_input) <- NULL
  rownames(Work_output) <- NULL
  Work_input$Start <- trimws(Work_input$Start)
  Work_input$Depth <- trimws(Work_input$Depth)
  Work_output$Start <- trimws(Work_output$Start)
  Work_output$Depth <- trimws(Work_output$Depth)
  Work_input <- Work_input[, c("Chr", "Start", "Depth", "Alt")]
  Work_output <- Work_output[, c("Chr", "Start", "Depth", "Alt")]

  if(Genotype_provided){
    commonCols<-c("Chr","Start","Depth","Alt","Genotype")
  }else{
    commonCols<-c("Chr","Start","Depth","Alt")
  }
  Keep<-apply(Work_input[,c("Chr","Start","Depth")],MARGIN = 1, function(z){
    t<-apply(Work_output[,c("Chr","Start","Depth")],MARGIN = 1,function(y){
      return(all(y == z))
      
    })
    
    if(sum(t)==0){
      return(FALSE)
    }
    else{
      return(TRUE)
    }
  })
  Cell<-list()
  for(i in 1:length(SNV_list)){
    Cell[[i]]<-SNV_list[[i]][Keep,]
  }
  result<-r
  for(i in 1:length(r$filtered.data)){
    result$filtered.data[[i]]<-merge(r$filtered.data[[i]], Cell[[i]],by = commonCols)
  }
  
  ### Extract new order of mutations:
  Order<-apply(result$filtered.data[[1]][,c("Chr","Start","Depth")],MARGIN = 1, function(z){
    t<-apply(r$filtered.data[[1]][,c("Chr","Start","Depth")],MARGIN = 1,function(y){
      all(z == y)
    })
    
    return(which(t))
  })
  result$cluster<-result$cluster[Order]
  result$EM.output$fik<-result$EM.output$fik[Order,]
  
  return(result)
}