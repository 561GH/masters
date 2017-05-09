scmc=function(tausteps=seq(from=0,to=10,length=100),priorsamples,constraint=0,input=suse,par=par,type=type){#odefn(suse,par,type)
  #timestep=10;prior_fun=odefn;post_fun=m_state_svec;input=suse;par=par;type=type

  theta     = priorsamples
  w         = priorsamples[,1]*0+1#The value of w   numeric(timestep)
  W         = w/dim(priorsamples)[1] #The value of W
  step=0
  for(tau in tausteps[-1]){
    step=step+1
    w = pnorm(tau*theta)/pnorm(taustep[step+1]*theta)#updata wj
    W = w*W#update Wj
    # W[,(ts+1)] = (W[,(ts+1)]-mean(W[,(ts+1)]))/sd(W[,(ts+1)]) #Normalize Wj
    ESS= sum(W^2)
    if(ESS<(length(W)/2)){
      theta = sample(ori.post,prob=W)#resample
      W     = 1/dim(priorsamples)[1] #The value of W
    }else{
      theta  = sample(ori.post[,1],prob=W[,(ts+1)]/sum(W[,(ts+1)]))#need to change mh step
    }
  }#end of time step
  after.post=matrix(theta[,timestep+1],byrow=F,ncol=ncol(post_fun),nrow=nrow(post_fun))
  return(after.post)
  b }
