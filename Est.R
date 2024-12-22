###############################################################################
#  estimation procedures
# ###############################################################################
#data: including the observations
# X, V: linear and nonlinear covariates in the latency components
# Z, W: linear and nonlinear covariates in the incidence components
# L, R: left and right endpoints of intervals
# status: status for the observations
# vknots,V.df, V.order: number of knots, degree of freedom for the nonlinear function phi(v)
# Vgrid: evaluation points for the phi(v)
# wknots, W.df, W.order:



BDS_est<- function(data, order = 2, knots= NULL, nknots = 10,grids=grids,
                   vknots = NULL,V.df, V.order, Vgrid= Vgrid,
                   wknots = NULL, W.df, W.order, Wgrid= Wgrid,
                   a_tau=1,  b_tau=1, r=1, s0=1,sig0= 10, 
                   coef_range,sig_theta =10,coef_range_alpha,coef_range_gcoef, 
                   niter=niter,burnin=niter/4, alpha0,beta0)
{
  
  library(splines)
  library(icenReg)
  library(dlm)
  library(bayestestR)
  
  
  Z = cbind(data$Z.1,data$Z.2); W = data$W; 
  X = cbind(data$X.1,data$X.2); V = data$V;
  L = data$L; R = data$R;
  status = data$status;
  L=matrix(L,ncol=1)
  R=matrix(R,ncol=1)
  R2=ifelse(is.finite(R),R,0)
  xcov=as.matrix(X)
  zcov= as.matrix(Z)
  wcov = as.matrix(W)
  vcov = as.matrix(V)
  p=ncol(xcov)
  q= ncol(zcov)
  Status=matrix(status,ncol=1)
  n=nrow(L)
  
  
  ## generate cumulative basis functions 
  if (is.null(knots)) {knots<-seq(min(c(L,R2)),max(c(L,R2)),length=nknots)} # default for equally-based according to the wang lianming
  k=length(knots)-2+order  # quantile  the first value bigger than second?
  if (is.null(grids)) {grids<-seq(min(c(L,R2)),max(c(L,R2)),length=100)} #a sequence of points where baseline survival function is to be estimated
  kgrids=length(grids)
  
  # basis function for f(v)
  if (is.null(vknots)) {
    vknots = quantile(V,seq(0,1,length.out=V.df))
    vknots = vknots[-c(1,length(vknots))] # inner knots 
  }
  # basis function for g(W)
  if (is.null(wknots)) {
    wknots = quantile(W,seq(0,1,length.out= W.df))
    #wknots =c(wknots[1]+0.01, wknots[2:(length(wknots)-1)],wknots[length(wknots)]-0.01)
    wknots = wknots[-c(1,length(wknots))]
  }
  
  
  bisL = Ispline(L,order,knots)  # knots * n
  bisR = Ispline(R2,order,knots) # knots * n 
  bgs = Ispline(grids,order,knots)
  bsV = bs(V, knots = vknots, degree = V.order)  # n*kv
  bsVg = bs(Vgrid,knots = vknots, degree = V.order)
  bsW = bs(W,knots= wknots, degree = W.order)
  bsWg =  bs(Wgrid,knots= wknots, degree = W.order)
  kf = ncol(bsV);kfe = nrow(bsVg)
  kg = ncol(bsW); kgg = nrow(bsWg)
  
  
  
  
  ## initial value
  beta=matrix(rep(0,p),ncol=1)    # fit particial linear for the failure subject to obtain the intial value of fcoef?
  icsp =  ic_sp(Surv(L,R,type="interval2")~X.1+X.2+bs(V, knots = vknots,degree = V.order), data=data[is.finite(data$R),])
  fprior= icsp$coefficients[c(-1,-2)]; # intial value of coef
  fcoef = fprior
  rho = c(beta, fcoef)           
  gamcoef=matrix(rgamma(k, 1, 1),ncol=k) # coef for cumulative hazard function
  tau=rgamma(1,a_tau,rate=b_tau)  # hyperparameter for gamma
  
  u = ifelse(status==3, 0, 1)
  alpha = rep(0,q)
  glm_fit= glm(u~Z.1+Z.2+ bs(W,knots= wknots, degree = W.order)-1, family = binomial(),data = data)
  gprior= glm_fit$coefficients[c(-1,-2)];
  gcoef = rep(0, kg)
  gcoef = gprior
  theta = c(alpha, gcoef)
  
  
  LambdatL=t(gamcoef%*%bisL) # nt x 1 updata Lambda based on the gamcoef from the last iteration 
  LambdatR=t(gamcoef%*%bisR) 
  
  
  
  #  pary = array(rep(0,niter*n),dim=c(niter,n))  # indicator of uncured subject with 1 for uncured
  parbeta=array(rep(0,niter*p),dim=c(niter,p))  # X^T*beta+ phi(V)
  parfcoef = array(rep(0,niter*kf),dim=c(niter,kf))
  parf = array(rep(0,niter*kfe),dim=c(niter,kfe))
  pargam=array(rep(0,niter*k),dim=c(niter,k))  
  parLamb0=array(rep(0,niter*kgrids),dim=c(niter,kgrids))
  partau=array(rep(0,niter),dim=c(niter,1))  # hyperparameter of gamma 
  paralpha = array(rep(0,niter*q), dim= c(niter, q)) 
  pargcoef = array(rep(0,niter*kg), dim= c(niter, kg))
  parg = array(rep(0,niter*kgg), dim= c(niter, kgg))
  
  ## iteration
  iter=1
  while (iter<niter+1)
  { 
    # sample v, vv, w and ww
    v=array(rep(0,n),dim=c(n,1)); w=v
    vv=array(rep(0,n*k),dim=c(n,k)); ww=vv
    
    for (i in 1:n){
      if (status[i]==1){    # for left censored data
        templam1=LambdatR[i]*exp(xcov[i,]%*%beta+ bsV[i,]%*%fcoef)  
        v[i]=positivepoissonrnd(templam1)
        vv[i,]=rmultinom(1,v[i],gamcoef*t(bisR[,i]))
      }else if (status[i]==2){   # for right censored data
        templam1=(LambdatR[i]-LambdatL[i])*exp(xcov[i,]%*%beta+bsV[i,]%*%fcoef)
        w[i]=positivepoissonrnd(templam1)
        ww[i,]=rmultinom(1,w[i],gamcoef*t(bisR[,i]-bisL[,i]))
      }
    }
    
    # sample rho, rho = (beta, phicoef)^T
    te1=v*as.numeric(status==1)+w*as.numeric(status==2)
    te2=(LambdatR*as.numeric(status==1)+LambdatR*as.numeric(status==2)+LambdatL*as.numeric(status==3)*u)

    
    # for beta N(0,sig)
    xx = as.matrix(cbind(xcov, bsV))
    for (j in 1:p){
      rho[j]<-arms(rho[j],beta_fun,ind_fun,1,j=j,beta=rho,xx=xx,te1=te1,te2=te2,mu=c(0,0),sig=sig0,coef_range=coef_range)
    }
    
    # for fcoef N(f',sig)
    for (j in (p+1):(p+kf)){
      rho[j]<-arms(rho[j],beta_fun,ind_fun,1,j=j,beta=rho,xx=xx,te1=te1,te2=te2,mu =fprior, sig=sig0*10,coef_range=coef_range)
    }
    
    # update beta and phi(V)
    beta = as.matrix(rho[1:p],ncol=1)
    fcoef = rho[(p+1):(p+kf)]
    fest =  bsVg%*%fcoef 
    
    parbeta[iter,] = beta
    parfcoef[iter,] = fcoef
    parf[iter,] = fest
    
    
    # sample gamcoef
    for (l in 1:k){
      tempa=1+sum(vv[,l]*as.numeric(status==1)*(bisR[l,]>0)+ww[,l]*as.numeric(status==2)*((bisR[l,]-bisL[l,])>0))
      tempb=tau+sum((bisR[l,]*as.numeric(status==1)+bisR[l,]*as.numeric(status==2)+bisL[l,]*as.numeric(status==3)*u)*exp(xcov%*%beta+bsV%*%fcoef)) 
      gamcoef[l]=rgamma(1,tempa,rate=tempb)
    }
    
    # updata  whole LambdaL and LambdaR used for v and w
    LambdatL=t(gamcoef%*%bisL) 
    LambdatR=t(gamcoef%*%bisR)
    Lambdatg=t(gamcoef%*%bgs)
    pargam[iter,]= gamcoef
    parLamb0[iter,]= Lambdatg
    
    #sample tau
    tau=rgamma(1,a_tau+k, rate=b_tau+sum(gamcoef))
    
    
    zz=cbind(zcov, bsW)
    ## #for alpha
    for (i in 1:q) {
      theta[i] = arms(theta[i],theta_fun,indFunc=ind_fun_theta,1,i=i,theta=theta,u=u,mu1= c(0,0),zz=zz,sig_theta=sig_theta,
                      coef_range_theta=coef_range_alpha)
    }
    alpha = theta[1:q]
    paralpha[iter,] = alpha
    
    
    # for g(W)
    for (i in (q+1):(q+kg)) {
      theta[i] = arms(theta[i],theta_fun,indFunc=ind_fun_theta,1,i=i,theta=theta,u=u,mu1= gprior,zz=zz,sig_theta=sig_theta,
                      coef_range_theta=coef_range_gcoef)
    }
    gcoef = theta[(q+1):(q+kg)]
    pargcoef[iter,] = gcoef
    parg[iter,] = bsWg%*%gcoef
    

    p_i = exp(zcov%*%alpha+bsW%*%gcoef)/(1+exp(zcov%*%alpha+bsW%*%gcoef))
    SL = exp(-LambdatL*exp(xcov%*%beta+bsV%*%fcoef))  # when beta is univariate
    pest = p_i*SL/(1-p_i+p_i*SL)
    u= ifelse(status==3, rbinom(n,1, pest), 1)

    
    iter=iter+1
    if (iter%%100==0){cat("iteration = ", iter, "\n")}
  } # end iteration
  
  
  # evaluate the estimate 
  thin =10; burnin = niter/4
 # print(thin)
  wbeta=as.matrix(parbeta[seq((burnin+thin),niter,by=thin),],ncol=p) # thinned beta samples
  wfcoef = as.matrix(parfcoef[seq((burnin+thin),niter,by=thin),],ncol=kv)
  wf = as.matrix(parf[seq((burnin+thin),niter,by=thin),],ncol=kve)
  wgam=as.matrix(pargam[seq((burnin+thin),niter,by=thin),],ncol=k)
  wparLamb0=as.matrix(parLamb0[seq((burnin+thin),niter,by=thin),],ncol=kgrids)
  
  
  walpha = as.matrix(paralpha[seq((burnin+thin),niter,by=thin),],ncol=q) 
  wgcoef = as.matrix(pargcoef[seq((burnin+thin),niter,by=thin),],ncol=kg)
  wg = as.matrix(parg[seq((burnin+thin),niter,by=thin),],ncol=kgg)

  
  # assessement of convergence
  # plot traceplot in each replication 
  mcmc_beta = mcmc(wbeta)
  mcmc_alpha = mcmc(walpha)
  

  # pdf("traceplot of alpha1 in heavy sce1.pdf",height=7)
  # traceplot(mcmc_alpha[,1])
  # dev.off()
  # 
  # pdf("traceplot of alpha2 in heavy sce1.pdf",height=7)
  # traceplot(mcmc_alpha[,2])
  # dev.off()
  # 
  # pdf("traceplot of beta1 in heavy sce1.pdf",height=7)
  # traceplot(mcmc_beta[,1])
  # dev.off()
  # 
  # pdf("traceplot of beta2 in heavy sce1.pdf",height=7)
  # traceplot(mcmc_beta[,2])
  # dev.off()
  
  ess.beta = effectiveSize(mcmc_beta)
  ess.alpha = effectiveSize(mcmc_alpha)
  print(ess.beta)
  print(ess.alpha)
  

  
  ebeta = apply(wbeta,2,mean) ; se.beta =  apply(wbeta,2,sd)
  # Using HPD method
  CI.beta = HPDinterval(mcmc_beta, prob = 0.95)
  
  efcoef = apply(wfcoef,2,mean)
  ef = apply(wf,2,mean); #phi= 
  egam = apply(wgam,2,mean)
  eparLamb0 = apply(wparLamb0,2,mean)  # Lambda
  
  # esigma = apply(wsigma,2,mean)
  ealpha = apply(walpha,2,mean);se.alpha = apply(walpha,2,sd)
  CI.alpha =  HPDinterval(mcmc_alpha, prob = 0.95)
  
  cp.alpha = is_in_ci(alpha0, CI.alpha)
  egcoef = apply(wgcoef,2,mean)
  eg = apply(wg,2,mean)
  
  
  # save the estimates
  est = list(ebeta = ebeta,se.beta = se.beta,cp.beta = cp.beta, Vgrid = Vgrid,efcoef =efcoef,ef = ef, egam = egam, grids = grids, eparLamb= eparLamb0,
             ealpha = ealpha,se.alpha= se.alpha,cp.alpha = cp.alpha,egcoef = egcoef,eg =eg )
  est
  
}

