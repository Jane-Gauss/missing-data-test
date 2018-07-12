    rm(list=ls())
    library(kyotil)
    library(robustrank)
    library(pracma)
    library(MASS)
    #########################################################################
    #N-cycle index
    #alternative-you can choose three ways "two.sided","less","greater" to calculate the p-value
    #distr-you can choose different kinds of distribution "normal","student","logistic","lognormal"
    #rho-correlation coefficient
    #m-Xpaired and Ypaired's number
    #n.x-Xextra's number
    #n.y-Yextra's number
    
     #generate random numbers
    sim.partially.matched1=function(m, n.x, n.y, distr=c("normal","logistic","student","mixnormal","gamma","lognormal","beta","uniform","hybrid1","hybrid2","doublexp"), params, seed){
    set.seed(seed)    
    N=m+n.x+n.y    
    distr=match.arg(distr)
    if (distr %in% c("normal","student","logistic","lognormal","doublexp")) {
        loc.2=params["loc.2"]
        scale.2=params["scale.2"]
        rho=params["rho"]; if (is.na(rho)) rho=0        
        if (distr=="normal") {
            dat.all=mvtnorm::rmvnorm(n = N, mean=c(0,loc.2),  sigma=matrix(c(1,scale.2*rho,scale.2*rho,scale.2^2),nrow=2))        
        } else if (distr=="lognormal") {
            dat.all=exp(mvtnorm::rmvnorm(n = N, mean=c(0,0),  sigma=matrix(c(1,scale.2*rho,scale.2*rho,scale.2^2),nrow=2)) )
            dat.all[,2]=dat.all[,2]+loc.2
        } else if (distr=="student") {        
            dat.all=mvtnorm::rmvt   (n = N, delta=c(0,loc.2), sigma=matrix(c(1,scale.2*rho,scale.2*rho,scale.2^2),nrow=2), df=3)       
        } else if (distr=="logistic") {
            dat.all=rbilogistic(N, 0, loc.2, scale.1=1, scale.2=scale.2, rho) 
        } else if (distr=="doublexp") { # double exponential also known as Laplace
            dat.all=rbidoublexp(N, 0, loc.2, scale.1=1, scale.2=scale.2, rho) 
        }
    } else if (distr=="mixnormal") {        
        p.1=params["p.1"]
        p.2=params["p.2"]
        sd.n=params["sd.n"]; if (is.na(sd.n)) sd.n=1
        sd1=params["sd1"];   if (is.na(sd1)) sd1=0.5
        sd2=params["sd2"];   if (is.na(sd2)) sd2=0.5
        dat.all=rnorm(N, sd=sd.n) + cbind(rmixnorm (N, mix.p=p.1, mu1=0, mu2=2, sd1=sd1, sd2=sd2), rmixnorm (N, mix.p=p.2, mu1=0, mu2=2, sd1=sd1, sd2=sd2))            
    } else if (distr=="gamma") {        
        loc.2=params["loc.2"]
        shape.1=params["shape.1"]
        shape.2=params["shape.2"]
        rate.1=params["rate.1"]
        rate.2=params["rate.2"]
        rho=params["rho"]
        dat.all=rbigamma(n=N, shape.1, shape.2, rate.1, rate.2, rho=rho) 
        dat.all[,2]=dat.all[,2]+loc.2
    } else if (distr=="beta") {        
        shape1.1=params["shape1.1"]
        shape2.1=params["shape2.1"]
        shape1.2=params["shape1.2"]
        shape2.2=params["shape2.2"]
        dat.all=cbind(rbeta(N,shape1.1,shape2.1),rbeta(N,shape1.2,shape2.2))
    } else if (distr=="uniform") {        
        loc.2=params["loc.2"]
        dat.all=cbind(runif(N),runif(N))
        dat.all[,2]=dat.all[,2]+loc.2
    } else if (distr=="hybrid1") {
        # X gamma Y logistic
        loc.2=params["loc.2"]
        toggle=params["toggle"]
        if (toggle==0) {
            dat.all=cbind(rgamma(N,shape=2,rate=0.5), rlogis(N)+loc.2)
        } else if (toggle==1) {
            dat.all=cbind(rlogis(N)+loc.2, rgamma(N,shape=2,rate=0.5))
        }
    } else if (distr=="hybrid2") {
        # X lognormal Y logistic
        loc.2=params["loc.2"]
        toggle=params["toggle"]
        if (toggle==0) {
            dat.all=cbind(exp(rnorm(N)), rlogis(N)+loc.2)
        } else if (toggle==1) {
            dat.all=cbind(rlogis(N)+loc.2, exp(rnorm(N)))
        }
        
    }    
    index.1=if(n.x==0) NULL else m+1:n.x
    index.2=if(n.y==0) NULL else m+n.x+1:n.y
    list(X=dat.all[1:m,1],            Y=dat.all[1:m,2], 
         Xprime=dat.all[index.1,1],   Ymissing=dat.all[index.1,2], 
         Xmissing=dat.all[index.2,1], Yprime=dat.all[index.2,2])         
}
    
   
    ranktests=function(N,alternative,distr,rho,m,n.x,n.y,loc.2){
    	C=matrix(0,nrow=N,ncol=4);#records of the results    
    dat=sim.partially.matched1(m*N,n.x*N,n.y*N,distr=distr,params=c("loc.2"=loc.2,"rho"=rho,"scale.2"=1),seed=1)  
    
    for(q in 1:N){
    s.1=(q-1)*m+1
    s.2=q*m
    Xpaired=dat$X[s.1:s.2]
    Ypaired=dat$Y[s.1:s.2]
    if(n.x==0){
    	Xextra=NULL
    	}else{
    s.11=(q-1)*n.x+1
    s.12=q*n.x
    Xextra=dat$Xprime[s.11:s.12]
    }
    s.21=(q-1)*n.y+1
    s.22=q*n.y
    Yextra=dat$Yprime[s.21:s.22]
    useC=FALSE
    correct=!useC# not implemented in .Call implemenation; partially implemented in R implementation 
    is.switched=FALSE #whether data conversion occurs
  
   # consolidate extra data if there are NA in paired data
    Yextra=c(Yextra,Ypaired[is.na(Xpaired)])
    Yextra=Yextra[!is.na(Yextra)]
    Xextra=c(Xextra,Xpaired[is.na(Ypaired)])
    Xextra=Xextra[!is.na(Xextra)]
    # remove NA from paired data
    notna=!is.na(Xpaired) & !is.na(Ypaired)
    Xpaired=Xpaired[notna]
    Ypaired=Ypaired[notna]
    # switch X and Y if Xextra is not null but Yextra is
    is.switched=FALSE
    if (length(Yextra)==0 & length(Xextra)>0) {
        cat("switch two samples \n")
        tmp=Xpaired; Xpaired=Ypaired; Ypaired=tmp
        Yextra=Xextra
        Xextra=NULL
        is.switched=TRUE
    }    
    
    if (length(Xextra)<=0 & length(Yextra)<=0) {
        cat("call wilcox.test since there are no unpaired samples\n")
        return(wilcox.test(Xpaired, Ypaired, paired=TRUE))
    } 
    
    m=length(Xpaired)
    n=length(Yextra)
    lx=length(Xextra) # so that it does not get mixed up with number 1
    N=m+n
    alpha=n/N
    beta=N/(m+N)
    
    if(!correct) .corr=0 else .corr=switch(alternative, "two.sided" = 2, "greater" = -1, "less" = 1)     
    if(is.switched & abs(.corr)==1) .corr=-.corr
   
    
   YYprime=c(Ypaired,Yextra)
    
    delta=Ypaired-Xpaired    
    abs.delta=abs(delta)
    sgn.delta=sign(delta)
    pos.delta=ifelse(delta>0,1,0)
    
    # empirical distribution function estimated from observed
    Fx=ecdf(Xpaired)
    Fy=ecdf(Ypaired)
    Fyprime=ecdf(Yextra)
    Fyyprime=ecdf(YYprime)
    Fplus=ecdf(abs.delta)
    Fxyyprime=ecdf(c(Xpaired,Ypaired,Yextra))
    cov.G.F=cov(Fy(Xpaired), Fx(Ypaired)) 
    cov.G.F.1=cov(Fyyprime(Xpaired), Fx(Ypaired)) 
    
    # compute h.1 for W.plus
    # Vectorize the following does not help speed b/c the loop is still implemented in R
    h.1.f=function(i) (sum(delta+delta[i]>0) - sum(delta[i]>0))/(length(delta)-1)
    h.1=numeric(m)
    for (i in 1:m) h.1[i]=h.1.f(i)
#    #using outer is actually slower probably b/c memory operations
#    delta.pair.plus=outer(delta, delta, FUN="+")
#    diag(delta.pair.plus)=NA # we don't want to count delta_i + delta_j
#    h.1.0=apply(delta.pair.plus>0, 1, mean, na.rm=TRUE)
#    all(h.1==h.1.0)
    
    # Wilcoxon signed rank stat and var
    r.1 <- rank(abs.delta)    
    W.plus = mean(pos.delta*r.1)
    if(correct) W.plus=W.plus-switch(alternative, "two.sided" = sign(W.plus-(m+1)/4) * 0.5/m, "greater" = ifelse(is.switched,-1,1)*0.5/m, "less" = -ifelse(is.switched,-1,1)*0.5/m)
    W.plus.2=mean(sgn.delta*r.1)
    var.W.plus = m * var(h.1) # var(W.plus.2) is 4 x var(W.plus)
    var.W.plus.0 = (m+1)*(2*m+1)/(24*m)
    
    # Wilcoxon-Mann-Whitney stat and var
    r.2 <- rank(c(Xpaired, Yextra))
    W.mw=sum(r.2[(m+1):N])/(N+1)
    var.W.mw=m * (alpha^2*var(Fyprime(Xpaired)) + alpha*(1-alpha)*var(Fx(Yextra)))
    var.W.mw.0=m*n/(12*(N+1))        
    if(correct) W.mw=W.mw-switch(alternative, "two.sided" = sign(W.mw-n/2) * 0.5/(N+1), "greater" = ifelse(is.switched,-1,1)*0.5/(N+1), "less" = -ifelse(is.switched,-1,1)*0.5/(N+1))

    # covariance between W.plus and W.mw    
    cov. = -m*alpha * cov(h.1, Fy(Xpaired))
    rho=cov./sqrt(var.W.plus*var.W.mw)
    rho.0 = - 6*sqrt(alpha) * cov(Fplus(abs.delta) * sgn.delta, Fx(Xpaired))
    cov.0 =  rho.0  * sqrt(var.W.plus.0*var.W.mw.0) # good
    # we may also put in the population mean in the estimate, it reduces bias but increases var greatly, not good
    rho.0a= - 6*sqrt(alpha) * mean(Fplus(abs.delta) * sgn.delta * Fx(Xpaired)) 
    cov.0a = rho.0a * sqrt(var.W.plus.0*var.W.mw.0) 
    # other choices tried that do not make big differences
    # cov.0d = (2*rho-rho.0) * sqrt(var.W.plus.0*var.W.mw.0) # 
    #rho.0b = - 6*sqrt(alpha) * cov(Fplus(abs.delta) * sgn.delta, Fyprime(Xpaired))
    #cov.0  = rho    * sqrt(var.W.plus.0*var.W.mw.0) # 
    
    # MW-MW1 and MW-MW2
    U.p=(sum(rank(c(Xpaired, Ypaired))[m+1:m]) - m*(m+1)/2)/(m*m)
    if(correct) U.p=U.p-switch(alternative, "two.sided" = sign(U.p-0) * 0.5/(m*m), "greater" = ifelse(is.switched,-1,1)*0.5/(m*m), "less" = -ifelse(is.switched,-1,1)*0.5/(m*m))
    U.mw=(W.mw*(N+1)-n*(n+1)/2)/(m*n)
    if(correct) U.mw=U.mw-switch(alternative, "two.sided" = sign(U.mw-0) * 0.5/(m*n), "greater" = ifelse(is.switched,-1,1)*0.5/(m*n), "less" = -ifelse(is.switched,-1,1)*0.5/(m*n))
    var.U.p.0=1/6-2*cov.G.F
    # more precise est of var.U.p.0
    #var.U.p.0.high=m*((U.p-0.5)/.Call("pair_wmw_test", Xpaired, Ypaired, .corr, as.integer(2), as.integer(1), as.integer(1)))^2
    #var.U.p.0=var.U.p.0.high
    var.U.mw.0=1/(12*alpha)
    cov.U.mw.p.0=1/12-cov.G.F
    var.U.p =  var(Fy(Xpaired)) + var(Fx(Ypaired)) - 2*cov.G.F
    var.U.mw=  var(Fyprime(Xpaired)) + (1/alpha-1)*var(Fx(Yextra))
    cov.U.mw.p=var(Fy(Xpaired)) - cov.G.F
    
    # MW-MW0, another way to get this is to combine MW-MW1 and MW-MW2
    T = sum(rank(c(Xpaired, Ypaired, Yextra))[(m+1):(m+N)])/(m+N+1)
    var.T.0=(m*N/(m+N))^2 * (1/12/m + 1/12/N - 2*cov.G.F.1/N) 
    var.T=  (m*N/(m+N))^2 * (var(Fyyprime(Xpaired))/m + var(Fx(YYprime))/N - 2*cov.G.F.1/N)
    #var.T=(m*N/(m+N))^2 * (var(Fyyprime(Xpaired))/m + var(Fx(YYprime))/N - 2*cov(Fyyprime(Xpaired), Fx(Ypaired))/N)
    #T.lin=m*N/(m+N)*(-mean(Fyyprime(Xpaired)) + mean(Fx(YYprime))) # linear expansion, note that mean(Fyyprime(Xpaired)) and mean(Fx(YYprime)) are perfectly correlated in data
    
    tests=NULL


       	#X，Y，X'，Y'   	        
        # extension to having both Yextra and Xextra        
        U.3=(sum(rank(c(Xextra, Ypaired))[lx+1:m]) - m*(m+1)/2)/(lx*m)
        #Continuous correction 
        if(correct) U.3=U.3-switch(alternative, "two.sided" = sign(U.3-0) * 0.5/(m*lx), "greater" = 0.5/(m*lx), "less" = -0.5/(m*lx))
        U.4=(sum(rank(c(Xextra, Yextra))[lx+1:n]) - n*(n+1)/2)/(lx*n)
        #Continuous correction 
        if(correct) U.4=U.4-switch(alternative, "two.sided" = sign(U.4-0) * 0.5/(n*lx), "greater" = 0.5/(n*lx), "less" = -0.5/(n*lx))
 
                
       ####################################################
       # SR-MW-fong
       
        vec=sqrt(m)*c((W.plus-(m+1)/4)/m, U.mw-1/2, U.3-1/2, U.4-1/2) 
        v.tmp=diag(c(var.W.plus.0/m, var.U.mw.0, ((m/lx+1)/12), ((m/lx+m/n)/12)))
        v.tmp[2,1]<-v.tmp[1,2]<- -1/2*cov(Fplus(abs.delta) * sgn.delta, Fx(Xpaired))
        v.tmp[3,1]<-v.tmp[1,3]<- 1/2*cov(Fplus(abs.delta) * sgn.delta, Fx(Ypaired))
        v.tmp[4,1]<-v.tmp[1,4]<- 0
        v.tmp[3,2]<-v.tmp[2,3]<- -cov.G.F
        v.tmp[4,2]<-v.tmp[2,4]<- m/(12*n)
        v.tmp[4,3]<-v.tmp[3,4]<- m/(12*lx)                               
        # linear combination
        comb=1/diag(v.tmp)
        #stat.1 = (comb %*% vec)^2 / (comb %*% v.tmp %*% comb)      
        #tests=rbind(tests, "sr.mw.20"=c(stat.1, pchisq(stat.1, df=1, lower.tail=FALSE)))    
        # use normal reference instead of chisq
        stat.1 = (comb %*% vec) / sqrt(comb %*% v.tmp %*% comb)      
        tests=rbind(tests, "sr.mw.fong"=c(stat.1, switch(alternative, "two.sided"=2*pnorm(abs(stat.1),lower.tail=FALSE), "less"=pnorm(stat.1,lower.tail=ifelse(!is.switched,FALSE,TRUE)), "greater"=pnorm(stat.1,lower.tail=ifelse(!is.switched,TRUE,FALSE)))  ))  
                
        
        ####################################################
        # SR-MW-2
        
        vec.nu=sqrt(m)*c((W.plus-(m+1)/4)/m, U.4-1/2) 
        v.tmp=diag(c(var.W.plus.0/m, ((m/lx+m/n)/12)))
        v.tmp[1,2]=v.tmp[2,1]=0;       
        # linear combination
        nu=c(1/v.tmp[1,1]/(1/v.tmp[1,1]+1/v.tmp[2,2]),1/v.tmp[2,2]/(1/v.tmp[1,1]+1/v.tmp[2,2]))
        #if(min(nu)<= 0){D1[q,2]=1}
        #if(min(nu)<= -0.1){D2[q,2]=1}        
        #if(min(nu)<= -0.2){D3[q,2]=1}
        #if(min(nu)<= -0.3){D4[q,2]=1}
        # use normal reference instead of chisq
        stat.opt.2 = (nu %*% vec.nu) / sqrt(nu %*% v.tmp %*% nu)      
        tests=rbind(tests, "sr.mw.2"=c(stat.opt.2, switch(alternative, "two.sided"=2*pnorm(abs(stat.opt.2),lower.tail=FALSE), "less"=pnorm(stat.opt.2,lower.tail=ifelse(!is.switched,FALSE,TRUE)), "greater"=pnorm(stat.opt.2,lower.tail=ifelse(!is.switched,TRUE,FALSE)))  ))
        
  
        ####################################################
        # MW-MW-fong  
          
        vec=sqrt(m)*c(U.p-1/2, U.mw-1/2, U.3-1/2, U.4-1/2)            
        v.tmp=diag(c(var.U.p.0, var.U.mw.0, ((m/lx+1)/12), ((m/lx+m/n)/12)))
        v.tmp[2,1]<-v.tmp[1,2]<- cov.U.mw.p.0
        v.tmp[3,1]<-v.tmp[1,3]<- cov.U.mw.p.0
        v.tmp[4,1]<-v.tmp[1,4]<- 0
        v.tmp[3,2]<-v.tmp[2,3]<- -cov.G.F
        v.tmp[4,2]<-v.tmp[2,4]<- m/(12*n)
        v.tmp[4,3]<-v.tmp[3,4]<- m/(12*lx)
        # linear combination
        comb=1/diag(v.tmp)
        #stat.1 = (comb %*% vec)^2 / (comb %*% v.tmp %*% comb)      
        #tests=rbind(tests, "mw.mw.20"=c(stat.1, pchisq(stat.1, df=1, lower.tail=FALSE)))
        # use normal reference instead of chisq
        stat.1 = (comb %*% vec) / sqrt(comb %*% v.tmp %*% comb)      
        tests=rbind(tests, "mw.mw.fong"=c(stat.1, switch(alternative, "two.sided"=2*pnorm(abs(stat.1),lower.tail=FALSE), "less"=pnorm(stat.1,lower.tail=ifelse(!is.switched,FALSE,TRUE)), "greater"=pnorm(stat.1,lower.tail=ifelse(!is.switched,TRUE,FALSE)))  ))
        

     ####################################################
        # MW-MW-2  
        
        vec.nu=sqrt(m)*c(U.p-1/2, U.4-1/2)            
        v.tmp=diag(c(var.U.p.0, ((m/lx+m/n)/12)))
        v.tmp[2,1]<-v.tmp[1,2]<-0      
        # linear combination
        nu=c(1/v.tmp[1,1]/(1/v.tmp[1,1]+1/v.tmp[2,2]),1/v.tmp[2,2]/(1/v.tmp[1,1]+1/v.tmp[2,2]))
        #if(min(nu)<= 0){D1[q,4]=1}
        #if(min(nu)<= -0.1){D2[q,4]=1}
        #if(min(nu)<= -0.2){D3[q,4]=1}
        #if(min(nu)<= -0.3){D4[q,4]=1}
        # use normal reference instead of chisq
        stat.opt.2 = (nu %*% vec.nu) / sqrt(nu %*% v.tmp %*% nu)      
        tests=rbind(tests, "mw.mw.2"=c(stat.opt.2, switch(alternative, "two.sided"=2*pnorm(abs(stat.opt.2),lower.tail=FALSE), "less"=pnorm(stat.opt.2,lower.tail=ifelse(!is.switched,FALSE,TRUE)), "greater"=pnorm(stat.opt.2,lower.tail=ifelse(!is.switched,TRUE,FALSE)))  ))
           
     C[q,1]=ifelse(tests[1,2]<0.05,1,0)
     C[q,2]=ifelse(tests[2,2]<0.05,1,0)
     C[q,3]=ifelse(tests[3,2]<0.05,1,0)
     C[q,4]=ifelse(tests[4,2]<0.05,1,0)
     
    }
    
    C[is.na(C)]=0
    colnames(C)=c("sr.mw.fong","sr.mw.2","mw.mw.fong","mw.mw.2")
    return(colMeans(C))	
    }
    

    ##################################################################   
    #normal--m=20--delta=0
    D1.1=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=5,n.y=10,loc.2=0) 
    D1.2=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=5,n.y=10,loc.2=0) 
    D1.3=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=5,n.y=10,loc.2=0)
    D1.4=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=10,n.y=10,loc.2=0) 
    D1.5=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=10,n.y=10,loc.2=0) 
    D1.6=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=10,n.y=10,loc.2=0)  
    D1.7=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=5,n.y=40,loc.2=0) 
    D1.8=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=5,n.y=40,loc.2=0) 
    D1.9=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=5,n.y=40,loc.2=0)
    
    #normal--m=20--delta=0.3
    D2.1=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=5,n.y=10,loc.2=0.3) 
    D2.2=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=5,n.y=10,loc.2=0.3) 
    D2.3=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=5,n.y=10,loc.2=0.3)
    D2.4=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=10,n.y=10,loc.2=0.3) 
    D2.5=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=10,n.y=10,loc.2=0.3) 
    D2.6=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=10,n.y=10,loc.2=0.3)  
    D2.7=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=5,n.y=40,loc.2=0.3) 
    D2.8=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=5,n.y=40,loc.2=0.3) 
    D2.9=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=5,n.y=40,loc.2=0.3)
      
    #normal--m=20--delta=0.6
    D3.1=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=5,n.y=10,loc.2=0.6) 
    D3.2=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=5,n.y=10,loc.2=0.6) 
    D3.3=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=5,n.y=10,loc.2=0.6)
    D3.4=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=10,n.y=10,loc.2=0.6) 
    D3.5=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=10,n.y=10,loc.2=0.6) 
    D3.6=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=10,n.y=10,loc.2=0.6)  
    D3.7=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=5,n.y=40,loc.2=0.6) 
    D3.8=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=5,n.y=40,loc.2=0.6) 
    D3.9=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=5,n.y=40,loc.2=0.6)
    
    #normal--m=40--delta=0   
    D4.1=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=5,n.y=10,loc.2=0) 
    D4.2=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=5,n.y=10,loc.2=0) 
    D4.3=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=5,n.y=10,loc.2=0)
    D4.4=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=10,n.y=10,loc.2=0) 
    D4.5=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=10,n.y=10,loc.2=0) 
    D4.6=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=10,n.y=10,loc.2=0)  
    D4.7=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=5,n.y=40,loc.2=0) 
    D4.8=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=5,n.y=40,loc.2=0) 
    D4.9=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=5,n.y=40,loc.2=0)
    
    #normal--m=40--delta=0.3
    D5.1=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=5,n.y=10,loc.2=0.3) 
    D5.2=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=5,n.y=10,loc.2=0.3) 
    D5.3=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=5,n.y=10,loc.2=0.3)
    D5.4=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=10,n.y=10,loc.2=0.3) 
    D5.5=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=10,n.y=10,loc.2=0.3) 
    D5.6=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=10,n.y=10,loc.2=0.3)  
    D5.7=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=5,n.y=40,loc.2=0.3) 
    D5.8=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=5,n.y=40,loc.2=0.3) 
    D5.9=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=5,n.y=40,loc.2=0.3)
       
    #normal--m=40--delta=0.6
    D6.1=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=5,n.y=10,loc.2=0.6) 
    D6.2=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=5,n.y=10,loc.2=0.6) 
    D6.3=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=5,n.y=10,loc.2=0.6)
    D6.4=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=10,n.y=10,loc.2=0.6) 
    D6.5=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=10,n.y=10,loc.2=0.6) 
    D6.6=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=10,n.y=10,loc.2=0.6)  
    D6.7=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=5,n.y=40,loc.2=0.6) 
    D6.8=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=5,n.y=40,loc.2=0.6) 
    D6.9=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=5,n.y=40,loc.2=0.6)

    
    ##################################################################   
    #student--m=20--delta=0(df=3)
    F1.1=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=5,n.y=10,loc.2=0) 
    F1.2=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=5,n.y=10,loc.2=0) 
    F1.3=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=5,n.y=10,loc.2=0)
    F1.4=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=10,n.y=10,loc.2=0) 
    F1.5=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=10,n.y=10,loc.2=0) 
    F1.6=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=10,n.y=10,loc.2=0)  
    F1.7=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=5,n.y=40,loc.2=0) 
    F1.8=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=5,n.y=40,loc.2=0) 
    F1.9=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=5,n.y=40,loc.2=0)
    
    #student--m=20--delta=0.3(df=3)
    F2.1=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=5,n.y=10,loc.2=0.3) 
    F2.2=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=5,n.y=10,loc.2=0.3) 
    F2.3=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=5,n.y=10,loc.2=0.3)
    F2.4=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=10,n.y=10,loc.2=0.3) 
    F2.5=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=10,n.y=10,loc.2=0.3) 
    F2.6=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=10,n.y=10,loc.2=0.3)  
    F2.7=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=5,n.y=40,loc.2=0.3) 
    F2.8=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=5,n.y=40,loc.2=0.3) 
    F2.9=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=5,n.y=40,loc.2=0.3)
        
    #student--m=20--delta=0.6(df=3)
    F3.1=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=5,n.y=10,loc.2=0.6) 
    F3.2=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=5,n.y=10,loc.2=0.6) 
    F3.3=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=5,n.y=10,loc.2=0.6)
    F3.4=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=10,n.y=10,loc.2=0.6) 
    F3.5=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=10,n.y=10,loc.2=0.6) 
    F3.6=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=10,n.y=10,loc.2=0.6)  
    F3.7=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=5,n.y=40,loc.2=0.6) 
    F3.8=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=5,n.y=40,loc.2=0.6) 
    F3.9=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=5,n.y=40,loc.2=0.6)
   
    #student--m=40--delta=0(df=3)
    F4.1=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=5,n.y=10,loc.2=0) 
    F4.2=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=5,n.y=10,loc.2=0) 
    F4.3=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=5,n.y=10,loc.2=0)
    F4.4=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=10,n.y=10,loc.2=0) 
    F4.5=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=10,n.y=10,loc.2=0) 
    F4.6=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=10,n.y=10,loc.2=0)  
    F4.7=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=5,n.y=40,loc.2=0) 
    F4.8=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=5,n.y=40,loc.2=0) 
    F4.9=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=5,n.y=40,loc.2=0)
    
    #student--m=40--delta=0.3(df=3)
    F5.1=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=5,n.y=10,loc.2=0.3) 
    F5.2=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=5,n.y=10,loc.2=0.3) 
    F5.3=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=5,n.y=10,loc.2=0.3)
    F5.4=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=10,n.y=10,loc.2=0.3) 
    F5.5=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=10,n.y=10,loc.2=0.3) 
    F5.6=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=10,n.y=10,loc.2=0.3)  
    F5.7=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=5,n.y=40,loc.2=0.3) 
    F5.8=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=5,n.y=40,loc.2=0.3) 
    F5.9=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=5,n.y=40,loc.2=0.3)
        
    #student--m=40--delta=0.6(df=3)
    F6.1=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=5,n.y=10,loc.2=0.6) 
    F6.2=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=5,n.y=10,loc.2=0.6) 
    F6.3=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=5,n.y=10,loc.2=0.6)
    F6.4=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=10,n.y=10,loc.2=0.6) 
    F6.5=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=10,n.y=10,loc.2=0.6) 
    F6.6=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=10,n.y=10,loc.2=0.6)  
    F6.7=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=5,n.y=40,loc.2=0.6) 
    F6.8=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=5,n.y=40,loc.2=0.6) 
    F6.9=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=5,n.y=40,loc.2=0.6)
 
    
    ##################################################################   
    #logistic--m=20--delta=0
    H1.1=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=5,n.y=10,loc.2=0) 
    H1.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=5,n.y=10,loc.2=0) 
    H1.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=5,n.y=10,loc.2=0)
    H1.4=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=10,n.y=10,loc.2=0) 
    H1.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=10,n.y=10,loc.2=0) 
    H1.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=10,n.y=10,loc.2=0)  
    H1.7=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=5,n.y=40,loc.2=0) 
    H1.8=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=5,n.y=40,loc.2=0) 
    H1.9=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=5,n.y=40,loc.2=0)
    
    #logistic--m=20--delta=0.3
    H2.1=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=5,n.y=10,loc.2=0.3) 
    H2.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=5,n.y=10,loc.2=0.3) 
    H2.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=5,n.y=10,loc.2=0.3)
    H2.4=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=10,n.y=10,loc.2=0.3) 
    H2.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=10,n.y=10,loc.2=0.3) 
    H2.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=10,n.y=10,loc.2=0.3)  
    H2.7=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=5,n.y=40,loc.2=0.3) 
    H2.8=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=5,n.y=40,loc.2=0.3) 
    H2.9=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=5,n.y=40,loc.2=0.3)
      
    #logistic--m=20--delta=0.6
    H3.1=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=5,n.y=10,loc.2=0.6) 
    H3.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=5,n.y=10,loc.2=0.6) 
    H3.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=5,n.y=10,loc.2=0.6)
    H3.4=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=10,n.y=10,loc.2=0.6) 
    H3.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=10,n.y=10,loc.2=0.6) 
    H3.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=10,n.y=10,loc.2=0.6)  
    H3.7=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=5,n.y=40,loc.2=0.6) 
    H3.8=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=5,n.y=40,loc.2=0.6) 
    H3.9=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=5,n.y=40,loc.2=0.6)
    
    #logistic--m=40--delta=0   
    H4.1=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=5,n.y=10,loc.2=0) 
    H4.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=5,n.y=10,loc.2=0) 
    H4.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=5,n.y=10,loc.2=0)
    H4.4=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=10,n.y=10,loc.2=0) 
    H4.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=10,n.y=10,loc.2=0) 
    H4.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=10,n.y=10,loc.2=0)  
    H4.7=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=5,n.y=40,loc.2=0) 
    H4.8=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=5,n.y=40,loc.2=0) 
    H4.9=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=5,n.y=40,loc.2=0)
    
    #logistic--m=40--delta=0.3
    H5.1=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=5,n.y=10,loc.2=0.3) 
    H5.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=5,n.y=10,loc.2=0.3) 
    H5.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=5,n.y=10,loc.2=0.3)
    H5.4=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=10,n.y=10,loc.2=0.3) 
    H5.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=10,n.y=10,loc.2=0.3) 
    H5.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=10,n.y=10,loc.2=0.3)  
    H5.7=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=5,n.y=40,loc.2=0.3) 
    H5.8=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=5,n.y=40,loc.2=0.3) 
    H5.9=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=5,n.y=40,loc.2=0.3)
       
    #logistic--m=40--delta=0.6
    H6.1=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=5,n.y=10,loc.2=0.6) 
    H6.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=5,n.y=10,loc.2=0.6) 
    H6.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=5,n.y=10,loc.2=0.6)
    H6.4=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=10,n.y=10,loc.2=0.6) 
    H6.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=10,n.y=10,loc.2=0.6) 
    H6.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=10,n.y=10,loc.2=0.6)  
    H6.7=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=5,n.y=40,loc.2=0.6) 
    H6.8=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=5,n.y=40,loc.2=0.6) 
    H6.9=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=5,n.y=40,loc.2=0.6)

    
    