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
    	B=matrix(0,nrow=N,ncol=4);#records of the results
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

          #X，Y，Y'          
        
         ####################################################
        # SR-MW tests, several versions of rho are tried
        
        #vec=c(W.plus-(m+1)/4, W.mw-n/2)
        #v.tmp=matrix(c(var.W.plus.0, cov.0, cov.0, var.W.mw.0),2)
        # linear combination
        #df.=Inf # degrees of freedom for Student's t: m or Inf
        #comb=c(1, var.W.plus.0/var.W.mw.0)
        #stat.1 = (comb %*% vec) / sqrt(comb %*% v.tmp %*% comb)      
        #tests=rbind(tests, "sr.mw.20"=c(stat.1, switch(alternative, "two.sided"=2*pt(abs(stat.1),df=df.,lower.tail=FALSE),  "less"=pt(stat.1,df=df.,lower.tail=ifelse(!is.switched,FALSE,TRUE)), "greater"=pt(stat.1,df=df.,lower.tail=ifelse(!is.switched,TRUE,FALSE)))  ))

        
        ####################################################
        # SR-MW tests, several versions of rho are tried
        # changed into W.plus and U.mw
        
        vec=sqrt(m)*c((W.plus-(m+1)/4)/m, U.mw-1/2) 
        v.tmp=matrix(c(var.W.plus.0/m,-1/2*cov(Fplus(abs.delta) * sgn.delta, Fx(Xpaired)),-1/2*cov(Fplus(abs.delta) * sgn.delta, Fx(Xpaired)), var.U.mw.0),2)                 
        # linear combination
        df.=Inf # degrees of freedom for Student's t: m or Inf
        comb=c(1, var.W.plus.0/m/var.U.mw.0)
        stat.1 =(comb %*% vec) / sqrt(comb %*% v.tmp %*% comb)      
        tests=rbind(tests, "sr.mw.fong"=c(stat.1, switch(alternative,"two.sided"=2*pt(abs(stat.1),df=df.,lower.tail=FALSE),"less"=pt(stat.1,df=df.,lower.tail=ifelse(!is.switched,FALSE,TRUE)),"greater"=pt(stat.1,df=df.,lower.tail=ifelse(!is.switched,TRUE,FALSE)))))              

        
        ####################################################
        # SR-MW-opt tests, several versions of rho are tried
            
        omega.0=c((v.tmp[2,2]-v.tmp[1,2])/(v.tmp[1,1]+v.tmp[2,2]-2*v.tmp[1,2]),
        (v.tmp[1,1]-v.tmp[1,2])/(v.tmp[1,1]+v.tmp[2,2]-2*v.tmp[1,2]))
        #if(min(omega.0)<0){D[q,1]=1}
        stat.opt.0=(omega.0 %*% vec)/ sqrt(omega.0 %*% v.tmp %*% omega.0)      
        tests=rbind(tests, "sr.mw.opt"=c(stat.opt.0, switch(alternative, "two.sided"=2*pt(abs(stat.opt.0),df=df.,lower.tail=FALSE), "less"=pt(stat.opt.0,df=df.,lower.tail=ifelse(!is.switched,FALSE,TRUE)), "greater"=pt(stat.opt.0,df=df.,lower.tail=ifelse(!is.switched,TRUE,FALSE)))))
        #use t-distribution can improve power
        #tests=rbind(tests, "sr.mw.opt"=c(stat.opt.0, switch(alternative, "two.sided"=2*pnorm(abs(stat.opt.0),lower.tail=FALSE), "less"=pnorm(stat.opt.0,lower.tail=ifelse(!is.switched,FALSE,TRUE)), "greater"=pnorm(stat.opt.0,lower.tail=ifelse(!is.switched,TRUE,FALSE)))  ))       
  
  
        ####################################################
         # MW-MW tests 
                  
       vec=sqrt(m)*c(U.p-1/2, U.mw-1/2) 
       v.tmp=matrix(c(var.U.p.0, cov.U.mw.p.0, cov.U.mw.p.0, var.U.mw.0),2) 
       # linear combination
        comb.0=c(1, var.U.p.0/var.U.mw.0)
        stat.1 = (comb.0 %*% vec) / sqrt(comb.0 %*% v.tmp %*% comb.0)      
        tests=rbind(tests, "mw.mw.fong"=c(stat.1, switch(alternative, "two.sided"=2*pnorm(abs(stat.1),lower.tail=FALSE), "less"=pnorm(stat.1,lower.tail=ifelse(!is.switched,FALSE,TRUE)), "greater"=pnorm(stat.1,lower.tail=ifelse(!is.switched,TRUE,FALSE)))  ))
 
             
        ####################################################
        # MW-MW-opt tests
         
        vec=sqrt(m)*c(U.p-1/2, U.mw-1/2)                        
        # optimal linear combination
        omega.0=c((v.tmp[2,2]-v.tmp[1,2])/(v.tmp[1,1]+v.tmp[2,2]-2*v.tmp[1,2]),
        (v.tmp[1,1]-v.tmp[1,2])/(v.tmp[1,1]+v.tmp[2,2]-2*v.tmp[1,2]))
        #if(min(omega.0)<0){D[q,2]=1}
        stat.opt.0=(omega.0 %*% vec) / sqrt(omega.0 %*% v.tmp %*% omega.0)      
        tests=rbind(tests, "mw.mw.opt"=c(stat.opt.0, switch(alternative, "two.sided"=2*pnorm(abs(stat.opt.0),lower.tail=FALSE), "less"=pnorm(stat.opt.0,lower.tail=ifelse(!is.switched,FALSE,TRUE)), "greater"=pnorm(stat.opt.0,lower.tail=ifelse(!is.switched,TRUE,FALSE)))  ))   
  

     B[q,1]=ifelse(tests[1,2]<0.05,1,0)
     B[q,2]=ifelse(tests[2,2]<0.05,1,0)
     B[q,3]=ifelse(tests[3,2]<0.05,1,0)
     B[q,4]=ifelse(tests[4,2]<0.05,1,0)       
    }
    
    B[is.na(B)]=0
    colnames(B)=c("sr.mw.fong","sr.mw.opt","mw.mw.fong","mw.mw.opt")
    return(colMeans(B))
    }
 
    #small sample 
    ##################################################################      
    #normal--m=5--delta=0
    A1.1=ranktests(10000,"two.sided","normal",rho=0,m=5,n.x=0,n.y=5,loc.2=0) 
    A1.2=ranktests(10000,"two.sided","normal",rho=0.5,m=5,n.x=0,n.y=5,loc.2=0) 
    A1.3=ranktests(10000,"two.sided","normal",rho=0.8,m=5,n.x=0,n.y=5,loc.2=0)
    A1.4=ranktests(10000,"two.sided","normal",rho=0,m=5,n.x=0,n.y=10,loc.2=0) 
    A1.5=ranktests(10000,"two.sided","normal",rho=0.5,m=5,n.x=0,n.y=10,loc.2=0) 
    A1.6=ranktests(10000,"two.sided","normal",rho=0.8,m=5,n.x=0,n.y=10,loc.2=0)
    
    #normal--m=5--delta=0.3
    A2.1=ranktests(10000,"two.sided","normal",rho=0,m=5,n.x=0,n.y=5,loc.2=0.3) 
    A2.2=ranktests(10000,"two.sided","normal",rho=0.5,m=5,n.x=0,n.y=5,loc.2=0.3) 
    A2.3=ranktests(10000,"two.sided","normal",rho=0.8,m=5,n.x=0,n.y=5,loc.2=0.3)
    A2.4=ranktests(10000,"two.sided","normal",rho=0,m=5,n.x=0,n.y=10,loc.2=0.3) 
    A2.5=ranktests(10000,"two.sided","normal",rho=0.5,m=5,n.x=0,n.y=10,loc.2=0.3) 
    A2.6=ranktests(10000,"two.sided","normal",rho=0.8,m=5,n.x=0,n.y=10,loc.2=0.3)  
    
    #normal--m=5--delta=0.6
    A3.1=ranktests(10000,"two.sided","normal",rho=0,m=5,n.x=0,n.y=5,loc.2=0.6) 
    A3.2=ranktests(10000,"two.sided","normal",rho=0.5,m=5,n.x=0,n.y=5,loc.2=0.6) 
    A3.3=ranktests(10000,"two.sided","normal",rho=0.8,m=5,n.x=0,n.y=5,loc.2=0.6)
    A3.4=ranktests(10000,"two.sided","normal",rho=0,m=5,n.x=0,n.y=10,loc.2=0.6) 
    A3.5=ranktests(10000,"two.sided","normal",rho=0.5,m=5,n.x=0,n.y=10,loc.2=0.6) 
    A3.6=ranktests(10000,"two.sided","normal",rho=0.8,m=5,n.x=0,n.y=10,loc.2=0.6)
    
    #normal--m=10--delta=0
    A4.1=ranktests(10000,"two.sided","normal",rho=0,m=10,n.x=0,n.y=5,loc.2=0) 
    A4.2=ranktests(10000,"two.sided","normal",rho=0.5,m=10,n.x=0,n.y=5,loc.2=0) 
    A4.3=ranktests(10000,"two.sided","normal",rho=0.8,m=10,n.x=0,n.y=5,loc.2=0)
    A4.4=ranktests(10000,"two.sided","normal",rho=0,m=10,n.x=0,n.y=10,loc.2=0) 
    A4.5=ranktests(10000,"two.sided","normal",rho=0.5,m=10,n.x=0,n.y=10,loc.2=0) 
    A4.6=ranktests(10000,"two.sided","normal",rho=0.8,m=10,n.x=0,n.y=10,loc.2=0)
    A4.7=ranktests(10000,"two.sided","normal",rho=0,m=10,n.x=0,n.y=20,loc.2=0) 
    A4.8=ranktests(10000,"two.sided","normal",rho=0.5,m=10,n.x=0,n.y=20,loc.2=0) 
    A4.9=ranktests(10000,"two.sided","normal",rho=0.8,m=10,n.x=0,n.y=20,loc.2=0)
       
    #normal--m=10--delta=0.3
    A5.1=ranktests(10000,"two.sided","normal",rho=0,m=10,n.x=0,n.y=5,loc.2=0.3) 
    A5.2=ranktests(10000,"two.sided","normal",rho=0.5,m=10,n.x=0,n.y=5,loc.2=0.3) 
    A5.3=ranktests(10000,"two.sided","normal",rho=0.8,m=10,n.x=0,n.y=5,loc.2=0.3)
    A5.4=ranktests(10000,"two.sided","normal",rho=0,m=10,n.x=0,n.y=10,loc.2=0.3) 
    A5.5=ranktests(10000,"two.sided","normal",rho=0.5,m=10,n.x=0,n.y=10,loc.2=0.3) 
    A5.6=ranktests(10000,"two.sided","normal",rho=0.8,m=10,n.x=0,n.y=10,loc.2=0.3)    
    A5.7=ranktests(10000,"two.sided","normal",rho=0,m=10,n.x=0,n.y=20,loc.2=0.3) 
    A5.8=ranktests(10000,"two.sided","normal",rho=0.5,m=10,n.x=0,n.y=20,loc.2=0.3) 
    A5.9=ranktests(10000,"two.sided","normal",rho=0.8,m=10,n.x=0,n.y=20,loc.2=0.3) 
    
    #normal--m=10--delta=0.6
    A6.1=ranktests(10000,"two.sided","normal",rho=0,m=10,n.x=0,n.y=5,loc.2=0.6) 
    A6.2=ranktests(10000,"two.sided","normal",rho=0.5,m=10,n.x=0,n.y=5,loc.2=0.6) 
    A6.3=ranktests(10000,"two.sided","normal",rho=0.8,m=10,n.x=0,n.y=5,loc.2=0.6)
    A6.4=ranktests(10000,"two.sided","normal",rho=0,m=10,n.x=0,n.y=10,loc.2=0.6) 
    A6.5=ranktests(10000,"two.sided","normal",rho=0.5,m=10,n.x=0,n.y=10,loc.2=0.6) 
    A6.6=ranktests(10000,"two.sided","normal",rho=0.8,m=10,n.x=0,n.y=10,loc.2=0.6)
    A6.7=ranktests(10000,"two.sided","normal",rho=0,m=10,n.x=0,n.y=20,loc.2=0.6) 
    A6.8=ranktests(10000,"two.sided","normal",rho=0.5,m=10,n.x=0,n.y=20,loc.2=0.6) 
    A6.9=ranktests(10000,"two.sided","normal",rho=0.8,m=10,n.x=0,n.y=20,loc.2=0.6)
   
 
    ##################################################################     
    #lognormal--m=5--delta=0
    B1.1=ranktests(10000,"two.sided","lognormal",rho=0,m=5,n.x=0,n.y=5,loc.2=0) 
    B1.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=5,n.x=0,n.y=5,loc.2=0) 
    B1.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=5,n.x=0,n.y=5,loc.2=0)
    B1.4=ranktests(10000,"two.sided","lognormal",rho=0,m=5,n.x=0,n.y=10,loc.2=0) 
    B1.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=5,n.x=0,n.y=10,loc.2=0) 
    B1.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=5,n.x=0,n.y=10,loc.2=0)
    
    #normal--m=5--delta=0.3
    B2.1=ranktests(10000,"two.sided","lognormal",rho=0,m=5,n.x=0,n.y=5,loc.2=0.3) 
    B2.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=5,n.x=0,n.y=5,loc.2=0.3) 
    B2.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=5,n.x=0,n.y=5,loc.2=0.3)
    B2.4=ranktests(10000,"two.sided","lognormal",rho=0,m=5,n.x=0,n.y=10,loc.2=0.3) 
    B2.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=5,n.x=0,n.y=10,loc.2=0.3) 
    B2.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=5,n.x=0,n.y=10,loc.2=0.3)  
    
    #normal--m=5--delta=0.6
    B3.1=ranktests(10000,"two.sided","lognormal",rho=0,m=5,n.x=0,n.y=5,loc.2=0.6) 
    B3.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=5,n.x=0,n.y=5,loc.2=0.6) 
    B3.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=5,n.x=0,n.y=5,loc.2=0.6)
    B3.4=ranktests(10000,"two.sided","lognormal",rho=0,m=5,n.x=0,n.y=10,loc.2=0.6) 
    B3.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=5,n.x=0,n.y=10,loc.2=0.6) 
    B3.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=5,n.x=0,n.y=10,loc.2=0.6)
    
    #normal--m=10--delta=0
    B4.1=ranktests(10000,"two.sided","lognormal",rho=0,m=10,n.x=0,n.y=5,loc.2=0) 
    B4.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=10,n.x=0,n.y=5,loc.2=0) 
    B4.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=10,n.x=0,n.y=5,loc.2=0)
    B4.4=ranktests(10000,"two.sided","lognormal",rho=0,m=10,n.x=0,n.y=10,loc.2=0) 
    B4.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=10,n.x=0,n.y=10,loc.2=0) 
    B4.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=10,n.x=0,n.y=10,loc.2=0)
    B4.7=ranktests(10000,"two.sided","lognormal",rho=0,m=10,n.x=0,n.y=20,loc.2=0) 
    B4.8=ranktests(10000,"two.sided","lognormal",rho=0.5,m=10,n.x=0,n.y=20,loc.2=0) 
    B4.9=ranktests(10000,"two.sided","lognormal",rho=0.8,m=10,n.x=0,n.y=20,loc.2=0)
       
    #normal--m=10--delta=0.3
    B5.1=ranktests(10000,"two.sided","lognormal",rho=0,m=10,n.x=0,n.y=5,loc.2=0.3) 
    B5.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=10,n.x=0,n.y=5,loc.2=0.3) 
    B5.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=10,n.x=0,n.y=5,loc.2=0.3)
    B5.4=ranktests(10000,"two.sided","lognormal",rho=0,m=10,n.x=0,n.y=10,loc.2=0.3) 
    B5.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=10,n.x=0,n.y=10,loc.2=0.3) 
    B5.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=10,n.x=0,n.y=10,loc.2=0.3)    
    B5.7=ranktests(10000,"two.sided","lognormal",rho=0,m=10,n.x=0,n.y=20,loc.2=0.3) 
    B5.8=ranktests(10000,"two.sided","lognormal",rho=0.5,m=10,n.x=0,n.y=20,loc.2=0.3) 
    B5.9=ranktests(10000,"two.sided","lognormal",rho=0.8,m=10,n.x=0,n.y=20,loc.2=0.3) 
    
    #normal--m=10--delta=0.6
    B6.1=ranktests(10000,"two.sided","lognormal",rho=0,m=10,n.x=0,n.y=5,loc.2=0.6) 
    B6.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=10,n.x=0,n.y=5,loc.2=0.6) 
    B6.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=10,n.x=0,n.y=5,loc.2=0.6)
    B6.4=ranktests(10000,"two.sided","lognormal",rho=0,m=10,n.x=0,n.y=10,loc.2=0.6) 
    B6.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=10,n.x=0,n.y=10,loc.2=0.6) 
    B6.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=10,n.x=0,n.y=10,loc.2=0.6)
    B6.7=ranktests(10000,"two.sided","lognormal",rho=0,m=10,n.x=0,n.y=20,loc.2=0.6) 
    B6.8=ranktests(10000,"two.sided","lognormal",rho=0.5,m=10,n.x=0,n.y=20,loc.2=0.6) 
    B6.9=ranktests(10000,"two.sided","lognormal",rho=0.8,m=10,n.x=0,n.y=20,loc.2=0.6)
   

    ##################################################################       
    #student--m=5--delta=0(df=3)
    C1.1=ranktests(10000,"two.sided","student",rho=0,m=5,n.x=0,n.y=5,loc.2=0) 
    C1.2=ranktests(10000,"two.sided","student",rho=0.5,m=5,n.x=0,n.y=5,loc.2=0) 
    C1.3=ranktests(10000,"two.sided","student",rho=0.8,m=5,n.x=0,n.y=5,loc.2=0)
    C1.4=ranktests(10000,"two.sided","student",rho=0,m=5,n.x=0,n.y=10,loc.2=0) 
    C1.5=ranktests(10000,"two.sided","student",rho=0.5,m=5,n.x=0,n.y=10,loc.2=0) 
    C1.6=ranktests(10000,"two.sided","student",rho=0.8,m=5,n.x=0,n.y=10,loc.2=0)
   
    #student--m=5--delta=0.3(df=3)
    C2.1=ranktests(10000,"two.sided","student",rho=0,m=5,n.x=0,n.y=5,loc.2=0.3) 
    C2.2=ranktests(10000,"two.sided","student",rho=0.5,m=5,n.x=0,n.y=5,loc.2=0.3) 
    C2.3=ranktests(10000,"two.sided","student",rho=0.8,m=5,n.x=0,n.y=5,loc.2=0.3)
    C2.4=ranktests(10000,"two.sided","student",rho=0,m=5,n.x=0,n.y=10,loc.2=0.3) 
    C2.5=ranktests(10000,"two.sided","student",rho=0.5,m=5,n.x=0,n.y=10,loc.2=0.3) 
    C2.6=ranktests(10000,"two.sided","student",rho=0.8,m=5,n.x=0,n.y=10,loc.2=0.3)  
    
    #student--m=5--delta=0.6(df=3)
    C3.1=ranktests(10000,"two.sided","student",rho=0,m=5,n.x=0,n.y=5,loc.2=0.6) 
    C3.2=ranktests(10000,"two.sided","student",rho=0.5,m=5,n.x=0,n.y=5,loc.2=0.6) 
    C3.3=ranktests(10000,"two.sided","student",rho=0.8,m=5,n.x=0,n.y=5,loc.2=0.6)
    C3.4=ranktests(10000,"two.sided","student",rho=0,m=5,n.x=0,n.y=10,loc.2=0.6) 
    C3.5=ranktests(10000,"two.sided","student",rho=0.5,m=5,n.x=0,n.y=10,loc.2=0.6) 
    C3.6=ranktests(10000,"two.sided","student",rho=0.8,m=5,n.x=0,n.y=10,loc.2=0.6)
    
    #student--m=10--delta=0(df=3)
    C4.1=ranktests(10000,"two.sided","student",rho=0,m=10,n.x=0,n.y=5,loc.2=0) 
    C4.2=ranktests(10000,"two.sided","student",rho=0.5,m=10,n.x=0,n.y=5,loc.2=0) 
    C4.3=ranktests(10000,"two.sided","student",rho=0.8,m=10,n.x=0,n.y=5,loc.2=0)
    C4.4=ranktests(10000,"two.sided","student",rho=0,m=10,n.x=0,n.y=10,loc.2=0) 
    C4.5=ranktests(10000,"two.sided","student",rho=0.5,m=10,n.x=0,n.y=10,loc.2=0) 
    C4.6=ranktests(10000,"two.sided","student",rho=0.8,m=10,n.x=0,n.y=10,loc.2=0)
    C4.7=ranktests(10000,"two.sided","student",rho=0,m=10,n.x=0,n.y=20,loc.2=0) 
    C4.8=ranktests(10000,"two.sided","student",rho=0.5,m=10,n.x=0,n.y=20,loc.2=0) 
    C4.9=ranktests(10000,"two.sided","student",rho=0.8,m=10,n.x=0,n.y=20,loc.2=0)
       
    #student--m=10--delta=0.3(df=3)
    C5.1=ranktests(10000,"two.sided","student",rho=0,m=10,n.x=0,n.y=5,loc.2=0.3) 
    C5.2=ranktests(10000,"two.sided","student",rho=0.5,m=10,n.x=0,n.y=5,loc.2=0.3) 
    C5.3=ranktests(10000,"two.sided","student",rho=0.8,m=10,n.x=0,n.y=5,loc.2=0.3)
    C5.4=ranktests(10000,"two.sided","student",rho=0,m=10,n.x=0,n.y=10,loc.2=0.3) 
    C5.5=ranktests(10000,"two.sided","student",rho=0.5,m=10,n.x=0,n.y=10,loc.2=0.3) 
    C5.6=ranktests(10000,"two.sided","student",rho=0.8,m=10,n.x=0,n.y=10,loc.2=0.3)    
    C5.7=ranktests(10000,"two.sided","student",rho=0,m=10,n.x=0,n.y=20,loc.2=0.3) 
    C5.8=ranktests(10000,"two.sided","student",rho=0.5,m=10,n.x=0,n.y=20,loc.2=0.3) 
    C5.9=ranktests(10000,"two.sided","student",rho=0.8,m=10,n.x=0,n.y=20,loc.2=0.3) 
    
    #student--m=10--delta=0.6(df=3)
    C6.1=ranktests(10000,"two.sided","student",rho=0,m=10,n.x=0,n.y=5,loc.2=0.6) 
    C6.2=ranktests(10000,"two.sided","student",rho=0.5,m=10,n.x=0,n.y=5,loc.2=0.6) 
    C6.3=ranktests(10000,"two.sided","student",rho=0.8,m=10,n.x=0,n.y=5,loc.2=0.6)
    C6.4=ranktests(10000,"two.sided","student",rho=0,m=10,n.x=0,n.y=10,loc.2=0.6) 
    C6.5=ranktests(10000,"two.sided","student",rho=0.5,m=10,n.x=0,n.y=10,loc.2=0.6) 
    C6.6=ranktests(10000,"two.sided","student",rho=0.8,m=10,n.x=0,n.y=10,loc.2=0.6)
    C6.7=ranktests(10000,"two.sided","student",rho=0,m=10,n.x=0,n.y=20,loc.2=0.6) 
    C6.8=ranktests(10000,"two.sided","student",rho=0.5,m=10,n.x=0,n.y=20,loc.2=0.6) 
    C6.9=ranktests(10000,"two.sided","student",rho=0.8,m=10,n.x=0,n.y=20,loc.2=0.6)
    
    ##################################################################       
    #logistic--m=5--delta=0
    G1.1=ranktests(10000,"two.sided","logistic",rho=0,m=5,n.x=0,n.y=5,loc.2=0) 
    G1.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=5,n.x=0,n.y=5,loc.2=0) 
    G1.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=5,n.x=0,n.y=5,loc.2=0)
    G1.4=ranktests(10000,"two.sided","logistic",rho=0,m=5,n.x=0,n.y=10,loc.2=0) 
    G1.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=5,n.x=0,n.y=10,loc.2=0) 
    G1.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=5,n.x=0,n.y=10,loc.2=0)
   
    #logistic--m=5--delta=0.3
    G2.1=ranktests(10000,"two.sided","logistic",rho=0,m=5,n.x=0,n.y=5,loc.2=0.3) 
    G2.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=5,n.x=0,n.y=5,loc.2=0.3) 
    G2.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=5,n.x=0,n.y=5,loc.2=0.3)
    G2.4=ranktests(10000,"two.sided","logistic",rho=0,m=5,n.x=0,n.y=10,loc.2=0.3) 
    G2.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=5,n.x=0,n.y=10,loc.2=0.3) 
    G2.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=5,n.x=0,n.y=10,loc.2=0.3)  
    
    #logistic--m=5--delta=0.6
    G3.1=ranktests(10000,"two.sided","logistic",rho=0,m=5,n.x=0,n.y=5,loc.2=0.6) 
    G3.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=5,n.x=0,n.y=5,loc.2=0.6) 
    G3.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=5,n.x=0,n.y=5,loc.2=0.6)
    G3.4=ranktests(10000,"two.sided","logistic",rho=0,m=5,n.x=0,n.y=10,loc.2=0.6) 
    G3.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=5,n.x=0,n.y=10,loc.2=0.6) 
    G3.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=5,n.x=0,n.y=10,loc.2=0.6)
    
    #logistic--m=10--delta=0
    G4.1=ranktests(10000,"two.sided","logistic",rho=0,m=10,n.x=0,n.y=5,loc.2=0) 
    G4.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=10,n.x=0,n.y=5,loc.2=0) 
    G4.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=10,n.x=0,n.y=5,loc.2=0)
    G4.4=ranktests(10000,"two.sided","logistic",rho=0,m=10,n.x=0,n.y=10,loc.2=0) 
    G4.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=10,n.x=0,n.y=10,loc.2=0) 
    G4.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=10,n.x=0,n.y=10,loc.2=0)
    G4.7=ranktests(10000,"two.sided","logistic",rho=0,m=10,n.x=0,n.y=20,loc.2=0) 
    G4.8=ranktests(10000,"two.sided","logistic",rho=0.5,m=10,n.x=0,n.y=20,loc.2=0) 
    G4.9=ranktests(10000,"two.sided","logistic",rho=0.8,m=10,n.x=0,n.y=20,loc.2=0)
       
    #logistic--m=10--delta=0.3
    G5.1=ranktests(10000,"two.sided","logistic",rho=0,m=10,n.x=0,n.y=5,loc.2=0.3) 
    G5.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=10,n.x=0,n.y=5,loc.2=0.3) 
    G5.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=10,n.x=0,n.y=5,loc.2=0.3)
    G5.4=ranktests(10000,"two.sided","logistic",rho=0,m=10,n.x=0,n.y=10,loc.2=0.3) 
    G5.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=10,n.x=0,n.y=10,loc.2=0.3) 
    G5.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=10,n.x=0,n.y=10,loc.2=0.3)    
    G5.7=ranktests(10000,"two.sided","logistic",rho=0,m=10,n.x=0,n.y=20,loc.2=0.3) 
    G5.8=ranktests(10000,"two.sided","logistic",rho=0.5,m=10,n.x=0,n.y=20,loc.2=0.3) 
    G5.9=ranktests(10000,"two.sided","logistic",rho=0.8,m=10,n.x=0,n.y=20,loc.2=0.3) 
    
    #logistic--m=10--delta=0.6
    G6.1=ranktests(10000,"two.sided","logistic",rho=0,m=10,n.x=0,n.y=5,loc.2=0.6) 
    G6.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=10,n.x=0,n.y=5,loc.2=0.6) 
    G6.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=10,n.x=0,n.y=5,loc.2=0.6)
    G6.4=ranktests(10000,"two.sided","logistic",rho=0,m=10,n.x=0,n.y=10,loc.2=0.6) 
    G6.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=10,n.x=0,n.y=10,loc.2=0.6) 
    G6.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=10,n.x=0,n.y=10,loc.2=0.6)
    G6.7=ranktests(10000,"two.sided","logistic",rho=0,m=10,n.x=0,n.y=20,loc.2=0.6) 
    G6.8=ranktests(10000,"two.sided","logistic",rho=0.5,m=10,n.x=0,n.y=20,loc.2=0.6) 
    G6.9=ranktests(10000,"two.sided","logistic",rho=0.8,m=10,n.x=0,n.y=20,loc.2=0.6)
    
    #big sample
    ##################################################################
    #normal--m=20--delta=0
    D1.1=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=0,n.y=10,loc.2=0) 
    D1.2=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=0,n.y=10,loc.2=0) 
    D1.3=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=0,n.y=10,loc.2=0)
    D1.4=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=0,n.y=40,loc.2=0) 
    D1.5=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=0,n.y=40,loc.2=0) 
    D1.6=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=0,n.y=40,loc.2=0)
    
    #normal--m=20--delta=0.3
    D2.1=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=0,n.y=10,loc.2=0.3) 
    D2.2=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=0,n.y=10,loc.2=0.3) 
    D2.3=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=0,n.y=10,loc.2=0.3)
    D2.4=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=0,n.y=40,loc.2=0.3) 
    D2.5=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=0,n.y=40,loc.2=0.3) 
    D2.6=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=0,n.y=40,loc.2=0.3)  
    
    #normal--m=20--delta=0.6
    D3.1=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=0,n.y=10,loc.2=0.6) 
    D3.2=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=0,n.y=10,loc.2=0.6) 
    D3.3=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=0,n.y=10,loc.2=0.6)
    D3.4=ranktests(10000,"two.sided","normal",rho=0,m=20,n.x=0,n.y=40,loc.2=0.6) 
    D3.5=ranktests(10000,"two.sided","normal",rho=0.5,m=20,n.x=0,n.y=40,loc.2=0.6) 
    D3.6=ranktests(10000,"two.sided","normal",rho=0.8,m=20,n.x=0,n.y=40,loc.2=0.6)
    
    #normal--m=40--delta=0
    D4.1=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=0,n.y=10,loc.2=0) 
    D4.2=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=0,n.y=10,loc.2=0) 
    D4.3=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=0,n.y=10,loc.2=0)
    D4.4=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=0,n.y=40,loc.2=0) 
    D4.5=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=0,n.y=40,loc.2=0) 
    D4.6=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=0,n.y=40,loc.2=0)
       
    #normal--m=40--delta=0.3
    D5.1=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=0,n.y=10,loc.2=0.3) 
    D5.2=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=0,n.y=10,loc.2=0.3) 
    D5.3=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=0,n.y=10,loc.2=0.3)
    D5.4=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=0,n.y=40,loc.2=0.3) 
    D5.5=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=0,n.y=40,loc.2=0.3) 
    D5.6=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=0,n.y=40,loc.2=0.3)    
    
    #normal--m=40--delta=0.6
    D6.1=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=0,n.y=10,loc.2=0.6) 
    D6.2=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=0,n.y=10,loc.2=0.6) 
    D6.3=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=0,n.y=10,loc.2=0.6)
    D6.4=ranktests(10000,"two.sided","normal",rho=0,m=40,n.x=0,n.y=40,loc.2=0.6) 
    D6.5=ranktests(10000,"two.sided","normal",rho=0.5,m=40,n.x=0,n.y=40,loc.2=0.6) 
    D6.6=ranktests(10000,"two.sided","normal",rho=0.8,m=40,n.x=0,n.y=40,loc.2=0.6)
    
    #normal--m=80--delta=0
    D7.1=ranktests(10000,"two.sided","normal",rho=0,m=80,n.x=0,n.y=10,loc.2=0) 
    D7.2=ranktests(10000,"two.sided","normal",rho=0.5,m=80,n.x=0,n.y=10,loc.2=0) 
    D7.3=ranktests(10000,"two.sided","normal",rho=0.8,m=80,n.x=0,n.y=10,loc.2=0)
    D7.4=ranktests(10000,"two.sided","normal",rho=0,m=80,n.x=0,n.y=40,loc.2=0) 
    D7.5=ranktests(10000,"two.sided","normal",rho=0.5,m=80,n.x=0,n.y=40,loc.2=0) 
    D7.6=ranktests(10000,"two.sided","normal",rho=0.8,m=80,n.x=0,n.y=40,loc.2=0)
    D7.7=ranktests(10000,"two.sided","normal",rho=0,m=80,n.x=0,n.y=80,loc.2=0) 
    D7.8=ranktests(10000,"two.sided","normal",rho=0.5,m=80,n.x=0,n.y=80,loc.2=0) 
    D7.9=ranktests(10000,"two.sided","normal",rho=0.8,m=80,n.x=0,n.y=80,loc.2=0)
       
    #normal--m=80--delta=0.3
    D8.1=ranktests(10000,"two.sided","normal",rho=0,m=80,n.x=0,n.y=10,loc.2=0.3) 
    D8.2=ranktests(10000,"two.sided","normal",rho=0.5,m=80,n.x=0,n.y=10,loc.2=0.3) 
    D8.3=ranktests(10000,"two.sided","normal",rho=0.8,m=80,n.x=0,n.y=10,loc.2=0.3)
    D8.4=ranktests(10000,"two.sided","normal",rho=0,m=80,n.x=0,n.y=40,loc.2=0.3) 
    D8.5=ranktests(10000,"two.sided","normal",rho=0.5,m=80,n.x=0,n.y=40,loc.2=0.3) 
    D8.6=ranktests(10000,"two.sided","normal",rho=0.8,m=80,n.x=0,n.y=40,loc.2=0.3)    
    D8.7=ranktests(10000,"two.sided","normal",rho=0,m=80,n.x=0,n.y=80,loc.2=0.3) 
    D8.8=ranktests(10000,"two.sided","normal",rho=0.5,m=80,n.x=0,n.y=80,loc.2=0.3) 
    D8.9=ranktests(10000,"two.sided","normal",rho=0.8,m=80,n.x=0,n.y=80,loc.2=0.3)  
   
    #normal--m=80--delta=0.6
    D9.1=ranktests(10000,"two.sided","normal",rho=0,m=80,n.x=0,n.y=10,loc.2=0.6) 
    D9.2=ranktests(10000,"two.sided","normal",rho=0.5,m=80,n.x=0,n.y=10,loc.2=0.6) 
    D9.3=ranktests(10000,"two.sided","normal",rho=0.8,m=80,n.x=0,n.y=10,loc.2=0.6)
    D9.4=ranktests(10000,"two.sided","normal",rho=0,m=80,n.x=0,n.y=40,loc.2=0.6) 
    D9.5=ranktests(10000,"two.sided","normal",rho=0.5,m=80,n.x=0,n.y=40,loc.2=0.6) 
    D9.6=ranktests(10000,"two.sided","normal",rho=0.8,m=80,n.x=0,n.y=40,loc.2=0.6)
    D9.7=ranktests(10000,"two.sided","normal",rho=0,m=80,n.x=0,n.y=80,loc.2=0.6) 
    D9.8=ranktests(10000,"two.sided","normal",rho=0.5,m=80,n.x=0,n.y=80,loc.2=0.6) 
    D9.9=ranktests(10000,"two.sided","normal",rho=0.8,m=80,n.x=0,n.y=80,loc.2=0.6)
    
    ##################################################################   
    #lognormal--m=20--delta=0
    E1.1=ranktests(10000,"two.sided","lognormal",rho=0,m=20,n.x=0,n.y=10,loc.2=0) 
    E1.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=20,n.x=0,n.y=10,loc.2=0) 
    E1.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=20,n.x=0,n.y=10,loc.2=0)
    E1.4=ranktests(10000,"two.sided","lognormal",rho=0,m=20,n.x=0,n.y=40,loc.2=0) 
    E1.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=20,n.x=0,n.y=40,loc.2=0) 
    E1.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=20,n.x=0,n.y=40,loc.2=0)
    
    #lognormal--m=20--delta=0.3
    E2.1=ranktests(10000,"two.sided","lognormal",rho=0,m=20,n.x=0,n.y=10,loc.2=0.3) 
    E2.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=20,n.x=0,n.y=10,loc.2=0.3) 
    E2.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=20,n.x=0,n.y=10,loc.2=0.3)
    E2.4=ranktests(10000,"two.sided","lognormal",rho=0,m=20,n.x=0,n.y=40,loc.2=0.3) 
    E2.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=20,n.x=0,n.y=40,loc.2=0.3) 
    E2.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=20,n.x=0,n.y=40,loc.2=0.3)  
    
    #lognormal--m=20--delta=0.6
    E3.1=ranktests(10000,"two.sided","lognormal",rho=0,m=20,n.x=0,n.y=10,loc.2=0.6) 
    E3.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=20,n.x=0,n.y=10,loc.2=0.6) 
    E3.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=20,n.x=0,n.y=10,loc.2=0.6)
    E3.4=ranktests(10000,"two.sided","lognormal",rho=0,m=20,n.x=0,n.y=40,loc.2=0.6) 
    E3.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=20,n.x=0,n.y=40,loc.2=0.6) 
    E3.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=20,n.x=0,n.y=40,loc.2=0.6)
    
    #lognormal--m=40--delta=0
    E4.1=ranktests(10000,"two.sided","lognormal",rho=0,m=40,n.x=0,n.y=10,loc.2=0) 
    E4.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=40,n.x=0,n.y=10,loc.2=0) 
    E4.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=40,n.x=0,n.y=10,loc.2=0)
    E4.4=ranktests(10000,"two.sided","lognormal",rho=0,m=40,n.x=0,n.y=40,loc.2=0) 
    E4.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=40,n.x=0,n.y=40,loc.2=0) 
    E4.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=40,n.x=0,n.y=40,loc.2=0)
       
    #lognormal--m=40--delta=0.3
    E5.1=ranktests(10000,"two.sided","lognormal",rho=0,m=40,n.x=0,n.y=10,loc.2=0.3) 
    E5.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=40,n.x=0,n.y=10,loc.2=0.3) 
    E5.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=40,n.x=0,n.y=10,loc.2=0.3)
    E5.4=ranktests(10000,"two.sided","lognormal",rho=0,m=40,n.x=0,n.y=40,loc.2=0.3) 
    E5.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=40,n.x=0,n.y=40,loc.2=0.3) 
    E5.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=40,n.x=0,n.y=40,loc.2=0.3)    
    
    #lognormal--m=40--delta=0.6
    E6.1=ranktests(10000,"two.sided","lognormal",rho=0,m=40,n.x=0,n.y=10,loc.2=0.6) 
    E6.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=40,n.x=0,n.y=10,loc.2=0.6) 
    E6.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=40,n.x=0,n.y=10,loc.2=0.6)
    E6.4=ranktests(10000,"two.sided","lognormal",rho=0,m=40,n.x=0,n.y=40,loc.2=0.6) 
    E6.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=40,n.x=0,n.y=40,loc.2=0.6) 
    E6.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=40,n.x=0,n.y=40,loc.2=0.6)

    #lognormal--m=80--delta=0
    E7.1=ranktests(10000,"two.sided","lognormal",rho=0,m=80,n.x=0,n.y=10,loc.2=0) 
    E7.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=80,n.x=0,n.y=10,loc.2=0) 
    E7.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=80,n.x=0,n.y=10,loc.2=0)
    E7.4=ranktests(10000,"two.sided","lognormal",rho=0,m=80,n.x=0,n.y=40,loc.2=0) 
    E7.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=80,n.x=0,n.y=40,loc.2=0) 
    E7.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=80,n.x=0,n.y=40,loc.2=0)
    E7.7=ranktests(10000,"two.sided","lognormal",rho=0,m=80,n.x=0,n.y=80,loc.2=0) 
    E7.8=ranktests(10000,"two.sided","lognormal",rho=0.5,m=80,n.x=0,n.y=80,loc.2=0) 
    E7.9=ranktests(10000,"two.sided","lognormal",rho=0.8,m=80,n.x=0,n.y=80,loc.2=0)
       
    #lognormal--m=80--delta=0.3
    E8.1=ranktests(10000,"two.sided","lognormal",rho=0,m=80,n.x=0,n.y=10,loc.2=0.3) 
    E8.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=80,n.x=0,n.y=10,loc.2=0.3) 
    E8.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=80,n.x=0,n.y=10,loc.2=0.3)
    E8.4=ranktests(10000,"two.sided","lognormal",rho=0,m=80,n.x=0,n.y=40,loc.2=0.3) 
    E8.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=80,n.x=0,n.y=40,loc.2=0.3) 
    E8.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=80,n.x=0,n.y=40,loc.2=0.3)    
    E8.7=ranktests(10000,"two.sided","lognormal",rho=0,m=80,n.x=0,n.y=80,loc.2=0.3) 
    E8.8=ranktests(10000,"two.sided","lognormal",rho=0.5,m=80,n.x=0,n.y=80,loc.2=0.3) 
    E8.9=ranktests(10000,"two.sided","lognormal",rho=0.8,m=80,n.x=0,n.y=80,loc.2=0.3)  
   
    #lognormal--m=80--delta=0.6
    E9.1=ranktests(10000,"two.sided","lognormal",rho=0,m=80,n.x=0,n.y=10,loc.2=0.6) 
    E9.2=ranktests(10000,"two.sided","lognormal",rho=0.5,m=80,n.x=0,n.y=10,loc.2=0.6) 
    E9.3=ranktests(10000,"two.sided","lognormal",rho=0.8,m=80,n.x=0,n.y=10,loc.2=0.6)
    E9.4=ranktests(10000,"two.sided","lognormal",rho=0,m=80,n.x=0,n.y=40,loc.2=0.6) 
    E9.5=ranktests(10000,"two.sided","lognormal",rho=0.5,m=80,n.x=0,n.y=40,loc.2=0.6) 
    E9.6=ranktests(10000,"two.sided","lognormal",rho=0.8,m=80,n.x=0,n.y=40,loc.2=0.6)
    E9.7=ranktests(10000,"two.sided","lognormal",rho=0,m=80,n.x=0,n.y=80,loc.2=0.6) 
    E9.8=ranktests(10000,"two.sided","lognormal",rho=0.5,m=80,n.x=0,n.y=80,loc.2=0.6) 
    E9.9=ranktests(10000,"two.sided","lognormal",rho=0.8,m=80,n.x=0,n.y=80,loc.2=0.6)

    ##################################################################   
    #student--m=20--delta=0(df=3)
    F1.1=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=0,n.y=10,loc.2=0) 
    F1.2=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=0,n.y=10,loc.2=0) 
    F1.3=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=0,n.y=10,loc.2=0)
    F1.4=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=0,n.y=40,loc.2=0) 
    F1.5=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=0,n.y=40,loc.2=0) 
    F1.6=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=0,n.y=40,loc.2=0)
    
    #student--m=20--delta=0.3(df=3)
    F2.1=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=0,n.y=10,loc.2=0.3) 
    F2.2=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=0,n.y=10,loc.2=0.3) 
    F2.3=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=0,n.y=10,loc.2=0.3)
    F2.4=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=0,n.y=40,loc.2=0.3) 
    F2.5=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=0,n.y=40,loc.2=0.3) 
    F2.6=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=0,n.y=40,loc.2=0.3)  
    
    #student--m=20--delta=0.6(df=3)
    F3.1=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=0,n.y=10,loc.2=0.6) 
    F3.2=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=0,n.y=10,loc.2=0.6) 
    F3.3=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=0,n.y=10,loc.2=0.6)
    F3.4=ranktests(10000,"two.sided","student",rho=0,m=20,n.x=0,n.y=40,loc.2=0.6) 
    F3.5=ranktests(10000,"two.sided","student",rho=0.5,m=20,n.x=0,n.y=40,loc.2=0.6) 
    F3.6=ranktests(10000,"two.sided","student",rho=0.8,m=20,n.x=0,n.y=40,loc.2=0.6)
    
    #student--m=40--delta=0(df=3)
    F4.1=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=0,n.y=10,loc.2=0) 
    F4.2=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=0,n.y=10,loc.2=0) 
    F4.3=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=0,n.y=10,loc.2=0)
    F4.4=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=0,n.y=40,loc.2=0) 
    F4.5=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=0,n.y=40,loc.2=0) 
    F4.6=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=0,n.y=40,loc.2=0)
       
    #student--m=40--delta=0.3(df=3)
    F5.1=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=0,n.y=10,loc.2=0.3) 
    F5.2=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=0,n.y=10,loc.2=0.3) 
    F5.3=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=0,n.y=10,loc.2=0.3)
    F5.4=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=0,n.y=40,loc.2=0.3) 
    F5.5=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=0,n.y=40,loc.2=0.3) 
    F5.6=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=0,n.y=40,loc.2=0.3)    
    
    #student--m=40--delta=0.6(df=3)
    F6.1=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=0,n.y=10,loc.2=0.6) 
    F6.2=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=0,n.y=10,loc.2=0.6) 
    F6.3=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=0,n.y=10,loc.2=0.6)
    F6.4=ranktests(10000,"two.sided","student",rho=0,m=40,n.x=0,n.y=40,loc.2=0.6) 
    F6.5=ranktests(10000,"two.sided","student",rho=0.5,m=40,n.x=0,n.y=40,loc.2=0.6) 
    F6.6=ranktests(10000,"two.sided","student",rho=0.8,m=40,n.x=0,n.y=40,loc.2=0.6)
    
    #student--m=80--delta=0(df=3)
    F7.1=ranktests(10000,"two.sided","student",rho=0,m=80,n.x=0,n.y=10,loc.2=0) 
    F7.2=ranktests(10000,"two.sided","student",rho=0.5,m=80,n.x=0,n.y=10,loc.2=0) 
    F7.3=ranktests(10000,"two.sided","student",rho=0.8,m=80,n.x=0,n.y=10,loc.2=0)
    F7.4=ranktests(10000,"two.sided","student",rho=0,m=80,n.x=0,n.y=40,loc.2=0) 
    F7.5=ranktests(10000,"two.sided","student",rho=0.5,m=80,n.x=0,n.y=40,loc.2=0) 
    F7.6=ranktests(10000,"two.sided","student",rho=0.8,m=80,n.x=0,n.y=40,loc.2=0)
    F7.7=ranktests(10000,"two.sided","student",rho=0,m=80,n.x=0,n.y=80,loc.2=0) 
    F7.8=ranktests(10000,"two.sided","student",rho=0.5,m=80,n.x=0,n.y=80,loc.2=0) 
    F7.9=ranktests(10000,"two.sided","student",rho=0.8,m=80,n.x=0,n.y=80,loc.2=0)
       
    #student--m=80--delta=0.3(df=3)
    F8.1=ranktests(10000,"two.sided","student",rho=0,m=80,n.x=0,n.y=10,loc.2=0.3) 
    F8.2=ranktests(10000,"two.sided","student",rho=0.5,m=80,n.x=0,n.y=10,loc.2=0.3) 
    F8.3=ranktests(10000,"two.sided","student",rho=0.8,m=80,n.x=0,n.y=10,loc.2=0.3)
    F8.4=ranktests(10000,"two.sided","student",rho=0,m=80,n.x=0,n.y=40,loc.2=0.3) 
    F8.5=ranktests(10000,"two.sided","student",rho=0.5,m=80,n.x=0,n.y=40,loc.2=0.3) 
    F8.6=ranktests(10000,"two.sided","student",rho=0.8,m=80,n.x=0,n.y=40,loc.2=0.3)    
    F8.7=ranktests(10000,"two.sided","student",rho=0,m=80,n.x=0,n.y=80,loc.2=0.3) 
    F8.8=ranktests(10000,"two.sided","student",rho=0.5,m=80,n.x=0,n.y=80,loc.2=0.3) 
    F8.9=ranktests(10000,"two.sided","student",rho=0.8,m=80,n.x=0,n.y=80,loc.2=0.3)  
   
    #student--m=80--delta=0.6(df=3)
    F9.1=ranktests(10000,"two.sided","student",rho=0,m=80,n.x=0,n.y=10,loc.2=0.6) 
    F9.2=ranktests(10000,"two.sided","student",rho=0.5,m=80,n.x=0,n.y=10,loc.2=0.6) 
    F9.3=ranktests(10000,"two.sided","student",rho=0.8,m=80,n.x=0,n.y=10,loc.2=0.6)
    F9.4=ranktests(10000,"two.sided","student",rho=0,m=80,n.x=0,n.y=40,loc.2=0.6) 
    F9.5=ranktests(10000,"two.sided","student",rho=0.5,m=80,n.x=0,n.y=40,loc.2=0.6) 
    F9.6=ranktests(10000,"two.sided","student",rho=0.8,m=80,n.x=0,n.y=40,loc.2=0.6)
    F9.7=ranktests(10000,"two.sided","student",rho=0,m=80,n.x=0,n.y=80,loc.2=0.6) 
    F9.8=ranktests(10000,"two.sided","student",rho=0.5,m=80,n.x=0,n.y=80,loc.2=0.6) 
    F9.9=ranktests(10000,"two.sided","student",rho=0.8,m=80,n.x=0,n.y=80,loc.2=0.6)
    
    #big sample
    ##################################################################
    #logistic--m=20--delta=0
    H1.1=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=0,n.y=10,loc.2=0) 
    H1.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=0,n.y=10,loc.2=0) 
    H1.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=0,n.y=10,loc.2=0)
    H1.4=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=0,n.y=40,loc.2=0) 
    H1.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=0,n.y=40,loc.2=0) 
    H1.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=0,n.y=40,loc.2=0)
    
    #logistic--m=20--delta=0.3
    H2.1=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=0,n.y=10,loc.2=0.3) 
    H2.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=0,n.y=10,loc.2=0.3) 
    H2.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=0,n.y=10,loc.2=0.3)
    H2.4=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=0,n.y=40,loc.2=0.3) 
    H2.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=0,n.y=40,loc.2=0.3) 
    H2.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=0,n.y=40,loc.2=0.3)  
    
    #logistic--m=20--delta=0.6
    H3.1=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=0,n.y=10,loc.2=0.6) 
    H3.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=0,n.y=10,loc.2=0.6) 
    H3.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=0,n.y=10,loc.2=0.6)
    H3.4=ranktests(10000,"two.sided","logistic",rho=0,m=20,n.x=0,n.y=40,loc.2=0.6) 
    H3.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=20,n.x=0,n.y=40,loc.2=0.6) 
    H3.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=20,n.x=0,n.y=40,loc.2=0.6)
    
    #logistic--m=40--delta=0
    H4.1=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=0,n.y=10,loc.2=0) 
    H4.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=0,n.y=10,loc.2=0) 
    H4.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=0,n.y=10,loc.2=0)
    H4.4=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=0,n.y=40,loc.2=0) 
    H4.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=0,n.y=40,loc.2=0) 
    H4.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=0,n.y=40,loc.2=0)
       
    #logistic--m=40--delta=0.3
    H5.1=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=0,n.y=10,loc.2=0.3) 
    H5.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=0,n.y=10,loc.2=0.3) 
    H5.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=0,n.y=10,loc.2=0.3)
    H5.4=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=0,n.y=40,loc.2=0.3) 
    H5.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=0,n.y=40,loc.2=0.3) 
    H5.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=0,n.y=40,loc.2=0.3)    
    
    #logistic--m=40--delta=0.6
    H6.1=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=0,n.y=10,loc.2=0.6) 
    H6.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=0,n.y=10,loc.2=0.6) 
    H6.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=0,n.y=10,loc.2=0.6)
    H6.4=ranktests(10000,"two.sided","logistic",rho=0,m=40,n.x=0,n.y=40,loc.2=0.6) 
    H6.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=40,n.x=0,n.y=40,loc.2=0.6) 
    H6.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=40,n.x=0,n.y=40,loc.2=0.6)
    
    #logistic--m=80--delta=0
    H7.1=ranktests(10000,"two.sided","logistic",rho=0,m=80,n.x=0,n.y=10,loc.2=0) 
    H7.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=80,n.x=0,n.y=10,loc.2=0) 
    H7.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=80,n.x=0,n.y=10,loc.2=0)
    H7.4=ranktests(10000,"two.sided","logistic",rho=0,m=80,n.x=0,n.y=40,loc.2=0) 
    H7.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=80,n.x=0,n.y=40,loc.2=0) 
    H7.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=80,n.x=0,n.y=40,loc.2=0)
    H7.7=ranktests(10000,"two.sided","logistic",rho=0,m=80,n.x=0,n.y=80,loc.2=0) 
    H7.8=ranktests(10000,"two.sided","logistic",rho=0.5,m=80,n.x=0,n.y=80,loc.2=0) 
    H7.9=ranktests(10000,"two.sided","logistic",rho=0.8,m=80,n.x=0,n.y=80,loc.2=0)
       
    #logistic--m=80--delta=0.3
    H8.1=ranktests(10000,"two.sided","logistic",rho=0,m=80,n.x=0,n.y=10,loc.2=0.3) 
    H8.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=80,n.x=0,n.y=10,loc.2=0.3) 
    H8.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=80,n.x=0,n.y=10,loc.2=0.3)
    H8.4=ranktests(10000,"two.sided","logistic",rho=0,m=80,n.x=0,n.y=40,loc.2=0.3) 
    H8.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=80,n.x=0,n.y=40,loc.2=0.3) 
    H8.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=80,n.x=0,n.y=40,loc.2=0.3)    
    H8.7=ranktests(10000,"two.sided","logistic",rho=0,m=80,n.x=0,n.y=80,loc.2=0.3) 
    H8.8=ranktests(10000,"two.sided","logistic",rho=0.5,m=80,n.x=0,n.y=80,loc.2=0.3) 
    H8.9=ranktests(10000,"two.sided","logistic",rho=0.8,m=80,n.x=0,n.y=80,loc.2=0.3)  
   
    #logistic--m=80--delta=0.6
    H9.1=ranktests(10000,"two.sided","logistic",rho=0,m=80,n.x=0,n.y=10,loc.2=0.6) 
    H9.2=ranktests(10000,"two.sided","logistic",rho=0.5,m=80,n.x=0,n.y=10,loc.2=0.6) 
    H9.3=ranktests(10000,"two.sided","logistic",rho=0.8,m=80,n.x=0,n.y=10,loc.2=0.6)
    H9.4=ranktests(10000,"two.sided","logistic",rho=0,m=80,n.x=0,n.y=40,loc.2=0.6) 
    H9.5=ranktests(10000,"two.sided","logistic",rho=0.5,m=80,n.x=0,n.y=40,loc.2=0.6) 
    H9.6=ranktests(10000,"two.sided","logistic",rho=0.8,m=80,n.x=0,n.y=40,loc.2=0.6)
    H9.7=ranktests(10000,"two.sided","logistic",rho=0,m=80,n.x=0,n.y=80,loc.2=0.6) 
    H9.8=ranktests(10000,"two.sided","logistic",rho=0.5,m=80,n.x=0,n.y=80,loc.2=0.6) 
    H9.9=ranktests(10000,"two.sided","logistic",rho=0.8,m=80,n.x=0,n.y=80,loc.2=0.6)

