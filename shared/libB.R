makeDataFrameML3=function(filename)
{
	indat<-read.csv(file=filename,header=TRUE)
	temp=data.frame(indat$Site,indat$Genderfactor,indat$CredCond,indat$mcdv1)
	badMat=as.matrix(is.na(temp))
	good=apply(badMat,1,mean)==0
	dat=temp[good,]
	colnames(dat)=c('lab','gender','cred','Y')
	dat$lab=as.numeric(dat$lab)
	dat$gender=as.numeric(dat$gender) #female=1, Male=2
	dat$cred=as.numeric(dat$cred) # Credentials=1, NoCredentials=2
	return(dat)
}

freqEst=function(dat)
{
	J=max(dat$lab)
	site=tapply(dat$Y[dat$cred==1],dat$lab[dat$cred==1],mean)
	ciSite= t(sapply(1:J, function(x) t.test(dat$Y[dat$lab==x & dat$cred==1])$conf.int))
	m=tapply(dat$Y,list(dat$lab,dat$gender),mean)
	eG=m[,2]-m[,1]  #men-women
	ciG <- t(-1*sapply(1:J, function(x) t.test(dat$Y[dat$lab==x] ~ dat$gender[dat$lab==x],var.equal=T)$conf.int))
	dm=subset(dat,gender==2)
	m=tapply(dm$Y,list(dm$lab,dm$cred),mean)
	eM=m[,1]-m[,2]  #Credential-noCredential
	ciM <- t(sapply(1:J, function(x) t.test(dm$Y[dm$lab==x] ~ dm$cred[dm$lab==x],var.equal=T)$conf.int))
	dw=subset(dat,gender==1)
	m=tapply(dw$Y,list(dw$lab,dw$cred),mean)
	eW=m[,1]-m[,2]  #Credential-noCredential
	ciW <- t(sapply(1:J, function(x) t.test(dw$Y[dw$lab==x] ~ dw$cred[dw$lab==x],var.equal=T)$conf.int))
	contrast=rep(c("site","gender","credMen","credWomen"),each=J)
	effect=c(site,eG,eM,eW)
	ci=rbind(ciSite,ciG,ciM,ciW)
	o=c(order(site),order(eG),order(eM),order(eW))
	out=data.frame(contrast,effect,ci,o)
	colnames(out)=c('contrast','center','upper','lower','order')
	return(out)	
}

doChain=function(dat,rScaleIndv=.3,rScaleMean=.5){
	library('BayesFactor')
	I=max(dat$lab)
	N=length(dat$lab)
	nu=1:I
	theta=(I+1):(2*I)
	alpha=(2*I+1):(3*I)
	beta=(3*I+1):(4*I)
	mu.theta=4*I+1
	mu.alpha=4*I+2
	mu.beta=4*I+3
	Xf=matrix(0,ncol=4*I+3,nrow=N)
	for (i in 1:N){
  		Xf[i,nu[dat$lab[i]]]=1
	 	Xf[i,theta[dat$lab[i]]]=dat$gender[i]-1.5
  		Xf[i,alpha[dat$lab[i]]]=(dat$cred[i]-1)*(dat$gender[i]-1)
  		Xf[i,beta[dat$lab[i]]]=(dat$cred[i]-1)*(1-(dat$gender[i]-1))
  		Xf[i,mu.theta]=dat$gender[i]-1.5
  		Xf[i,mu.alpha]=(dat$cred[i]-1)*(dat$gender[i]-1)
  		Xf[i,mu.beta]=(dat$cred[i]-1)*(1-(dat$gender[i]-1))}
	gMap=rep(0:6,c(I,I,I,I,1,1,1))
	samp=nWayAOV(dat$Y,Xf,gMap=gMap,rscale=rep(c(rScaleIndv,rScaleMean),c(4,3)),posterior=T)
	return(samp)
}

bayesEst=function(samp){
	nVar=dim(samp)[2]
	I=((nVar-8)/4)-1
	alpha=(2*I+1):(3*I)
	beta=(3*I+1):(4*I)
	theta=(I+1):(2*I)
	mu.alpha=4*I+2
	mu.beta=4*I+3
	mu.theta=4*I+1
	mu=1
	nu=1:I
	mcmc=samp[,nu+1]+samp[,mu]
	estSite=apply(mcmc,2,mean)
	qSite=apply(mcmc,2,quantile,p=c(.975,.025))	
	mcmc=samp[,theta+1]+samp[,mu.theta+1]
	estG=apply(mcmc,2,mean)
	qG=apply(mcmc,2,quantile,p=c(.975,.025))
	mcmc=samp[,alpha+1]+samp[,mu.alpha+1]
	estM=apply(mcmc,2,mean)
	qM=apply(mcmc,2,quantile,p=c(.975,.025))
	mcmc=samp[,beta+1]+samp[,mu.beta+1]
	estW=apply(mcmc,2,mean)
	qW=apply(mcmc,2,quantile,p=c(.975,.025))
	contrast=rep(c("site","gender","credMen","credWomen"),each=I)
	effect=c(estSite,estG,-estM,-estW)
	ci=rbind(t(qSite),t(qG),t(-qM),t(-qW))
	o=c(order(estSite),order(estG),order(-estM),order(-estW))
	out=data.frame(contrast,effect,ci,o)
	colnames(out)=c('contrast','center','upper','lower','order')
	return(out)	
}

plotter=function(f,b,axisText,range=1.5,main){
	J=length(f$center)
	o=f$order
	plot(f$center[o],1:J, xlim=c(-range,range), axes=F, xlab="", typ='n', ylab="Site", ylim=c(0,J))
#	abline(v=mean(f$center),lty=2,col="grey30")
	arrows(f$lower[o],1:J,f$upper[o],1:J,code=3,angle=90,length=.05)
	abline(v=0,lty=2)
	points(f$center[o],1:J,pch=21,bg='white',cex=1.1)
	axis(2,at=c(1,J))
#	mtext(side=3,line=-2,"A.",cex=1.4,adj=.05)
	mtext(side=3,line=-3,"Sample",cex=1.2,adj=.05)
	mtext(side=3,line=1,main,cex=1.4,adj=.5)
	plot(b$center[o],1:J,xlim=c(-range,range),axes=F,xlab="",typ='n',ylab="Site",ylim=c(0,J))
#	abline(v=mean(b$center),lty=2,col="grey30")
	arrows(b$lower[o],1:J,b$upper[o],1:J,code=3,angle=90,length=.05)
	abline(v=0,lty=2)
	points(f$center[o],1:J,pch=4,cex=1.1)
	points(b$center[o],1:J,pch=21,bg='black',cex=1.1)
	axis(1,at=seq(-range,range,1),label=c(axisText[1], seq(-range+1,range-1,1),axisText[2]))
	axis(2,at=c(1,J))
#	mtext(side=3,line=-2,"B.",cex=1.4,adj=.05)
	mtext(side=3,line=-3,"Hierarchical",cex=1.2,adj=.05)
}

bfTable=function(dat,rScaleMean,rScaleIndv){

	I=max(dat$lab)
	N=length(dat$lab)
	Xnu=matrix(0,ncol=I,nrow=N)
	for (i in 1:N) Xnu[i,dat$lab[i]]=1  #0/1 coded
	Xgen=matrix(0,ncol=I,nrow=N)
	for (i in 1:N) Xgen[i,dat$lab[i]]=dat$gender[i]-1.5 #-.5/.5 coded
	XcredM=matrix(0,ncol=I,nrow=N)
	for (i in 1:N) XcredM[i,dat$lab[i]]=(dat$cred[i]-1)*(dat$gender[i]-1)  #0/1 coded
	XcredW=matrix(0,ncol=I,nrow=N)
	for (i in 1:N) XcredW[i,dat$lab[i]]=(dat$cred[i]-1)*(1-(dat$gender[i]-1))  #0/1 coded
	Xgen0=dat$gender-1.5
	XcredM0=(dat$cred-1)*(dat$gender-1)
	XcredW0=(dat$cred-1)*(1-(dat$gender-1))
	Xcred=XcredM+XcredW
	Xcred0=XcredM0+XcredW0

	#prior probability
	library('MCMCpack')
	M=100000
	mn=rcauchy(M,0,rScaleMean)
	s2=rinvgamma(M,.5,.5*rScaleIndv^2)
	count=0
	for (m in 1:M){
		prior=rnorm(I,mn[m],sqrt(s2[m]))
		if (mean(prior>0)==1) count=count+1}
	prior.pos=count/M

labels=c(
	'Common Site + Common Gender + Common Credential',
	'Common Site + Common Gender + Common Men Credential',
	'Common Site + Common Gender + Common 2 Credentials',
	'Common Site + Common Gender' ,
	'Common Site + Commen Credential',
	'Positive Site + Common Gender + Common Credential',
	'Unconstrained Site + Common Gender + Common Credential',
	'Common Site + Positive Gender + Common Credential',
	'Common Site + Common Gender + Positive Credential',
	'Common Site + Common Gender + Unconstrained Credential',
	'Comon Site + Unconstrained Credential + Common Credential'
	)

bf=1:length(labels)
bf[1]=nWayAOV(dat$Y,cbind(Xgen0,Xcred0),gMap=rep(0:1,c(1,1)),rscale=rep(rScaleMean,2))$bf
bf[2]=nWayAOV(dat$Y,cbind(Xgen0,XcredM0),gMap=rep(0:1,c(1,1)),rscale=rep(rScaleMean,2))$bf
bf[3]=nWayAOV(dat$Y,cbind(Xgen0,XcredM0,XcredW0),gMap=rep(0:2,c(1,1,1)),rscale=rep(rScaleMean,3))$bf
bf[4]=nWayAOV(dat$Y,cbind(Xgen0),gMap=rep(0,1),rscale=rep(rScaleMean,1))$bf
bf[5]=nWayAOV(dat$Y,cbind(Xcred0),gMap=rep(0,1),rscale=rep(rScaleMean,1))$bf

#site intercept effects
bf[7]=nWayAOV(dat$Y,cbind(Xnu,Xgen0,Xcred0),gMap=rep(0:2,c(I,1,1)),rscale=rep(c(rScaleIndv,rScaleMean),c(1,2)))$bf
samp=nWayAOV(dat$Y,cbind(Xnu,Xgen0,Xcred0),gMap=rep(0:2,c(I,1,1)),rscale=rep(c(rScaleIndv,rScaleMean),c(1,2)),posterior=T)
post.pos=mean(apply(samp[,1]+samp[,1+(1:I)]>0,1,mean)==1)
bf[6]=bf[7]+log(post.pos)-log(prior.pos)

#gender effects
bf[9]=nWayAOV(dat$Y,cbind(Xgen,Xgen0,Xcred0),gMap=rep(0:2,c(I,1,1)),rscale=rep(c(rScaleIndv,rScaleMean),c(1,2)))$bf
samp=nWayAOV(dat$Y,cbind(Xgen,Xgen0,Xcred0),gMap=rep(0:2,c(I,1,1)),rscale=rep(c(rScaleIndv,rScaleMean),c(1,2)),posterior=T)
post.pos=mean(apply(samp[,1+(I+1)]+samp[,1+(1:I)]>0,1,mean)==1)
bf[8]=bf[9]+log(post.pos)-log(prior.pos)

#credential effects
bf[11]=nWayAOV(dat$Y,cbind(Xcred,Xgen0,Xcred0),gMap=rep(0:2,c(I,1,1)),rscale=rep(c(rScaleIndv,rScaleMean),c(1,2)))$bf
samp=nWayAOV(dat$Y,cbind(Xcred,Xgen0,Xcred0),gMap=rep(0:2,c(I,1,1)),rscale=rep(c(rScaleIndv,rScaleMean),c(1,2)),posterior=T)
post.pos=mean(apply(samp[,1+(I+2)]+samp[,1+(1:I)]>0,1,mean)==1)
bf[10]=bf[11]+log(post.pos)-log(prior.pos)

printBF=paste('1-to-',round(1/exp(bf-bf[1]),1),sep='')
myTab=data.frame(labels,printBF)
colnames(myTab)=c("Model","Bayes Factor")
return(myTab)
}
