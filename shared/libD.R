library('devtools')
library('BayesFactor')
library('MCMCpack')

process=function(dat,site,fac=c("cong","inc"),start=0){
	y=dat$rt
	sub=as.integer(dat$sub)
	cond=dat$cond
	mrt=as.data.frame.table(tapply(y,list(sub,cond),mean))
	colnames(mrt)=c('sub','cond','y')
	mrt$cond=as.factor(mrt$cond)
	mrt$sub=as.integer(mrt$sub)+start
	levels(mrt$cond)=fac
	site=rep(site,length(mrt$y))
	mrt=data.frame(site,mrt)
	return(mrt)
}

makeDatHR=function(){
	
	SourceURL <- "https://raw.githubusercontent.com/PerceptionCognitionLab/data0/master/contexteffects/FlankerStroopSimon/cleaning.R"
	source_url(SourceURL)
	stroop$sub <- stroop$ID
	dat1=process(stroop,site=1,start=0)

	SourceURL <- "https://raw.githubusercontent.com/PerceptionCognitionLab/data0/master/contexteffects/StroopSimonAPP2010/cleaning.R"
	source_url(SourceURL)
	dat.stroop.p1$cond=1-dat.stroop.p1$cond
	dat.stroop.p2$cond=1-dat.stroop.p2$cond
	dat2=process(dat.stroop.p1,2,start=121)
	dat3=process(dat.stroop.p2,3,start=159)
	dat=rbind(dat1,dat2,dat3)
	return(dat)
	}


doBayesFactor=function(dat,rscale)
{
	#figure sigma at 30ms for rscale
	#rscale=c(intercept,effect,intercept.site,effect.site,effect.site.mean)	
	N=length(dat$site)
	I=max(dat$site)
	J=max(dat$sub)
	cond=as.integer(dat$cond=='inc')
	sub=dat$sub
	site=dat$site
	#Y_ijk=mu+alpha_j+x_k*theta_j
	X=matrix(nrow=N,ncol=2*J+2*I+1,0)
	for (r in 1:N){
		X[r,sub[r]]=1
		X[r,J+sub[r]]=cond[r]
	 	X[r,2*J+site[r]]=1
	 	X[r,2*J+I+site[r]]=cond[r]
	 	X[r,2*J+2*I+1]=cond[r]}
	gMap=rep(0:4,c(J,J,I,I,1))
	bfFull=nWayAOV(dat$y,X,gMap=gMap,rscale=rscale)$bf
	samples=nWayAOV(dat$y,X,gMap=gMap,rscale=rscale,posterior=T)

	X1=matrix(nrow=N,ncol=2*J+I+1,0)
	for (r in 1:N){
		X1[r,sub[r]]=1
	 	X1[r,J+sub[r]]=cond[r]
	 	X1[r,2*J+site[r]]=1
	 	X1[r,2*J+I+1]=cond[r]}
	gMap=rep(0:3,c(J,J,I,1))
	bfOne=nWayAOV(dat$y,X1,gMap=gMap,rscale=rscale[1:4])$bf

	X0=matrix(nrow=N,ncol=2*J+I,0)
	for (r in 1:N){
		 X0[r,sub[r]]=1
		 X0[r,J+sub[r]]=cond[r]
		 X0[r,2*J+site[r]]=1}
	gMap=rep(0:2,c(J,J,I))
	bfNone=nWayAOV(dat$y,X0,gMap=gMap,rscale=rscale[1:3])$bf


	pm=apply(samples,2,mean)
	effect=samples[,(2*J+I+2):(2*J+2*I+1)]+samples[,2*J+2*I+2]
	post.pos=mean(apply(effect>0,1,mean)==1)
	Mprior=10000
	gm=rinvgamma(Mprior,.5,.5*.2^2)
	m=rnorm(Mprior,0,sqrt(gm))
	g=rinvgamma(Mprior,.5,.5*.2^2)
	a1=1:Mprior
	for (i in 1:Mprior)
	{
		a=rnorm(I,m[i],sqrt(g[i]))
		a1[i]=(mean(a>0)==1)
	}
	prior.pos=mean(a1)
	post.pos/prior.pos
	bf=c(exp(bfFull-bfNone),exp(bfOne-bfNone),exp(bfFull-bfNone)*post.pos/prior.pos)
        names(bf)=c("F0","10","P0")
	return(list("bf"=bf,"effect"=effect))	
}


plotter=function(samp.eff,effect){
	m=apply(effect,2,mean)
	I=length(m)
	o=order(m)
    q=apply(effect,2,quantile,p=c(.025,.975))
    plot(m[o],1:I,xlim=c(min(q)-.01,max(q)+.01),typ='n',axes=F,ylab="Study",xlab="Stroop Effect (sec)",ylim=c(.6,I+.4))
    arrows(q[1,o],1:3,q[2,o],1:3,code=3,angle=90,length=.1)
    axis(1)
    axis(2,at=1:I,las=1,label=NA)
    points(m[o],1:I,pch=21,bg='black',cex=1.1)
    points(samp.eff[o],1:I,pch=4,cex=1.1)}


