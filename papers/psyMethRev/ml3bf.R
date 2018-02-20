#pdf('ml3bf.pdf')
source('../../shared/libB.R')
library('BayesFactor')
rScaleIndv=.3
rScaleMean=.5
dat=makeDataFrameML3('../../shared/ML3AllSites.csv')
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


