# Library for Wagenmakers AND ML3 datasets.


###
extractW=function(dat,i){
#Extracts info from a single wagenmakers csv file
dat[is.na(dat)] <- -1 
comprehension=dat[,16]
aware=dat[,17]
madeFace=dat[,8]>2
goodRating=dat[,14] %in% 0:10
keep=(comprehension==1 & aware==0 & madeFace & goodRating)

cond=dat[keep,3]   # 0= POUT an 1= SMILE
ratings=as.matrix(dat[keep,11:14])
class(ratings)="numeric"	
madeSpecificFace=dat[keep,4:7]==1
ratings=ifelse(madeSpecificFace,ratings,NA_character_)
colnames(ratings)=1:4
class(ratings)="numeric"	
ave=apply(ratings,1,mean,na.rm=T)
out=data.frame(rep(i,sum(keep)),cond,ratings,ave)
colnames(out)=c('lab','cond','c1','c2','c3','c4','Y')
return(out)
}

makeDataFrameW=function(fileroot) {
#Merges all .csv files into a data frame	
	fname=paste(fileroot,c('Albohn', 'Allard', 'Benning', 'Bulnes', 'Capaldi', 'Chasten', 'Holmes', 'Koch', 'Korb', 'Lynott', 'Oosterwijk', 'Ozdogru', 'Pacheco-Unguetti', 'Talarico', 'Wagenmakers', 'Wayand', 'Zeelenberg'),"_Data.csv",sep='')
	J=length(fname)	
	j=1
	indat=read.csv(fname[j],header = TRUE, stringsAsFactors = FALSE)
	dat=extractW(indat,j)
	for (j in 2:J){
		indat=read.csv(fname[j],header = TRUE, stringsAsFactors = FALSE)
		datTemp=extractW(indat,j)
		dat=rbind(dat,datTemp)}	
	return(dat)	
}


freqEst=function(dat)
{
	J=max(dat$lab)
	m=tapply(dat$Y,list(dat$lab,dat$cond),mean)
	effect=m[,2]-m[,1]  #Smile - Pout
	ci <- t(-1*sapply(1:J, function(x) t.test(dat$Y[dat$lab==x] ~ dat$cond[dat$lab==x],var.equal=T)		$conf.int))
	out=data.frame(effect,ci,order(effect))
	colnames(out)=c('center','upper','lower','order')
	return(out)	
}

bayesEst=function(dat,rScaleIndv=.3,rScaleMean=.5)
{
	library('BayesFactor')
	N=length(dat$lab)
	J=max(dat$lab)
	cond=(dat$cond-.5)
	alpha=1:J
	beta=(J+1):(2*J)
	mu=2*J+1
	X=matrix(nrow=N,ncol=2*J+1,0)
	for (i in 1:N)
	{
		X[i,alpha[dat$lab[i]]]=1
		X[i,beta[dat$lab[i]]]=cond[i]
		X[i,mu]=cond[i]
	}
	gMap=rep(0:2,c(J,J,1))
	samples=nWayAOV(dat$Y,X,gMap,rscale=c(rScaleIndv,rScaleIndv,rScaleMean),posterior=T)
	mcmc=samples[,beta+1]+samples[,mu+1]
	est=apply(mcmc,2,mean)
	q=apply(mcmc,2,quantile,p=c(.975,.025))
	out=data.frame(est,t(q),order(est))
	colnames(out)=c('center','upper','lower','order')
	return(out)	
}


bayesBF=function(dat,rScaleMean,rScaleIndv)
{
	library('BayesFactor')
	library('MCMCpack')
	N=length(dat$lab)
	J=max(dat$lab)
	cond=(dat$cond-.5)
	alpha=1:J
	beta=(J+1):(2*J)
	mu=2*J+1
	X=matrix(nrow=N,ncol=2*J+1,0)
	for (i in 1:N)
	{
		X[i,alpha[dat$lab[i]]]=1
		X[i,beta[dat$lab[i]]]=cond[i]
		X[i,mu]=cond[i]
	}
	gMap=rep(0:2,c(J,J,1))
	samples=nWayAOV(dat$Y,X,gMap,rscale=c(rScaleIndv,rScaleIndv,rScaleMean),posterior=T)
	bfFull=nWayAOV(dat$Y,X,gMap,rscale=c(rScaleIndv,rScaleIndv,rScaleMean),posterior=F)$bf
	bfNull=nWayAOV(dat$Y,X[,1:J],gMap[1:J],rscale=c(rScaleIndv),posterior=F)$bf
	bfOne=nWayAOV(dat$Y,X[,c(1:J,(2*J+1))],gMap=rep(0:1,c(J,1),),rscale=c(rScaleIndv,rScaleMean),posterior=F)$bf
	
	effect=samples[,beta+1]+samples[,mu+1]
	post.pos=mean(apply(effect>0,1,mean)==1)
	
	Mprior=10000
	gm=rinvgamma(Mprior,.5,.5*rScaleMean^2)
	m.1=rnorm(Mprior,0,sqrt(gm))
	g=rinvgamma(Mprior,.5,.5*rScaleIndv^2)
	a1=1:Mprior
	for (m in 1:Mprior) a1[m]=mean(rnorm(J,m.1[m],sqrt(g[m]))>0)
	prior.pos=mean(a1==1)
	bf=c(exp(bfFull-bfNull),exp(bfOne-bfNull),exp(bfFull-bfNull)*post.pos/prior.pos)
	names(bf)=c("F0","10","P0")
	return(bf)
}


bfPlot=function(bf)
{
library(diagram)
myCol=c('salmon','white','white','white')
names=c('Unconstrained','Positive\n Effects','Common\n Effects','Null')
bfPlot=c(bf[1],bf[3],bf[2],1)
winner=order(bfPlot,decreasing=T)[1]
bfPlot=bfPlot/bfPlot[winner]
other=(1:4)[!((1:4) ==winner)]
M <- matrix(nrow=length(names),byrow=F,ncol=length(names),data=0)
M<-as.data.frame(M)
#for (i in 1:3) M[winner,i]=round(1/bfPlot[other[i]],0)
M[winner,1]<-679
M[winner,2]<-'~1.2e6'
M[winner,3]<-14
G=plotmat(M,c(1,2,1),name=names,curve=0,box.type="circle",box.size=.15,box.prop=.7,box.col=myCol,arr.length=.6,arr.type='triangle',arr.pos=.4,cex=.8)
}



plotter=function(f,b,bf,axisText){
	J=length(f$center)
	o=f$order
	layout(matrix(c(1,2,0,3,3,0),ncol=2,nrow=3),heights=c(1,1,.2),widths=c(1,1))
	par(mar=c(1,4,.5,0),cex=1.1,mgp=c(1,1,0))
	plot(f$center[o],1:J, xlim=c(-1.5,1.5), axes=F, xlab="", typ='n', ylab="Site", ylim=c(0,J))
#	abline(v=mean(f$center),lty=2,col='red')
	arrows(f$lower[o],1:J,f$upper[o],1:J,code=3,angle=90,length=.05)
	abline(v=0,lty=2)
	points(f$center[o],1:J,pch=21,bg='white',cex=1.1)
	axis(2,at=c(1,J))
	mtext(side=3,line=-2,"A.",cex=1.4,adj=.05)
	mtext(side=3,line=-3,"Sample",cex=1.2,adj=.05)
	plot(b$center[o],1:J,xlim=c(-1.5,1.5),axes=F,xlab="",typ='n',ylab="Site",ylim=c(0,J))
#	abline(v=mean(b$center),lty=2,col='red')
	arrows(b$lower[o],1:J,b$upper[o],1:J,code=3,angle=90,length=.05)
	abline(v=0,lty=2)
	points(f$center[o],1:J,pch=4,cex=1.1)
	points(b$center[o],1:J,pch=21,bg='black',cex=1.1)
	axis(1,at=seq(-1.5,1.5,.5),label=c(axisText[1],seq(-1,1,.5),axisText[2]))
	axis(2,at=c(1,J))
	mtext(side=3,line=-2,"B.",cex=1.4,adj=.05)
	mtext(side=3,line=-3,"Hierarchical",cex=1.2,adj=.05)
#	par(mar=c(0,0,.5,0),cex=1.0)
 #   bfPlot(bf)	
#	mtext(side=3,line=-2,"C.",cex=1.4,adj=.95)
#	mtext(side=3,line=-3,"Bayes Factor",cex=1.1,adj=.95)
	
}


#test case
# dat=makeDataFrameW("../../shared/wagenmakersRawData/")
# f=freqEst(dat)
# b=bayesEst(dat)
# bf=bayesBF(dat)
# plotter(f,b,bf,axisText=c("Pout","Smile"))
# dev.off()

