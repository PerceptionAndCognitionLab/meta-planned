
makeDataFrame2=function(){
	library(RCurl)
	x <- getURL("https://dl.dropboxusercontent.com/s/zsxa2a7vt5dy4vy/music_big5_anon.csv?dl=0", ssl.verifypeer = FALSE)
	indat <- read.csv(text = x)
	N=length(indat$site_name)
	subNum=1:N
	lab=indat$site_name
	site=as.numeric(rep(lab,5))
	sub=rep(subNum,5)
	char=rep(1:5,each=N)
	isAgree=(char==1)
	isCons=(char==2)
	isExtra=(char==3)
	isNeuro=(char==4)
	isOpen=(char==5)	
	val=c(indat$agreeb,indat$cons,indat$extra,indat$neuro,indat$open)
	good=!is.na(val)
	dat=data.frame(site,sub,char,isAgree,isCons,isExtra,isNeuro,isOpen,val)
	return(dat[good,])
}


plotter=function(f,b,...){
		o=order(b[,1])
        m=mean(b[,1])
        xl=c(m-.5,m+.5)
		I=length(f)
        plot(b[o,1],1:I,xlim=xl,axes=F,xlab="",typ='n',ylab="",ylim=c(0,I),...)
        abline(v=m,lty=2,col='grey30')
        arrows(b[o,2],1:I,b[o,3],1:I,code=3,angle=90,length=.05,...)
        par(xpd=T)
        points(f[o],1:I,pch=4,cex=.9,...)
        par(xpd=F)
        points(b[o,1],1:I,pch=19,cex=.9,...)
        text(m,.5,paste("M=",round(m,1),sep=''),adj=-.1,...)
        axis(1,at=m+seq(-.5,.5,.5),label=seq(-.5,.5,.5),...)
}


bf=function(dat,rM,rSD)
{
	N=length(dat$site)
	I=max(dat$site)
	Xsite=Xa=Xc=Xe=Xn=Xo=matrix(nrow=N,ncol=I,0)
	Xsite[cbind(1:N,dat$site)]=1
	Xa[cbind(1:N,dat$site)]=as.numeric(dat$isAgree)
	Xc[cbind(1:N,dat$site)]=as.numeric(dat$isCons)
	Xe[cbind(1:N,dat$site)]=as.numeric(dat$isExtra)
	Xn[cbind(1:N,dat$site)]=as.numeric(dat$isNeuro)
	Xo[cbind(1:N,dat$site)]=as.numeric(dat$isOpen)
	XaM=as.numeric(dat$isAgree)
	XcM=as.numeric(dat$isCons)
	XeM=as.numeric(dat$isExtra)
	XnM=as.numeric(dat$isNeuro)
	XoM=as.numeric(dat$isOpen)

	#Full
	X=cbind(Xsite,Xa,Xc,Xe,Xn,Xo,XaM,XcM,XeM,XnM,XoM)
	g=rep(0:1,c(6*I,5))
	samp=nWayAOV(dat$val,X,gMap=g,rscale=c(rSD,rM),posterior=T)
	bfF=nWayAOV(dat$val,X,gMap=g,rscale=c(rSD,rM))$bf
	
	#Random Intercept
	X=cbind(Xsite,XaM,XcM,XeM,XnM,XoM)
	g=rep(0:1,c(I,5))
	bf1=nWayAOV(dat$val,X,gMap=g,rscale=c(rSD,rM))$bf
	
	#Fixed Intercept
	X=cbind(XaM,XcM,XeM,XnM,XoM)
	g=rep(0,5)
	bf0=nWayAOV(dat$val,X,gMap=g,rscale=c(rM))$bf
	bf=c(bfF,bf1,bf0)
	names(bf)=c('bfF','bf1','bf0')
		
	return(list('samp'=samp,'bf'=bf))
}	

