# Rscript --vanilla tools/plot-gradient-scaling.R OUT
out<-commandArgs(TRUE)[1]
t<-read.csv(file.path(out,'timings.csv'),stringsAsFactors=FALSE)
m<-read.csv(file.path(out,'memory.csv'),stringsAsFactors=FALSE)
p<-read.csv(file.path(out,'pairs.csv'),stringsAsFactors=FALSE)
metadata<-unique(t[,c('id','task','n','d','k','g','feedback','observed','direct','events_used')])
p<-merge(p,metadata,by=c('id','task'));p$parameters<-3L*p$g
write.csv(p,file.path(out,'comparisons.csv'),row.names=FALSE)
png(file.path(out,'scaling.png'),width=1800,height=1250,res=150)
par(mfrow=c(2,2),mar=c(4.3,4.3,3,1),oma=c(2.2,0,0,0),las=1)
colours<-c('#B05A37','#256BA0')
draw<-function(x,numeric,analytic,xlab,title) {
  use<-order(x);x<-x[use];numeric<-numeric[use];analytic<-analytic[use]
  plot(x,numeric,type='b',pch=16,col=colours[1],log='y',xaxt='n',xlab=xlab,
    ylab='Value + score (seconds)',ylim=range(c(numeric,analytic)),main=title)
  axis(1,at=x,labels=x);lines(x,analytic,type='b',pch=17,col=colours[2])
  legend('topleft',c('Numeric','Analytic'),col=colours,pch=c(16,17),lty=1,bty='n',cex=.85)
}
x<-subset(p,task=='score'&n==10000&d==8&k==16&feedback=='diagonal'&observed&!direct)
draw(x$parameters,x$seconds_numeric,x$seconds_analytic,'Fitted parameters','8 dimensions, 16 kernels, 10,000 events')
x<-subset(p,task=='score'&n==10000&d==k&g==1&feedback=='diagonal'&observed&!direct)
draw(x$d,x$seconds_numeric,x$seconds_analytic,'Dimensions (= kernel columns)','Shared parameters: P = 3')
x<-subset(p,task=='score'&n==10000&d==4&g==1&feedback=='diagonal'&observed&!direct)
draw(x$k,x$seconds_numeric,x$seconds_analytic,'Kernel columns','4 dimensions, P = 3')
x<-subset(m,task=='score'&n==10000&d==8&k==16&feedback=='diagonal'&observed&!direct)
a<-x[x$method=='analytic',];b<-x[x$method=='numeric',];a<-a[order(a$g),];b<-b[order(b$g),]
stopifnot(identical(a$g,b$g))
plot(3*b$g,b$allocated_bytes/2^20,type='b',pch=16,col=colours[1],xaxt='n',xlab='Fitted parameters',
  ylab='Recorded R allocation (MiB)',ylim=range(x$allocated_bytes/2^20),main='8 dimensions, 16 kernels, 10,000 events')
axis(1,at=3*b$g,labels=3*b$g);lines(3*a$g,a$allocated_bytes/2^20,type='b',pch=17,col=colours[2])
legend('topleft',c('Numeric','Analytic'),col=colours,pch=c(16,17),lty=1,bty='n',cex=.85)
mtext('Times: medians of 3 isolated runs. R allocations exclude native standard-library buffers.',side=1,outer=TRUE,cex=.8)
invisible(dev.off())
