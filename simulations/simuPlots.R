library(ape)
library(mlesky)

# Bottleneck
rm(list=ls())
set.seed(2)
alphaFun=function(x){if (x<2005 || x>2010) 10 else 1}
sampleDates=seq(2000,2020,0.1)
tree=simCoal(sampleDates,alphaFun,alphaMin = 0.1)

res=optim_res_aic(tree,ncpu=6,res=1:50,model=2)
print(res)
fit=mlskygrid(tree,res=res,tau=NULL,tau_lower = 0.001,tau_upper = 10000,model = 2,ncpu=6)

pdf('simuBottle.pdf',7,10)
par(mfrow=c(3,1),mar=c(4,4,1,4))
if (!is.null(tree$root.time)) from=tree$root.time else from=-max(dist.nodes(tree)[Ntip(tree)+1,])
to=from+max(dist.nodes(tree)[Ntip(tree)+1,])
xs=seq(from=from,to=to,length.out=100)
ys=xs
for (i in 1:length(ys)) ys[i]=alphaFun(ys[i])
plot(xs,ys,type='l',xlab='',ylab='Population size', bty='l',ylim=c(0,1.05*max(ys)))
plot(tree,show.tip.label = F)
axisPhylo(1,backward = F)
plot(fit)
dev.off()


distance=function(alphaFun,fit) {
	mean((fit$ne-alphaFun(fit$time+max(sampleDates)))^2)
}

print(distance(alphaFun,fit))

# Const
rm(list=ls())
set.seed(1)
alphaFun=function(x){20}
sampleDates=seq(2000,2020,0.1)
t=simCoal(sampleDates,alphaFun,alphaMin = 0.1)
pdf('simuConst.pdf',7,7)
plotBoth(t,alphaFun)
dev.off()

res=optim_res_aic(t,ncpu=6,res=1:20,model=2)
print(res)

pdf('simuConstResult.pdf',7,10)
par(mfrow=c(3,3))
for (res in c(5,20,50)) for (tau in c(1,10,20)) {
	fit=mlskygrid(t,res=res,tau=tau,model = 2,ncpu=6)
	plot(fit,logy = F) #,ylim=c(0,100)
}
dev.off()

# Covar
# rm(list=ls())
# betas=sigmas=rep(NA,100)
# for (b in seq(1,length(betas),1)) {
# 	set.seed(b)
# 	print(b)
# 	x=seq(1990,2020,1/12)#monthly
# 	driver=-((x-2005))^2/200+0.5
# 	sigma=(floor((b-1)/10)+1)*0.2
# 	rho=driver+rnorm(length(driver),0,sigma)
# 	ne=rep(10,length(x))
# 	for (i in 2:length(x)) ne[i]=ne[i-1]*(1+rho[i-1]*(x[2]-x[1]))
# 	
# 	alphaFun=function(newx){approx(x,ne,newx,rule=2)$y}
# 	
# 	sampleDates=seq(2000,2020,0.1)
# 	tree=simCoal(sampleDates,alphaFun,alphaMin = max(1,min(ne)))
# 	
# 	if (b==11) {
# 		pdf('simuCovar.pdf',7,10)
# 		par(mfrow=c(4,1))
# 		plot(x,driver,type = 'l',xlab = '',ylab='Covariate')
# 		plot(x,rho,type = 'l',xlab = '',ylab='Growth rate')
# 		plot(x,ne,type='l',xlab = '',ylab='Effective population size')
# 		plot(tree,show.tip.label = F)
# 		axisPhylo(1,backward = F)
# 		dev.off()
# 	}
# 	
# 	covar=cbind(x,driver)
# 	colnames(covar)<-c('time','var')
# 	covar=as.data.frame(covar)
# 	sampleTimes=tree$root.time+dist.nodes(tree)[Ntip(tree)+1,]
# 	fit=mlskygrid(tree,res=30,tau=10,model = 2,ncpu=6,sampleTimes=sampleTimes,formula=~var,data=covar) #,formula_order=2
# 	betas[b]=fit$beta
# 	sigmas[b]=sigma
# }
# 
# pdf('simuCovar2.pdf',10,10)
# boxplot(betas~sigmas,ylab='Association coefficient (beta)',xlab='Noise in growth rate relative to the covariate data',outline=F)
# dev.off()

# Sinus
rm(list=ls())
set.seed(3)
alphaFun=function(x){sin(x)*10+12}
sampleDates=seq(2000,2020,0.1)
t=simCoal(sampleDates,alphaFun,alphaMin = 0.1)
pdf('simuSinus.pdf',7,7)
plotBoth(t,alphaFun)
dev.off()

fit=as.list(rep(NA,3))
for (m in c(1,2,3)) {
	fit[[m]]=mlskygrid(t,res=20,tau=NULL,tau_lower = 0.001,tau_upper = 10000,model = m,ncpu=6)
}

pdf('simuSinusResult.pdf',7,7)
par(mfrow=c(3,1))
for (m in c(2,3,1)) {
	plot(fit[[m]],logy = F) #,ylim=c(0,50)
}
dev.off()

distance=function(alphaFun,fit) {
	mean((fit$ne-alphaFun(fit$time+max(sampleDates)))^2)
}

distances=c()
for (m in c(2,3,1)) {
	distances=c(distances,distance(alphaFun,fit[[m]]))
}
print(distances)