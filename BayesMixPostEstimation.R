rm(list=ls(all=T))
setwd('~/Documents/GitHub/capture-het/Code/')
if (!require('pacman')){install.packages('pacman')}
pacman::p_load(R.matlab,coda,reshape2,ggplot2,scales,gridExtra,tidyr,dplyr)

dat <- readMat('Output/betamix.mat')
nmc <- length(dat$Nhat)

df <- data.frame(iter=seq(nmc),A=dat$A,B=dat$B,Nhat=dat$Nhat,eta=dat$eta)
p <- ncol(df)-1

ES <- matrix(0,p,1)
for (j in 1:p) {
  ES[j] <- effectiveSize(as.mcmc(df[,j+1]))
}

df <- melt(df,id='iter')
ggplot(df,aes(x=value)) + geom_histogram() + facet_grid(~variable,scales='free')

dat <- readMat('Output/betamix2.mat')
nmc <- length(dat$Nhat)

df <- data.frame(iter=seq(nmc),A=dat$A,B=dat$B,Nhat=dat$Nhat,eta=dat$eta)
p <- ncol(df)-1

ES <- matrix(0,p,1)
for (j in 1:p) {
  ES[j] <- effectiveSize(as.mcmc(df[,j+1]))
}

df <- melt(df,id='iter')
ggplot(df,aes(x=value)) + geom_histogram() + facet_grid(~variable,scales='free')

dat <- readMat('Output/betamix_full.mat')
nmc <- length(dat$Nhat)

df <- data.frame(iter=seq(nmc),A=dat$A,B=dat$B,Nhat=dat$Nhat,eta=dat$eta)
p <- ncol(df)-1

ES <- matrix(0,p,1)
for (j in 1:p) {
  ES[j] <- effectiveSize(as.mcmc(df[,j+1]))
}

df <- melt(df,id='iter')
ggplot(df,aes(x=iter,y=value)) + geom_point() + facet_wrap(~variable,scales='free')


# simulations for paper
nmc <- 100000
nsim <- 5
N <- 500
pob <- c(1,1,1,1)
pob <- cbind(pob,c(1,1,.7,.7))
#pob <- cbind(pob,c(0.999582, 0.999163, 0.995694, 0.990052))
pob <- cbind(pob,c(0.981601, 0.973948, 0.94118, 0.915762))
#pob <- cbind(pob,c(0.954623, 0.932935, 0.868657, 0.831888))
pob <- cbind(pob,c(0.991681, 0.983389, 0.918063, 0.838953))
#pob <- cbind(pob,c(0.999165, 0.998326, 0.991488, 0.982593))
pob <- cbind(pob,c(0.999165, 0.998326, 0.991488, 0.982593))


for (j in 1:5) {
  dat <- readMat(paste('Output/bayes_results_case_',j,'.mat',sep=''))
  Nhat <- cbind(dat$Nhat,dat$Nhata)
  Nhat.atom <- cbind(dat$Nhat.atom,dat$Nhata.atom)
  if (j==1) {
    Nhat.all <- Nhat
    Nhat.atom.all <- Nhat.atom
    nms <- rep(j,nmc)
  } else {
    Nhat.all <- rbind(Nhat.all,Nhat)
    Nhat.atom.all <- rbind(Nhat.atom.all,Nhat.atom)
    nms <- c(nms,rep(j,nmc))
  }
}

df.N <- data.frame(Nhat.all)
#names(df.N) <- c('N',paste('Na-',c(0.005,0.01,0.05,0.1),sep=''))
#names(df.N) <- c('N',as.character(c(0.005,0.01,0.05,0.1)))
names(df.N) <- as.character(c(0,0.005,0.01,0.05,0.1))
df.N$sim <- nms
df.N <- melt(df.N,id='sim')

df.N.atom <- data.frame(Nhat.atom.all)
#names(df.N.atom) <- c('N',paste('Na-',c(0.005,0.01,0.05,0.1),sep=''))
#names(df.N.atom) <- c('N',as.character(c(0.005,0.01,0.05,0.1)))
names(df.N.atom) <- as.character(c(0,0.005,0.01,0.05,0.1))
df.N.atom$sim <- nms
df.N.atom <- melt(df.N.atom,id='sim')

df.Ntr <- t(pob*N)
df.Ntr <- cbind(N,df.Ntr)
df.Ntr <- data.frame(df.Ntr)
#names(df.Ntr) <- c('N',paste('Na-',c(0.005,0.01,0.05,0.1),sep=''))
#names(df.Ntr) <- c('N',as.character(c(0.005,0.01,0.05,0.1)))
names(df.Ntr) <- as.character(c(0,0.005,0.01,0.05,0.1))
df.Ntr$sim <- seq(nsim)
df.Ntr <- melt(df.Ntr,id='sim')

df.N <- df.N %>% group_by(variable,sim) %>% mutate(q975=quantile(value,0.975),q025=quantile(value,0.025))
df.qs <- df.N %>% group_by(variable,sim) %>% summarize(q975=quantile(value,0.975),q025=quantile(value,0.025))


png('Figures/bayes-beta-n-4.png',width=500,height=300)
ggplot(df.N[df.N$variable=='0' & df.N$sim==4,],aes(x=value)) + geom_histogram() + #xlim(c(0,5000)) + 
  geom_vline(xintercept=rep(N,5),col=I(2)) + scale_x_log10() +
  geom_vline(data=df.qs[df.qs$variable=='N' & df.qs$sim==4,],aes(xintercept=q975)) +
  geom_vline(data=df.qs[df.qs$variable=='N' & df.qs$sim==4,],aes(xintercept=q025)) + 
  theme(text=element_text(size=30)) + xlab('N')
#xlim(0,2000)
dev.off()

png('Figures/bayes-beta-nalpha-4.png',width=500,height=300)
ggplot(df.N[df.N$variable=='0.01' & df.N$sim==4,],aes(x=value)) + geom_histogram() + #xlim(c(0,5000)) + 
  geom_vline(xintercept=rep(N,5),col=I(2)) + scale_x_log10() +
  geom_vline(data=df.qs[df.qs$variable=='0.01' & df.qs$sim==4,],aes(xintercept=q975)) +
  geom_vline(data=df.qs[df.qs$variable=='0.01' & df.qs$sim==4,],aes(xintercept=q025)) + 
  theme(text=element_text(size=30)) + xlim(0,1500) + xlab(expression(N[alpha]))
#xlim(0,2000)
dev.off()


df.N$value <- log10(df.N$value)
df.qs$q975 <- log10(df.qs$q975)
df.qs$q025 <- log10(df.qs$q025)
df.Ntr$value <- log10(df.Ntr$value)
png('Figures/bayes-beta-nalpha.png',width=900,height=400)
ggplot(df.N,aes(x=variable,y=value)) + 
  geom_violin(draw_quantiles = c(0.025,0.975)) + facet_wrap(~sim,scales='free',nrow=1) +
  geom_point(data=df.qs,aes(x=variable,y=q975),shape=2,size=3) +
  geom_point(data=df.qs,aes(x=variable,y=q025),shape=2,size=3) +
  geom_point(data=df.Ntr,aes(x=variable,y=value),size=3) + 
  theme(text=element_text(size=24),axis.text.x=element_text(angle=90)) +
  ylab(expression(log[10](N[alpha]))) + xlab(expression(alpha)) #+
  #ylim(c(2,4))
  #scale_y_log10(labels=trans_format('log10',math_format(10^.x)),limits=c(100,10000))
dev.off()


df.N.atom <- df.N.atom %>% group_by(variable,sim) %>% mutate(q975=quantile(value,0.975),q025=quantile(value,0.025))
df.qs.atom <- df.N.atom %>% group_by(variable,sim) %>% summarize(q975=quantile(value,0.975),q025=quantile(value,0.025))

png('Figures/bayes-atom-n-1.png',width=900,height=400)
ggplot(df.N.atom[df.N.atom$variable=='N' & df.N.atom$sim==1,],aes(x=value)) + geom_histogram() + #facet_grid(~sim,scales='free') +
  theme(text=element_text(size=30),axis.text.x=element_text(angle=90)) + #xlim(c(0,5000)) + 
  geom_vline(xintercept=rep(N,5),col=I(2)) + #scale_x_log10()
  geom_vline(data=df.qs.atom[df.qs.atom$variable=='N' & df.qs$sim==1,],aes(xintercept=q975)) +
  geom_vline(data=df.qs.atom[df.qs.atom$variable=='N' & df.qs$sim==1,],aes(xintercept=q025))
#xlim(0,2000)
dev.off()

df.N.atom$value <- log10(df.N.atom$value)
df.qs.atom$q975 <- log10(df.qs.atom$q975)
df.qs.atom$q025 <- log10(df.qs.atom$q025)

png('Figures/bayes-atom-nalpha.png',width=900,height=400)
#ggplot(df.N.atom[df.N.atom$variable!='N',],aes(x=variable,y=value)) +
ggplot(df.N.atom,aes(x=variable,y=value)) + 
  geom_violin(draw_quantiles = c(0.025,0.975)) + facet_wrap(~sim,scales='free',nrow=1) +
  geom_point(data=df.qs.atom,aes(x=variable,y=q975),shape=2,size=3) +
  geom_point(data=df.qs.atom,aes(x=variable,y=q025),shape=2,size=3) +
  geom_point(data=df.Ntr,aes(x=variable,y=value),size=3) + 
  theme(text=element_text(size=24),axis.text.x=element_text(angle=90)) +
  ylab(expression(log[10](N[alpha]))) + xlab(expression(alpha)) #+
  #scale_y_log10(labels=trans_format('log10',math_format(10^.x)))
  #ylim(0,1250)
dev.off()

# png('Figures/bayes-beta-n.png',width=600,height=400)
# ggplot(df.N[df.N$variable=='N',],aes(x=log(value,10))) + geom_histogram() + facet_grid(~sim,scales='free') +
#   theme(text=element_text(size=24),axis.text.x=element_text(angle=90)) + #xlim(c(0,5000)) + 
#   geom_vline(xintercept=rep(log(N,10),5),col=I(2)) +#+ scale_x_log10(breaks=c(100,1000,10000,100000),labels=trans_format('log10',math_format(10^.x))) +
#   xlab(expression(log[10](N))) 
#   #xlim(0,2000)
# dev.off()
# 
# png('Figures/bayes-atom-n.png',width=600,height=400)
# ggplot(df.N.atom[df.N.atom$variable=='N',],aes(x=value)) + geom_histogram() + facet_grid(~sim,scales='free') +
#   theme(text=element_text(size=24),axis.text.x=element_text(angle=90)) + #xlim(c(0,5000)) + 
#   geom_vline(xintercept=rep(N,5),col=I(2)) + xlab("N")#+ scale_x_log10()
#   #xlim(0,2000)
# dev.off()



