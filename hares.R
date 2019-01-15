rm(list=ls(all=T))
setwd('~/Documents/GitHub/capture-het/Code')
if (!require('pacman')) {install.packages('pacman')} 
pacman::p_load(ggplot2,reshape2,Rcapture,R.matlab,dplyr,tidyr,xtable,scales)

data(hare)

n <- apply(hare,1,sum)
n <- table(n)
nmc <- 100000

# get results
dat <- readMat('Output/bayes_results_hares.mat')

df.beta <- data.frame(Nhat=dat$Nhat,Nhata=dat$Nhata)
names(df.beta) <- c(0,0.005,0.01,0.05,0.1)
df.beta$iter <- seq(nmc)
df.beta <- melt(df.beta,id='iter')
names(df.beta) <- c('iter','alpha','value')
beta.ptest <- df.beta %>% group_by(alpha) %>% summarise(Nhat=mean(value))
beta.ptest

df.atom <- data.frame(Nhat=dat$Nhat.atom,Nhata=dat$Nhata.atom)
names(df.atom) <- c(0,0.005,0.01,0.05,0.1)
df.atom$iter <- seq(nmc)
df.atom <- melt(df.atom,id='iter')
names(df.atom) <- c('iter','alpha','value')
atom.ptest <- df.atom %>% group_by(alpha) %>% summarise(Nhat=mean(value))
atom.ptest

ptest.tab <- cbind(beta.ptest,atom.ptest[,2])
colnames(ptest.tab) <- c('$\\alpha$','$\\wh N_{\\alpha}$ - Beta','$\\wh N_{\\alpha}$ - Discrete Mixture')
print(xtable(ptest.tab,caption = 'Point estimates of $N_{\\alpha}$ for snowshoe hare data',label = 'tab:hare-nhat',digits = 0),
      sanitize.text.function=function(x){x},include.rownames = FALSE,file='Figures/hare-nhat.tex',caption.placement='top')


df.beta.qs <- df.beta %>% group_by(alpha) %>% summarize(q975=quantile(value,0.975),q025=quantile(value,0.025))
df.atom.qs <- df.atom %>% group_by(alpha) %>% summarize(q975=quantile(value,0.975),q025=quantile(value,0.025))

df.beta$value <- log10(df.beta$value)
df.atom$value <- log10(df.atom$value)
df.beta.qs$q975 <- log10(df.beta.qs$q975)
df.beta.qs$q025 <- log10(df.beta.qs$q025)
df.atom.qs$q025 <- log10(df.atom.qs$q025)
df.atom.qs$q975 <- log10(df.atom.qs$q975)

ptest.tab <- data.frame(ptest.tab)
names(ptest.tab) <- c('alpha','Nhat_beta','Nhat_mix')

png('Figures/hares-beta-nalpha.png',width=500,height=300)
ggplot(df.beta,aes(x=alpha,y=value)) + geom_violin() +
  theme(text=element_text(size=32)) +
  xlab(expression(alpha)) + ylab(expression(log[10](N[alpha]))) +   #+ ggtitle('Beta')
  geom_point(data=ptest.tab,aes(x=alpha,y=log10(Nhat_beta)),size=3)
dev.off()

png('Figures/hares-atom-nalpha.png',width=500,height=300)
ggplot(df.atom,aes(x=alpha,y=value)) + geom_violin() +
  theme(text=element_text(size=32)) +
  xlab(expression(alpha)) + ylab(expression(log[10](N[alpha]))) + #+ ggtitle('Discrete mixture') + ylim(c(25,400))
  geom_point(data=ptest.tab,aes(x=alpha,y=log10(Nhat_mix)),size=3)
dev.off()

mnval.beta <- log(mean(df.beta[df.beta$alpha=='0',]$value))
# 
# png('Figures/hares-beta-n.png',width=500,height=300)
# ggplot(df.beta[df.beta$alpha=='0',],aes(x=value)) + geom_histogram() +
#   theme(text=element_text(size=32)) + scale_x_log10(labels=trans_format('log10',math_format(10^.x))) + 
#   geom_vline(data=df.beta.qs[df.beta.qs$alpha=='0',],aes(xintercept = q975)) +
#   geom_vline(data=df.beta.qs[df.beta.qs$alpha=='0',],aes(xintercept = q025)) + 
#   xlab('N') #+ xlim(0,10^5)
# dev.off()
# 
# png('Figures/hares-atom-n.png',width=500,height=300)
# ggplot(df.atom[df.atom$alpha=='0',],aes(x=value)) + geom_histogram() +
#   theme(text=element_text(size=32)) + scale_x_log10(labels=trans_format('log10',math_format(10^.x))) +
#   geom_vline(data=df.atom.qs[df.atom.qs$alpha=='0',],aes(xintercept = q975)) +
#   geom_vline(data=df.atom.qs[df.atom.qs$alpha=='0',],aes(xintercept = q025)) +
#   xlab('N') #+ xlim(0,300)
# dev.off()
