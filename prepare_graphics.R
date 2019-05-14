err.bar <- round(3*sqrt(.05*.95/(n.replications)),3) # confidence limits of p-value
.limits <- 0.05 + c(-err.bar,0,err.bar)

pvals.1.9 %<>%  as.data.table()
pvals.melt <- melt(pvals.1.9[,!'replication'], 
                   id.vars=c("effect"), 
                   variable.name='statistic') 
pvals.melt[,c("reject","effect.factor"):=list(as.numeric(value <= 0.05), as.factor(effect)),] 

# Reorder levels
pvals.melt$statistic <- factor(pvals.melt$statistic)

# Filter statistics
pvals.melt <- pvals.melt[!is.na(statistic),,]

plot.3 <- pvals.melt %>% 
  ggplot(aes(y=reject, x=statistic, group=effect.factor, shape=effect.factor, color=effect.factor)) +
  theme_bw(base_size = 20)+
  theme(legend.position="none")+
  # ggtitle("Fixed signal, Gaussian Noise")+
  ylab('Power')+ 
  xlab('')+
  ylim(0,1)+
  stat_summary(fun.y='mean', geom="point", cex=4) +
  geom_hline(yintercept=0.05, lty=2)+
  geom_vline(xintercept=.limits, lty=c(3,2,3))+
  coord_flip()
plot(plot.3)

pvals.melt[,.(Power=mean(value<0.05)),.(effect,statistic)]
