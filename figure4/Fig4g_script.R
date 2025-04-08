# load packages 
library(readr)
library(ggplot2)
library(stringr)
library(dplyr)
library(gridExtra)

#load data

TSB.df <- read_delim("TSB_significance.txt",delim = "\t")

genotype_level <- c('WT','Xpa','Xpc','Polk','Xpa_Polk',"Xpc_Polk")

#for setting significance groups
limits2 <- function(pval) {
  sig <- vector(mode='numeric',length=length(pval))
  for (i in seq(1,length(pval))){
    if (pval[i]>0.05 ){sig[i] <- 0} 
    else if (pval[i]>0.01 & pval[i]<=0.05 ){sig[i] <- 0.25}		# P < 0.05
    else if (pval[i]>0.001 & pval[i]<=0.01 ){sig[i] <- 0.5}	  # P < 0.01
    else if (pval[i]>0.0001 & pval[i]<=0.001 ){sig[i] <- 0.5} # P < 0.001
    else if (pval[i]<=0.0001 ){sig[i]=1}	                    # P < 0.0001
  }
  return(sig)
}


tra <- TSB.df %>% mutate(significance_p=limits2(p_value))
tra$log2_ratio <- log2(tra$real_ratio)
tra$log2_ratio[is.na(tra$log2_ratio)] <- 0

tra$genotype <- factor(tra$genotype,levels = rev(genotype_level))

t1.rect1 <- data.frame (xmin=1.5, xmax=2.5, ymin=0.5, ymax=6.5)
t2.rect1 <- data.frame (xmin=3.5, xmax=4.5, ymin=0.5, ymax=6.5)
t3.rect1 <- data.frame (xmin=5.5, xmax=6.5, ymin=0.5, ymax=6.5)

spot.theme.1 <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 13)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 12)),
  theme(axis.line=element_blank()),
  theme(title = element_text(size = 10)),
  theme(text = element_text(size = 10)),
  theme(panel.background = element_rect(fill = 'white')),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  theme(legend.box.background = element_rect(color='white')),
  scale_size_continuous(range = c(-1, 8)),
  scale_colour_gradient(low = "white",high = "red",guide = "colourbar",aesthetics = "colour"),
  scale_x_discrete(position = "top"))

spot.theme.2 <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 13)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 12)),
  theme(axis.line=element_blank()),
  theme(title = element_text(size = 10)),
  theme(text = element_text(size = 10)),
  theme(panel.background = element_rect(fill = 'white')),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  theme(legend.box.background = element_rect(color='white')),
  scale_size_continuous(range = c(-1, 8)),
  scale_colour_gradient(low = "white",high = "darkred",guide = "colourbar",aesthetics = "colour"),
  scale_fill_discrete(type = c("blue","red")),
  scale_x_discrete(position = "top"))

p1 <- ggplot(tra[tra$signature=="SBS-Polk",],aes(x=mutation_type, y=genotype)) + 
  geom_point(aes(colour = log2_ratio, size = significance_p))+
  labs(x="",y="",title="SBS POLK")+
  spot.theme.1+ 
  geom_rect(data=t1.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t2.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t3.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_hline(yintercept=1.5,color = "white",linewidth=1)+  geom_hline(yintercept=2.5,color = "white",linewidth=1)+geom_hline(yintercept=3.5,color = "white",linewidth=1)+
  geom_hline(yintercept=4.5,color = "white",linewidth=1)+  geom_hline(yintercept=4.5,color = "white",linewidth=1)+geom_hline(yintercept=5.5,color = "white",linewidth=1)+
  geom_hline(yintercept=6.5,color = "white",linewidth=1)+  geom_hline(yintercept=7.5,color = "white",linewidth=1)+geom_hline(yintercept=8.5,color = "white",linewidth=1)+
  geom_hline(yintercept=9.5,color = "white",linewidth=1)+  geom_hline(yintercept=10.5,color = "white",linewidth=1)+
  geom_hline(yintercept=11.5,color = "white",linewidth=1)+  geom_hline(yintercept=12.5,color = "white",linewidth=1)+geom_hline(yintercept=13.5,color = "white",linewidth=1)+
  geom_hline(yintercept=14.5,color = "white",linewidth=1)+  geom_hline(yintercept=15.5,color = "white",linewidth=1)+geom_hline(yintercept=16.5,color = "white",linewidth=1)+
  geom_hline(yintercept=17.5,color = "white",linewidth=1)+  geom_hline(yintercept=18.5,color = "white",linewidth=1)+geom_hline(yintercept=19.5,color = "white",linewidth=1)+
  geom_hline(yintercept=20.5,color = "white",linewidth=1)+  geom_hline(yintercept=21.5,color = "white",linewidth=1)+
  geom_vline(xintercept=0.5,color ="gray",linewidth=1)+geom_vline(xintercept=6.5,color = "gray",linewidth=1)


p2 <- ggplot(tra[tra$signature=="SBS-NER",],aes(x=mutation_type, y=genotype)) + 
  geom_point(aes(colour = log2_ratio, size = significance_p))+
  spot.theme.1+ 
  labs(x="",y="",title="SBS NER")+
  geom_rect(data=t1.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t2.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_rect(data=t3.rect1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='black', alpha=0.1, inherit.aes = FALSE) +
  geom_hline(yintercept=1.5,color = "white",linewidth=1)+  geom_hline(yintercept=2.5,color = "white",linewidth=1)+geom_hline(yintercept=3.5,color = "white",linewidth=1)+
  geom_hline(yintercept=4.5,color = "white",linewidth=1)+  geom_hline(yintercept=4.5,color = "white",linewidth=1)+geom_hline(yintercept=5.5,color = "white",linewidth=1)+
  geom_hline(yintercept=6.5,color = "white",linewidth=1)+  geom_hline(yintercept=7.5,color = "white",linewidth=1)+geom_hline(yintercept=8.5,color = "white",linewidth=1)+
  geom_hline(yintercept=9.5,color = "white",linewidth=1)+  geom_hline(yintercept=10.5,color = "white",linewidth=1)+
  geom_hline(yintercept=11.5,color = "white",linewidth=1)+  geom_hline(yintercept=12.5,color = "white",linewidth=1)+geom_hline(yintercept=13.5,color = "white",linewidth=1)+
  geom_hline(yintercept=14.5,color = "white",linewidth=1)+  geom_hline(yintercept=15.5,color = "white",linewidth=1)+geom_hline(yintercept=16.5,color = "white",linewidth=1)+
  geom_hline(yintercept=17.5,color = "white",linewidth=1)+  geom_hline(yintercept=18.5,color = "white",linewidth=1)+geom_hline(yintercept=19.5,color = "white",linewidth=1)+
  geom_hline(yintercept=20.5,color = "white",linewidth=1)+  geom_hline(yintercept=21.5,color = "white",linewidth=1)+
  geom_vline(xintercept=0.5,color ="gray",linewidth=1)+geom_vline(xintercept=6.5,color = "gray",linewidth=1)


pp <- grid.arrange(p1,p2,ncol=1,nrow=2)

ggsave("TSB_Significance.pdf",plot = pp,device = "pdf",width = 7,height = 8.5,units = "in")



