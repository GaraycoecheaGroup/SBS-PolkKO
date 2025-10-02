#!/usr/bin/env Rscript
#Yang Jiang 30.07.2025
#for plot the mpileup output from MutSeq analysis pipeline
#need readr, ggplot, reshape2 installed
#could be used in command line

library(readr)
library(tibble)
library(reshape2)
library(ggplot2)
library(stringr)
library(dplyr)
library(cowplot)
library(grid)
library(gridExtra)

Sys.setenv(VROOM_CONNECTION_SIZE=100000000)

#summerise all mpileup files ####
count_chars <- function(line) {
  cleaned_line <- gsub("[+-]?\\d+", "", line)
  cleaned_line <- gsub("\\*", "X", cleaned_line)
  cleaned_line <- gsub(",",".", cleaned_line)
  chars <- strsplit(cleaned_line, NULL)[[1]]
  
  char_count <- table(chars)
  df <- as.data.frame(as.list(char_count))
  
  return(df)
}

fill_empty <- function(df) {
  missing_cols <- setdiff(all_chars, names(df))
  df[missing_cols] <- 0
  df <- df[all_chars] 
  return(df)
}

file.path <- "/path/to/fig5f_mpileup"
position.modification <- 3287

file.list <- list.files(path = file.path,pattern = ".txt.gz",full.names = TRUE)
sample_names <- list.files(path = file.path,pattern = ".txt.gz",full.names = FALSE)
#sample_names <- str_sub(sample_names,1,-19)

sample_names <- gsub(".mpileup.txt.gz", "", sample_names)
df.correction <- read_delim("/path/to/Fig5f_correction_factors.tsv",delim = "\t",col_names = TRUE)


plot.order <- c("pM_AA_Mock_180_exp2","pM_AA_DeltaPolK_180_exp2",
                "pM_ACR_Mock_180_exp2","pM_ACR_DeltaPolK_180_exp2",
                "pCtrl_Mock","pCtrl_Polk")


df.mut.all <- data.frame()
df.plot.all <- data.frame()
df.depth <- data_frame()

for (i in 1:length(file.list)){
  
#for per sample analysis:
df <- read_delim(file.list[i],delim = "\t",col_names = FALSE,show_col_types = FALSE)
colnames(df) <- c("chr","position","REF","depth","pile","score")

list_of_dfs <- lapply(df$pile, count_chars)
all_chars <- unique(unlist(lapply(list_of_dfs, names)))
list_of_dfs <- lapply(list_of_dfs, fill_empty)
char_counts_df <- do.call(rbind, list_of_dfs)
full_count <- cbind(df[,1:4],char_counts_df)
colnames(full_count)[c(5,6)] <- c("same_as_ref","indel")

char_percent_df <- as_tibble(100*(char_counts_df/rowSums(char_counts_df)))
char_percent_df_no_ref <- char_percent_df[,colnames(char_counts_df)!="."]

df.plot <- as.data.frame(1:nrow(char_percent_df))
df.plot <- cbind(df.plot,df[,2:3])
df.plot <- cbind(df.plot,char_percent_df_no_ref)
colnames(df.plot)[c(1,2)] <- c("position","chr_position")
colnames(df.plot)[colnames(df.plot)=="X"] <- "indel"
 
df.plot <- melt(df.plot,id.vars = 1:3,variable.name = "type",value.name = "percent")
df.plot$percent_corrected <- df.plot$percent/df.correction$Correction_factor[i]
df.plot$sample <- sample_names[i]

df.plot$type <- factor(df.plot$type,levels = c("A","C","G","T","indel"))
df.plot <- df.plot[order(df.plot$type),]

df.depth.temp <- data.frame(round(sum(df$depth)/nrow(df)))
colnames(df.depth.temp) <- "mean_depth"
df.depth.temp$sample <- sample_names[i]
  
if (i==1) {
  df.plot.all <- df.plot
  df.depth <- df.depth.temp
} else { 
  df.plot.all <- rbind(df.plot.all,df.plot)
  df.depth <- rbind(df.depth,df.depth.temp)
}

#for grouped analysis:
df.plot.distal <- df.plot[df.plot$chr_position<(position.modification-20) | df.plot$chr_position>(position.modification+20),]
df.plot.proximal <- df.plot[df.plot$chr_position>=(position.modification-20) & df.plot$chr_position<=(position.modification+20),]

df.mut.pre.5 <- df.plot.proximal[df.plot.proximal$type!="indel" & df.plot.proximal$chr_position<position.modification & df.plot.proximal$chr_position>=position.modification-5,]
df.mut.pre.5$position <- "proximal.pre.1-5bp"
df.mut.pre.15 <- df.plot.proximal[df.plot.proximal$type!="indel" & df.plot.proximal$chr_position<position.modification-5,]
df.mut.pre.15$position <- "proximal.pre.6-20bp"
df.mut.mod <- df.plot.proximal[df.plot.proximal$type!="indel" & df.plot.proximal$chr_position==position.modification,]
df.mut.mod$position <- "modified"
df.mut.post <- df.plot.proximal[df.plot.proximal$type!="indel" & df.plot.proximal$chr_position>position.modification,]
df.mut.post$position <- "proximal.post"
df.mut.distal <- df.plot.distal[df.plot.proximal$type!="indel",]
df.mut.distal$position <- "distal"

df.mut <- rbind(df.mut.pre.5,df.mut.pre.15,df.mut.mod,df.mut.post,df.mut.distal)

#write_delim(df.mut,paste(file.path,sample_names[i],"_perbase_mut_percentage.tsv",sep = ""),delim = "\t",col_names = TRUE,)

df.mut$mut <- paste(df.mut$REF,">",df.mut$type,sep = "")

df.mut$mut2 <- df.mut$mut
df.mut$mut2.count <- 1

df.mut$mut2[df.mut$mut2=="G>A"] <- "C>T"
df.mut$mut2[df.mut$mut2=="G>C"] <- "C>G"
df.mut$mut2[df.mut$mut2=="G>T"] <- "C>A"
df.mut$mut2[df.mut$mut2=="A>C"] <- "T>G"
df.mut$mut2[df.mut$mut2=="A>G"] <- "T>C"
df.mut$mut2[df.mut$mut2=="A>T"] <- "T>A"

df.mut.simple.a <- aggregate(percent ~ mut2 + position, data = df.mut, sum)
df.mut.simple.b <- aggregate(mut2.count  ~ mut2 + position, data = df.mut, sum)
df.mut.simple.c <- aggregate(percent_corrected ~ mut2 + position, data = df.mut, sum)

#df.mut.simple <- left_join(df.mut.simple.a,df.mut.simple.b)
df.mut.simple <- left_join(df.mut.simple.a,df.mut.simple.c)
df.mut.simple <- left_join(df.mut.simple,df.mut.simple.b)

df.mut.simple$sample <- sample_names[i]
df.mut.simple$position <- factor(df.mut.simple$position, levels = c("proximal.pre.1-5bp","proximal.pre.6-20bp","modified","proximal.post","distal"))
df.mut.simple$mut2 <- factor(df.mut.simple$mut2, levels = c("C>A","C>G","C>T","T>A","T>C","T>G","A>indel","C>indel","G>indel","T>indel"))
df.mut.simple <- df.mut.simple[!is.na(df.mut.simple$mut2),]
df.mut.simple <- df.mut.simple[order(df.mut.simple$position),]

df.mut.simple <- df.mut.simple %>%
  group_by(position) %>%
  mutate(SumPercent = sum(percent)) 

#percentage of the mutation type in each group
df.mut.simple$compostion <- 100*df.mut.simple$percent/df.mut.simple$SumPercent

if (i==1) {df.mut.all <- df.mut.simple} 
else { 
  df.mut.all <- rbind(df.mut.all,df.mut.simple)
}
}

#write_delim(df.plot.all,paste(file.path,"PerBase_mut_percentage_allSamples.tsv",sep = ""),delim = "\t",col_names = TRUE,)
#write_delim(df.mut.all,paste(file.path,"PerGroup_mut_percentage_allSamples.tsv",sep = ""),delim = "\t",col_names = TRUE,)


df.plot.proximal.10 <- df.plot.all[df.plot.all$chr_position>=(position.modification-10) & df.plot.all$chr_position<=(position.modification+10),]
df.max <- aggregate(percent ~ sample+position, data = df.plot.proximal.10, sum)
yMax <- ceiling(max(df.max$percent)*100)/100
df.max_cor <- aggregate(percent_corrected ~ sample+position, data = df.plot.proximal.10, sum)
yMax_cor <- ceiling(max(df.max_cor$percent_corrected)*100)/100


p.l.10_withID <- list()
p.l.10_withoutID <- list()
p.l.10_withID_cor <- list()
p.l.10_withoutID_cor <- list()

#Plot
for (i in 1:length(plot.order)){

p.l.10_withID[[i]] <- ggplot(df.plot.proximal.10[df.plot.proximal.10$sample==plot.order[i],],aes(x=position,y=percent,fill=type)) +
  geom_col(position = "stack",) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=15),axis.title.x = element_text(size = 15))+
  scale_x_continuous(n.breaks = 5)+
  labs(title = paste(plot.order[i], "with InDel"),x=paste("average depth per base:",df.depth$mean_depth[df.depth$sample==plot.order[i]],sep = " "),
       fill=NULL,y="mutation rate (%)")+
  scale_y_continuous(limits = c(-0.1 ,yMax))+
  annotate("text", x = df.plot.proximal.10$position, y =-0.1, label = df.plot.proximal.10$REF, vjust = 1) +
  scale_fill_manual(values = c("A" = "#25A0DA", "C" = "#AFAFAF", "G" = "#040404","T"="#D1352E", "indel"="#FFD14A")) +
  geom_segment(aes(x = 113.5, xend = 115.5, y = -0.05, yend = -0.05), color = "purple", linewidth = 1.5)

ggsave(filename = paste(file.path,"/withID/",plot.order[i],"_withID_10bp_flanking.pdf",sep = ""),plot = p.l.10_withID[[i]],device = "pdf",width = 5,height = 5,units = "in",dpi = 600,create.dir = TRUE,bg = "white")

p.l.10_withoutID[[i]] <- ggplot(df.plot.proximal.10[df.plot.proximal.10$sample==plot.order[i] & df.plot.proximal.10$type!="indel",],aes(x=position,y=percent,fill=type)) +
  geom_col(position = "stack",) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=15),axis.title.x = element_text(size = 15))+
  scale_x_continuous(n.breaks = 5)+
  labs(title = paste(plot.order[i], "w/o InDel"),x=paste("average depth per base:",df.depth$mean_depth[df.depth$sample==plot.order[i]],sep = " "),
       fill=NULL,y="mutation rate (%)")+
  scale_y_continuous(limits = c(-0.1 ,yMax))+
  annotate("text", x = df.plot.proximal.10$position, y =-0.1, label = df.plot.proximal.10$REF, vjust = 1) +
  scale_fill_manual(values = c("A" = "#25A0DA", "C" = "#AFAFAF", "G" = "#040404","T"="#D1352E", "indel"="#FFD14A")) +
  geom_segment(aes(x = 113.5, xend = 115.5, y = -0.05, yend = -0.05), color = "purple", linewidth = 1.5)

ggsave(filename = paste(file.path,"/withoutID/",plot.order[i],"_withoutID_10bp_flanking.pdf",sep = ""),plot = p.l.10_withoutID[[i]],device = "pdf",width = 5,height = 5,units = "in",dpi = 600,create.dir = TRUE,bg = "white")

p.l.10_withID_cor[[i]] <- ggplot(df.plot.proximal.10[df.plot.proximal.10$sample==plot.order[i],],aes(x=position,y=percent_corrected,fill=type)) +
  geom_col(position = "stack",) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=15),axis.title.x = element_text(size = 15))+
  scale_x_continuous(n.breaks = 5)+
  labs(title = paste("Corrected",plot.order[i], "with InDel"),x=paste("average depth per base:",df.depth$mean_depth[df.depth$sample==plot.order[i]],sep = " "),
       fill=NULL,y="corrected mutation rate (%)")+
  scale_y_continuous(limits = c(-0.1 ,yMax_cor))+
  annotate("text", x = df.plot.proximal.10$position, y =-0.1, label = df.plot.proximal.10$REF, vjust = 1) +
  scale_fill_manual(values = c("A" = "#25A0DA", "C" = "#AFAFAF", "G" = "#040404","T"="#D1352E", "indel"="#FFD14A")) +
  geom_segment(aes(x = 113.5, xend = 115.5, y = -0.05, yend = -0.05), color = "purple", linewidth = 1.5)

ggsave(filename = paste(file.path,"/withID_corrected/",plot.order[i],"_withID_10bp_flanking_corrected.pdf",sep = ""),plot = p.l.10_withID_cor[[i]],device = "pdf",width = 5,height = 5,units = "in",dpi = 600,create.dir = TRUE,bg = "white")

p.l.10_withoutID_cor[[i]] <- ggplot(df.plot.proximal.10[df.plot.proximal.10$sample==plot.order[i] & df.plot.proximal.10$type!="indel",],aes(x=position,y=percent_corrected,fill=type)) +
  geom_col(position = "stack",) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=15),axis.title.x = element_text(size = 15))+
  scale_x_continuous(n.breaks = 5)+
  labs(title = paste("Corrected",plot.order[i], "w/o InDel"),x=paste("average depth per base:",df.depth$mean_depth[df.depth$sample==plot.order[i]],sep = " "),
       fill=NULL,y="corrected mutation rate (%)")+
  scale_y_continuous(limits = c(-0.1 ,yMax_cor))+
  annotate("text", x = df.plot.proximal.10$position, y =-0.1, label = df.plot.proximal.10$REF, vjust = 1) +
  scale_fill_manual(values = c("A" = "#25A0DA", "C" = "#AFAFAF", "G" = "#040404","T"="#D1352E", "indel"="#FFD14A")) +
  geom_segment(aes(x = 113.5, xend = 115.5, y = -0.05, yend = -0.05), color = "purple", linewidth = 1.5)

ggsave(filename = paste(file.path,"/withoutID_corrected/",plot.order[i],"_withoutID_10bp_flanking_corrected.pdf",sep = ""),plot = p.l.10_withoutID_cor[[i]],device = "pdf",width = 5,height = 5,units = "in",dpi = 600,create.dir = TRUE,bg = "white")


}
