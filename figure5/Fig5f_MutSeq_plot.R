#!/usr/bin/env Rscript
#Yang Jiang 30.07.2025
#to plot the mpileup files


library(readr)
library(tibble)
library(reshape2)
library(ggplot2)
library(stringr)
library(dplyr)

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

file.path <- "/Users/yangjiang/Library/CloudStorage/OneDrive-HubrechtInstitute/Bioinformatic/MutSeq/202507_Koichi/fig5f_mpileup"
position.modification <- 3287

file.list <- list.files(path = file.path,pattern = ".txt.gz",full.names = TRUE)
sample_names <- list.files(path = file.path,pattern = ".txt.gz",full.names = FALSE)
#sample_names <- str_sub(sample_names,1,-19)

sample_names <- gsub(".mpileup.txt.gz", "", sample_names)
df.correction <- read_delim("/Users/yangjiang/Library/CloudStorage/OneDrive-HubrechtInstitute/Bioinformatic/MutSeq/202507_Koichi/Fig5f_correction_factors.tsv",delim = "\t",col_names = TRUE)

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
df.plot$percent_corrected <- df.plot$percent/as.numeric(df.correction[df.correction$sample_names==sample_names[i],2])
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

}


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

