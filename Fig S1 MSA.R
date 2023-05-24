#load necessary packages
#install.packages("tibble")
library(tibble)
#install.packages("magrittr")
library(magrittr)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("readxl")
library(readxl)
#install.packages("dplyr")
library(dplyr)
meme_intron <- readRDS("MEME_intron.RDS")

#strsplit(meme_intron$Column2, split = "")%>%unlist%>%enframe%>%view
intron_msa <- list()
for(i in 1:4){
intron_msa[[i]] <- ggplot(meme_intron[[i]], aes(x = position, y = Column1, color = nucleotide, label = nucleotide))+
  geom_text(fontface = "bold" )+
  xlab("")+
  ylab("Species")+
  scale_color_discrete(type = c("red","blue","orange2", "green4"))+
  migtheme+
  theme(legend.position="none")+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks = element_blank())
}
#intron_msa[[1]]

intron_msa[[1]]


meme_promoter <- readRDS("MEME_promoter.RDS")

promoter_msa <- list()
for(i in 1:4){
  promoter_msa[[i]] <- ggplot(meme_promoter[[i]], aes(x = position, y = Column1, color = nucleotide, label = nucleotide))+
    geom_text(fontface = "bold" )+
    xlab("")+
    ylab("Species")+
    scale_color_discrete(type = c("red","blue","orange2", "green4"))+
    migtheme+
    theme(legend.position="none")+
    theme(axis.text.x = element_blank())+
    theme(axis.ticks = element_blank())
}
#promoter_msa[[1]]
