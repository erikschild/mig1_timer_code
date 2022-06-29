source("E:/OneDrive - Hubrecht Institute//Work/Functions/all_functional_functions.R")


meme_intron <- readRDS("E:/OneDrive - Hubrecht Institute/Work/mig-1_paper/figures/github version figure code/MEME_intron.RDS")

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


meme_promoter <- write_rds("E:/OneDrive - Hubrecht Institute/Work/mig-1_paper/figures/github version figure code/MEME_promoter.RDS")


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