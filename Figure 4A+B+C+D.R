#load necessary packages
#install_div_mutants.packages("tibble")
library(tibble)
#install_div_mutants.packages("magrittr")
library(magrittr)
#install_div_mutants.packages("ggplot2")
library(ggplot2)
#install_div_mutants.packages("readxl")
library(readxl)
#install_div_mutants.packages("dplyr")
library(dplyr)

#set ggplot figure theme
migtheme <-   theme(axis.line = element_blank(),
                    axis.text = element_text(colour = "black", size = 12),
                    strip.text = element_text(colour = "black", size = 12),
                    legend.text = element_text(colour = "black", size = 12),
                    axis.title = element_text(colour = "black", size = 12),
                    axis.ticks = element_line(colour = "black", size = 0.3),
                    axis.ticks.length = unit(-0.1, "cm"),
                    legend.title = element_text(colour = "black", size = 12),
                    panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
                    panel.grid = element_blank(),
                    plot.background = element_blank(),
                    panel.background = element_blank(),
                    legend.background = element_blank(),
                    legend.key = element_rect(fill = "white", colour = "white"),
                    strip.background = element_blank(),
                    text = element_text(colour = "black", family = "sans"))


#load data
pig <- read_xlsx("pig-1_smFISH.xlsx")%>%
  add_column(label = "pig-1(gm344)")
pig$Image%<>%as.numeric

ced <- read_xlsx("ced-3_smFISH.xlsx")%>%
  add_column(label =  "ced-3(n717)")
ced$Image%<>%as.numeric

SV1009 <- read_xlsx("SV1009_smFISH.xlsx")%>%
  add_column(label =  "Control")
SV1009

nuclear_size <- read_xlsx("nuclear_size.xlsx")


all_div_mutants <- bind_rows(pig, ced, SV1009)

all_div_mutants$label <- factor(all_div_mutants$label, levels =  c("Control", "ced-3(n717)", "pig-1(gm344)"))

fig.4A <- ggplot(filter(all_div_mutants, Cell %in% c("QR.p", "QR.pa", "QR.pp")) , aes(x = Position, y = smFISH.Counts, colour = Cell))+
  geom_point(size = 2)+
  scale_x_continuous(labels = c("H2","V1", "V2", "V3", "V4", "V5"), limits = c(-1,4))+
  scale_y_continuous(limits = c(0,45))+
  geom_vline(xintercept = 1, colour = "black", lty = 2)+
  scale_color_manual(values = c("green3", "magenta", "darkorange1"))+
  labs(x = "position along AP axis", y = "# mRNA spots",colour = "Cell")+
  facet_wrap(~label, ncol = 3, scales = "free")+
  migtheme
fig.4A


fig.4B <- ggplot(filter(all_div_mutants, Cell =="QR.pa" & Ratio != 0) , aes(y = Ratio, x = label))+
  xlab(label = "")+
  ylab("QR.pa : QR.pp size ratio")+
  geom_boxplot(colour = "black", fill = "Gray80", size = 0.3)+
  stat_boxplot(geom = 'errorbar')+
  geom_jitter(colour = "black", width = 0.06, height = 0, size = 0.1)+
  migtheme
fig.4B

fig.4C <- ggplot(filter(pig, Cell =="QR.pa" & Ratio != 0) , aes(x = Ratio, y = smFISH.Counts))+
  geom_point(size = 2)+
  geom_smooth(method = "lm", formula = y~x, color = "black")+
  scale_color_manual(values = c("magenta"))+
  labs(x = "pig-1(gm344) QR.pa : QR.pp size ratio", y = "# mRNA spots",colour = "Cell")+
  theme(legend.position = "none")+
  migtheme
fig.4C

fig.4D <- ggplot(nuclear_size, aes(x = Cell, y = Area))+
  geom_boxplot(colour = "black", fill = "Gray80", size = 0.3)+
  stat_boxplot(geom = 'errorbar')+
  geom_jitter(colour = "black", width = 0.06, height = 0, size = 0.1)+
  scale_y_continuous(name = "Nuclear area (AU)",limits = c(0,0.9))+
  migtheme
fig.4D

#Statistical testing
# mig-1 expression level in QR.pa of:
#pig-1 mutant
t.test(filter(SV1009, Cell == "QR.pa")$smFISH.Counts, filter(pig, Cell == "QR.pa")$smFISH.Counts) #p = 7.034e-7
#ced-3 mutant
t.test(filter(SV1009, Cell == "QR.pa")$smFISH.Counts, filter(ced, Cell == "QR.pa")$smFISH.Counts) #p = 0.278
#Pearson's correlation of QR.pa:QR.pp size ratio and smFISH expression in pig-1 mutant
cor(filter(pig, Cell =="QR.pa" & Ratio != 0)$Ratio, filter(pig, Cell =="QR.pa" & Ratio != 0)$smFISH.Counts, method = "pearson") #R = 0.47
#nuclear sizes QR.p and QR.pa
t.test(filter(nuclear_size, Cell == "QR.p")$Area,filter(nuclear_size, Cell == "QR.pa")$Area) #p = 0.2225


