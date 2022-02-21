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
KN3071 <- read_xlsx("KN3071_smFISH.xlsx")%>%
  add_column(label = "mig-1ΔIntron 1")

KN3078 <- read_xlsx("KN3078_smFISH.xlsx")%>%
  add_column(label = "mig-1ΔIntron 2")

SV1009 <- read_xlsx("SV1009_smFISH.xlsx")%>%
  add_column(label = "Control")

introndels <- bind_rows(KN3071, KN3078, SV1009)%>%
    filter(Cell %in% c("QR", "QR.p", "QR.pa"))
introndels$label <- factor(introndels$label, levels = c("Control", "mig-1ΔIntron 1", "mig-1ΔIntron 2"))

#Create figure panels
fig.2B <- ggplot(introndels , aes(x = Position, y = smFISH.Counts, colour = Cell))+
  geom_point(size = 2)+
  scale_x_continuous(labels = c("H2","V1", "V2", "V3", "V4", "V5"), limits = c(-1,4))+
  scale_y_continuous(limits = c(0,45))+
  geom_vline(xintercept = 1, colour = "black", lty = 2)+
  scale_color_manual(values = c("blue","green3", "magenta"))+
  labs(x = "position along AP axis", y = "# mRNA spots",colour = "Cell")+
  facet_wrap(~label, ncol = 3, scales = "fixed")+
  migtheme
fig.2B

fig.2C <- ggplot(introndels%>%filter(Cell == "QR"), aes(x = label, y = smFISH.Counts))+
  geom_boxplot(fill = "gray80", size = 0.3)+
  stat_boxplot(geom = 'errorbar')+
  geom_jitter(colour = "black", width = 0.06, height = 0, size = 0.1)+
  labs(x = "", y = "# mRNA spots QR")+
  migtheme
fig.2C

#Statistical testing
# mig-1 expression level in QR of:
#Δintron 1
t.test(filter(SV1009, Cell == "QR")$smFISH.Counts,filter(KN3071, Cell == "QR")$smFISH.Counts) #2.2e-16
#Δintron 2
t.test(filter(SV1009, Cell == "QR")$smFISH.Counts,filter(KN3078, Cell == "QR")$smFISH.Counts) #p = 0.3054

# mig-1 expression level in QR.pa of:
#Δintron 1
t.test(filter(SV1009, Cell == "QR.pa")$smFISH.Counts,filter(KN3071, Cell == "QR.pa")$smFISH.Counts) #p = 0.1411
#Δintron 2
t.test(filter(SV1009, Cell == "QR.pa")$smFISH.Counts,filter(KN3078, Cell == "QR.pa")$smFISH.Counts) #p = 0.6813