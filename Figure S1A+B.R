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
SV1009 <- read_xlsx("SV1009_pax.xlsx")%>%
  add_column(Condition = "Control")

KN3133 <- read_xlsx("KN3133_pax.xlsx")%>%
  add_column(Condition = "mig-1ΔPromoter 1")

KN3134 <- read_xlsx("KN3134_pax.xlsx")%>%
  add_column(Condition = "mig-1ΔPromoter 2")

KN3071 <- read_xlsx("KN3071_pax.xlsx")%>%
  add_column(Condition = "mig-1ΔIntron 1")

KN3078 <- read_xlsx("KN3078_pax.xlsx")%>%
  add_column(Condition = "mig-1ΔIntron 2")

KN1500 <- read_xlsx("KN1500_pax.xlsx")%>%
  add_column(Condition = "mig-1 loss-of-function")

KN3222 <- read_xlsx("KN3222_pax.xlsx")%>%
  add_column(Condition = "mig-1ΔPromoter 1-2")

all_QR.pap <- bind_rows(SV1009, KN1500, KN3133, KN3134, KN3071, KN3078, KN3222)
all_QR.pap$Condition <- factor(all_QR.pap$Condition, levels = rev(c("Control", "mig-1ΔIntron 1", "mig-1ΔIntron 2",  "mig-1ΔPromoter 1", "mig-1ΔPromoter 2", "mig-1ΔPromoter 1-2", "mig-1 loss-of-function")))

#Create figure panels
fig.S1A <- ggplot(filter(all_QR.pap, Condition %in% c("Control", "mig-1ΔIntron 1", "mig-1ΔIntron 2", "mig-1 loss-of-function") & Cell == "QR.pap"), aes(x = Condition, y = Position))+
  geom_boxplot(colour = "gray10", fill = "gray50")+
  geom_jitter(colour = "magenta", width = 0.1, height = 0, size = 1, alpha = 0.6)+
  scale_y_continuous(limits = c(0.5,2),breaks = c(0,1,2), labels = c("H2", "V1.p", "V2.p"),name = "final position QR.pap on AP axis")+
  coord_flip()+
  migtheme
fig.S1A

fig.S1B <- ggplot(filter(all_QR.pap, Condition %in% c("Control", "mig-1ΔPromoter 1", "mig-1ΔPromoter 2",  "mig-1ΔPromoter 1-2", "mig-1 loss-of-function") & Cell == "QR.pap"), aes(x = Condition, y = Position))+
  geom_boxplot(colour = "gray10", fill = "gray50")+
  geom_jitter(colour = "magenta", width = 0.1, height = 0, size = 1, alpha = 0.6)+
  scale_y_continuous(limits = c(0.5,2),breaks = c(0,1,2), labels = c("H2", "V1.p", "V2.p"), name = "final position QR.pap on AP axis")+
  coord_flip()+
  migtheme
fig.S1B

#statistical testing
#final position QR.pap of:
#mig-1ΔIntron 1
t.test(x= filter(SV1009, Cell == 'QR.pap')$Position, y = filter(KN3071, Cell == 'QR.pap')$Position) #p = 0.2527
#mig-1ΔIntron 2
t.test(x= filter(SV1009, Cell == 'QR.pap')$Position, y = filter(KN3078, Cell == 'QR.pap')$Position) #p = 0.1723
#mig-1ΔPromoter 1
t.test(x= filter(SV1009, Cell == 'QR.pap')$Position, y = filter(KN3133, Cell == 'QR.pap')$Position) #p = 0.001743
#mig-1ΔPromoter 2
t.test(x= filter(SV1009, Cell == 'QR.pap')$Position, y = filter(KN3134, Cell == 'QR.pap')$Position) #p =  1.101E-6
#Pmig 1234
t.test(x= filter(SV1009, Cell == 'QR.pap')$Position, y = filter(KN3222, Cell == 'QR.pap')$Position) #p =  3.482E-13
#mig-1 loss-of-function mutant
t.test(x= filter(SV1009, Cell == 'QR.pap')$Position, y = filter(KN1500, Cell == 'QR.pap')$Position) #p = 2.039E-5

