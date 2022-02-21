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
KN3133 <- read_xlsx("KN3133_smFISH.xlsx")%>%
  add_column(label = "mig-1ΔPromoter 1")

KN3134 <- read_xlsx("KN3134_smFISH.xlsx")%>%
  add_column(label = "mig-1ΔPromoter 2")

SV1009 <- read_xlsx("SV1009_smFISH.xlsx")%>%
  add_column(label = "Control")

KN3222 <- read_xlsx("KN3222_smFISH.xlsx")%>%
  add_column(label = "mig-1ΔPromoter 1-2")


promoterdels <- bind_rows(KN3133,KN3134, KN3222, SV1009)%>%
    filter(Cell %in% c("QR", "QR.p", "QR.pa"))
promoterdels$label <- factor(promoterdels$label, levels = c("Control", "mig-1ΔPromoter 1", "mig-1ΔPromoter 2", "mig-1ΔPromoter 1-2"))

#Create figure panel
fig.3C <- ggplot(promoterdels , aes(x = Position, y = smFISH.Counts, colour = Cell))+
  geom_point(size = 2)+
  scale_x_continuous(labels = c("","V1", "V2", "V3", "V4", "V5"), limits = c(-1,4))+
  scale_y_continuous(limits = c(0,45))+
  geom_vline(xintercept = 1, colour = "black", lty = 2)+
  scale_color_manual(values = c("blue","green3", "magenta"))+
  labs(x = "position along AP axis", y = "# mRNA spots",colour = "Cell")+
  facet_wrap(~label, ncol = 2, scales = "free")+
  theme(legend.position = "none")+
  migtheme
fig.3C


#Statistical testing
# mig-1 expression level in QR of:
#Δpromoter 1
t.test(filter(SV1009, Cell == "QR")$smFISH.Counts,filter(KN3133, Cell == "QR")$smFISH.Counts) #p = 0.5915
#Δpromoter 2
t.test(filter(SV1009, Cell == "QR")$smFISH.Counts,filter(KN3134, Cell == "QR")$smFISH.Counts) #p = 0.06716
#Δpromoter 1-2
t.test(filter(SV1009, Cell == "QR")$smFISH.Counts,filter(KN3222, Cell == "QR")$smFISH.Counts) #p = 0.8895

# mig-1 expression level in QR.pa of:
#Δpromoter 1
t.test(filter(SV1009, Cell == "QR.pa")$smFISH.Counts,filter(KN3133, Cell == "QR.pa")$smFISH.Counts) #p = 0.4381
#Δpromoter 2
t.test(filter(SV1009, Cell == "QR.pa")$smFISH.Counts,filter(KN3134, Cell == "QR.pa")$smFISH.Counts) #p = 2.2e-16
#Δpromoter 1-2
t.test(filter(SV1009, Cell == "QR.pa")$smFISH.Counts,filter(KN3222, Cell == "QR.pa")$smFISH.Counts) #p = 2.2e-16

