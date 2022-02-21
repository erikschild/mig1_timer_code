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
AIDctrl <- read_xlsx("AIDctrl_smFISH.xlsx")%>%
  add_column(label = "AID::CDK-1, no auxin")

AID <- read_xlsx("AID_smFISH.xlsx")%>%
  filter(Cell != "QR.pa")

AID$label <- "AID::CDK-1 + auxin"

SV <- read_xlsx("SV1009_smFISH.xlsx")

all_AID_smFISH <- bind_rows(AIDctrl, AID)%>%
  filter(Cell != "unknown")
all_AID_smFISH$label <- factor(all_AID_smFISH$label, levels = c("AID::CDK-1, no auxin", "AID::CDK-1 + auxin"))

AID_DAPI <- read_xlsx("AID_DAPIquant.xlsx")%>%
  add_column(label = "AID::CDK-1")

#create figure panels
fig.1B <- ggplot(filter(SV, Cell %in% c("QR", "QR.p", "QR.pa", "QR.pp")) , aes(x = Position, y = smFISH.Counts, colour = Cell))+
  geom_point(size = 2)+
  scale_x_continuous(labels = c("H2","V1", "V2", "V3", "V4", "V5"), limits = c(-1,4))+
  scale_y_continuous(limits = c(0,45))+
  geom_vline(xintercept = 1, colour = "black", lty = 2)+
  scale_color_manual(values = c("blue", "green3", "magenta", "darkorange1"))+
  labs(x = "position along AP axis", y = "# mRNA spots",colour = "Cell")+
  migtheme
fig.1B


fig.1C <- ggplot(filter(AID_DAPI, Cell == "QR.p*" | Cell == 'seam' ), aes(x = Cell,y = Ratio))+
  xlab(label = "")+
  ylab("DNA content")+
  geom_boxplot(colour = "black", fill = "Gray80", size = 0.3)+
  geom_jitter(colour = "black", width = 0.06, height = 0, size = 0.1)+
  scale_y_continuous(limits = c(0.5,4),labels = c("","2n", "4n", "", "8n"))+
  migtheme
fig.1C

fig.1E <- ggplot(filter(all_AID_smFISH, Cell %in% c("QR.p", "QR.pa")) , aes(x = Position, y = smFISH.Counts, colour = Cell))+
  geom_point(size = 2)+
  scale_x_continuous(labels = c("H2","V1", "V2", "V3", "V4", "V5"), limits = c(-1,4))+
  scale_y_continuous(limits = c(0,45))+
  geom_vline(xintercept = 1, colour = "black", lty = 2)+
  scale_color_manual(values = c("green3", "magenta"))+
  labs(x = "position along AP axis", y = "# mRNA spots",colour = "Cell")+
  facet_wrap(~label, ncol = 2, scales = "free")+
  migtheme
fig.1E

#statistical testing
#test if DNA content AID QR.p is less than 8n
t.test(filter(AID_DAPI, Cell == "QR.p*")$Ratio, mu = 4, alternative = "less") # p = 2.2e-16, yes
#test if DNA contentAID QR.p is greater than 4n
t.test(filter(AID_DAPI, Cell == "QR.p*")$Ratio, mu = 2, alternative = "greater") # p= 0.696, no

