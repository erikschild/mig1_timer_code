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
#install.packages("tidyr")
library(tidyr)

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
KN2651 <- read_excel("KN2651_smFISH.xlsx")%>%
  add_column(label = "bar-1 gof")

KN1739 <- read_excel("KN1739_smFISH.xlsx")%>%
  add_column(label = "bar-1 lof")

SV1009 <- read_excel("SV1009_smFISH.xlsx")%>%
  add_column(label = "Control")%>%
  filter(Date %in% c("20210415", "20190813", "20191015")) #To ensure precise comparison for this figure, we only use the wild type data which was staged in parallel with the bar-1 strains

variances <- read_xlsx("Variances.xlsx")%>%
  slice(1)%>% #we plot only the most inclusive case, i.e. minimum 10 mig-1 expressed
  pivot_longer(cols = c(2:4))
variances$name <- factor(variances$name, levels = c("Control", "bar-1 lof", "bar-1 gof"))

all_bar1 <- bind_rows(KN2651, KN1739, SV1009)
all_bar1$label <- factor(all_bar1$label, levels = c("Control", "bar-1 lof", "bar-1 gof"))


#Create figure panels
fig.5A <- ggplot(filter(all_bar1, Cell %in% c("QR.p", "QR.pa")) , aes(x = Position, y = smFISH.Counts, colour = Cell))+
  geom_point(size = 2)+
  scale_x_continuous(labels = c("H2","V1", "V2", "V3", "V4", "V5"), limits = c(-1,4))+
  scale_y_continuous(limits = c(0,45))+
  geom_vline(xintercept = 1, colour = "black", lty = 2)+
  scale_color_manual(values = c("green3", "magenta", "darkorange1"))+
  labs(x = "position along AP axis", y = "# mRNA spots",colour = "Cell")+
  facet_wrap(~label, ncol = 2, scales = "free")+
  migtheme
fig.5A

fig.5B <- ggplot(variances, aes(x = name, y = value))+
  geom_col(fill = "Gray80", color = "Gray80")+
  scale_y_continuous(limits = c(0,0.45))+
  ylab("Variance")+
  xlab("")+
  migtheme
fig.5B

fig.5sup1 <- ggplot(filter(all_bar1, Cell == "QR.pa") , aes(x = label, y = smFISH.Counts))+
  xlab(label = "")+
  ylab("# mRNA spots")+
  geom_boxplot(colour = "black", fill = "Gray80", size = 0.3)+
  stat_boxplot(geom = 'errorbar')+
  geom_jitter(colour = "black", width = 0.06, height = 0, size = 0.1)+
  migtheme
fig.5sup1

#statistical testing
# mig-1 expression level in QR of:
#bar-1 gain-of-function
t.test(filter(KN2651, Cell == "QR.pa")$smFISH.Counts,filter(SV1009, Cell == "QR.pa")$smFISH.Counts)#p = 0.8176
#bar-1 loss-of-function
t.test(filter(KN1739, Cell == "QR.pa")$smFISH.Counts,filter(SV1009,Cell == "QR.pa")$smFISH.Counts)#p = 0.01789

