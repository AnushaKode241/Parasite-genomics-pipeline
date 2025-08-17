
#All the boxplots in this paper are generated using the below basic R script:

# load libraries
library(viridis)
library(ggplot2)
library(ggsignif)
library(ggpubr)

## load the file
library(readxl)

WGA_comparison <- read_excel("path_to_inputfile.xlsx", sheet = "Sheet1")
View(WGA_comparison)

## here is one way to do color fill in ggplot: add an aesthetic to geom_point and indicate the color using scale_fill_manual
p <- ggplot(WGA_comparison,aes(x=WGA_Kits,y=concn_WGA))+
  geom_point(aes(x=WGA_Kits,y=concn_WGA,fill=Person_ID),shape=21,size=3,position=position_jitter(width=0.1, height = 0))+
  geom_boxplot(fill=NA,outlier.color=NA)+scale_fill_manual(values=c("#b3e2cd","#fdcdac","#cbd5e8", "#f4cae4","#e6f5c9", "#fff2ae")) +
  stat_compare_means(method = "anova", label.y = 100.5) +
  stat_compare_means(label = "p.signif", method = "anova", label.x = 1.53) 
  #stat_compare_means(label = "p.signif",method = "t.test", label.y = 0.65, ref.group = ".all.")
 
p + labs(x = "WGA kits", y = "DNA yield post WGA", fill = "Person ID")

