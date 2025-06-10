# This file takes output data from other code files and replicates the figures presented in Kraft et al

# load packages
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(gtable)
library(grid)
library(cowplot)
library(ggpubr)
library(ggrepel)

# load the now complete data file
load("analysis_data.RData")


# color palette
cpal <- c("#999999", "#0072B2")


### Figure 2: Barplots of constituent quantities
all$Species <- factor(all$Species, levels=c("orangutan", "gorilla", "chimpanzee", "Human (Hadza)", "Human (Tsimane)"))

cost_bar <- ggplot(all, aes(fill=Sex, y=cost, x=Species)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=low.cost, ymax = hi.cost), stat = "identity", position= position_dodge(width = 0.9), width=0.2) +
  labs(y = expression("Subsistence cost (E"[f]*"), kcal/d")) +
  theme_classic(base_size = 18) +
  scale_fill_manual(values=cpal, labels=c("Female", "Male")) +
  scale_x_discrete(labels= c("orangutan" = "Orangutan", "gorilla" = "Gorilla", "chimpanzee" = "Chimp", "Human (Hadza)" = "Hadza", "Human (Tsimane)" = "Tsimane")) +
  theme(legend.title=element_blank(), legend.position=c(0.2, .80), axis.title.x=element_blank(), axis.text.x=element_text(size=16), axis.text.y=element_text(size=18), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
  
time_bar <- ggplot(all, aes(fill=Sex, y=time_forage_and_moving, x=Species)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=low.time_forage_and_moving, ymax = hi.time_forage_and_moving), stat = "identity", position= position_dodge(width = 0.9), width=0.2) +
  labs(y = expression("Time spent on subsistence (T"[f]*"), hrs")) +
  theme_classic(base_size = 18) +
  scale_fill_manual(values=cpal, labels=c("Female", "Male")) +
  scale_x_discrete(labels= c("orangutan" = "Orangutan", "gorilla" = "Gorilla", "chimpanzee" = "Chimp", "Human (Hadza)" = "Hadza", "Human (Tsimane)" = "Tsimane")) +
  theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_text(size=16), axis.text.y=element_text(size=18))

prod_bar <- ggplot(all, aes(fill=Sex, y=prod, x=Species)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=low.prod, ymax = hi.prod), stat = "identity", position= position_dodge(width = 0.9), width=0.2) +
  labs(y = expression("Production (E"[a]*"), kcal/d")) + 
  theme_classic(base_size = 18) +
  scale_fill_manual(values=cpal, labels=c("Female", "Male")) +
  scale_x_discrete(labels= c("orangutan" = "Orangutan", "gorilla" = "Gorilla", "chimpanzee" = "Chimp", "Human (Hadza)" = "Hadza", "Human (Tsimane)" = "Tsimane")) +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_text(size=16), axis.text.y=element_text(size=18))

Ei_bar <- ggplot(all, aes(fill=Sex, y=E_i, x=Species)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=E_i_low, ymax = E_i_high), stat = "identity", position= position_dodge(width = 0.9), width=0.2) +
  labs(y = expression("Net energy intake (E"[i]*"), kcal/day/kg"^0.75)) +
  theme_classic(base_size = 18) +
  scale_fill_manual(values=cpal, labels=c("Female", "Male")) +
  scale_x_discrete(labels= c("orangutan" = "Orangutan", "gorilla" = "Gorilla", "chimpanzee" = "Chimp", "Human (Hadza)" = "Hadza", "Human (Tsimane)" = "Tsimane")) +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_text(size=16), axis.text.y=element_text(size=18), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

g1 <- ggplotGrob(prod_bar)
g2 <- ggplotGrob(cost_bar)
g3 <- ggplotGrob(time_bar)
g4 <- ggplotGrob(Ei_bar)

gr1 <- do.call(cbind, c(list(g1, g2), size = "first"))
gr2 <- do.call(cbind, c(list(g3, g4), size = "first"))
gall <- do.call(rbind, c(list(gr1,gr2), size="first"))

pdf(file="fig1_barplot.pdf", height=12, width=13)
grid.newpage()
grid.draw(gall)
dev.off()



### Figure 3: Efficiency and return rate
## Prepare cross-cultural data
names(comp)[which(names(comp) == "Ea (kcal/day)")] <- "Ea..kcal."
names(comp)[which(names(comp) == "Ef (kcal/day)")] <- "Ef..kcal."
names(comp)[which(names(comp) == "Tf (hours/day)")] <- "Tf..hours."
names(comp)[which(names(comp) == "Rg (kcal/hr)")] <- "Rg..kcal.hr."

comp$Ea..kcal. <- as.numeric(comp$Ea..kcal.)
comp <- comp[which(comp$use == 1),]
comp$subsistence_mode <- dplyr::recode(comp$subsistence_mode, "forager-horticulturalist" = "horticulturalist")

# rename column
names(comp)[which(names(comp) == "F")] <- "EE"

comp$population <- factor(comp$population, levels = sort(unique(comp$population)))
comp2 <- comp %>% group_by(population, sex, subsistence_mode) %>% summarise_if(is.numeric, mean, na.rm=T)

## Calculate some summaries on the cross cultural data
human_comp_summary <- suppressWarnings(
  comp2 %>% group_by(subsistence_mode, sex) %>% summarise(mu_Rg = mean(Rg..kcal.hr., na.rm=T), mu_F=mean(EE, na.rm=T), sd_Rg = sd(Rg..kcal.hr., na.rm=T), sd_F = sd(EE, na.rm=T), n_Rg=sum(!is.na(Rg..kcal.hr.)), n_F=sum(!is.na(EE)))  %>% 
  mutate(se_Rg = sd_Rg/sqrt(n_Rg), se_F = sd_F/sqrt(n_F)) %>% ungroup()
)

human_comp_summary$Rg_low <- human_comp_summary$mu_Rg - human_comp_summary$se_Rg*1.96
human_comp_summary$Rg_hi <- human_comp_summary$mu_Rg + human_comp_summary$se_Rg*1.96

human_comp_summary$EE_low <- human_comp_summary$mu_F - human_comp_summary$se_F*1.96
human_comp_summary$EE_hi <- human_comp_summary$mu_F + human_comp_summary$se_F*1.96

human_comp_summary$Species <- c(rep("horticulture", 3), rep("HG", 3)) 
human_comp_summary_RR <- human_comp_summary %>% select(Species, sex, mu_Rg, Rg_hi, Rg_low, mu_F, EE_low, EE_hi) %>% rename(Rg = "mu_Rg", EE = "mu_F", Sex = "sex") 
human_comp_summary_RR$Sex[human_comp_summary_RR$Sex == "M"] <- "male"
human_comp_summary_RR$Sex[human_comp_summary_RR$Sex == "F"] <- "female"

human_comp_summary_RR <- bind_rows(all, human_comp_summary_RR)
human_comp_summary_RR$Species <- factor(human_comp_summary_RR$Species, levels= c("orangutan", "gorilla", "chimpanzee", "Human (Hadza)", "HG", "Human (Tsimane)", "horticulture"))
human_comp_summary_RR$Sex <- factor(human_comp_summary_RR$Sex, levels= c("female", "male", "Both"))


cpal2 <- c(cpal, "black")
main_fig1 <- ggplot(all, aes(fill=Sex, y=Rn, x=Species)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=Rn_low, ymax = Rn_hi), stat = "identity", position= position_dodge(width = 0.9), width=0.2) +
  labs(y = expression("Return rate (R"["n"]*"), kcal/hr")) +
  theme_classic(base_size = 18) +
  scale_fill_manual(values=cpal, labels=c("Female", "Male")) +
  scale_x_discrete(labels= c("chimpanzee" = "Chimp", "gorilla" = "Gorilla", "orangutan" = "Orangutan", "Human (Hadza)" = "Hadza", "Human (Tsimane)" = "Tsimane")) +
  theme(legend.title=element_blank(), legend.position=c(0.2, .80), axis.title.x=element_blank(), axis.text.x=element_text(size=14, angle = 45, hjust = 1), axis.text.y=element_text(size=18), axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))+
  ggtitle("(B)")

main_fig2 <- ggplot(all, aes(fill=Sex, y=EE, x=Species)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=EE_low, ymax = EE_hi), stat = "identity", position= position_dodge(width = 0.9), width=0.2) +
  labs(y = expression("Foraging efficiency (F)")) +
  theme_classic(base_size = 18) +
  scale_fill_manual(values=cpal2, labels=c("Female", "Male")) +
  scale_x_discrete(labels= c("chimpanzee" = "Chimp", "gorilla" = "Gorilla", "orangutan" = "Orangutan", "Human (Hadza)" = "Hadza", "Human (Tsimane)" = "Tsimane")) +
  theme(legend.position = "hidden", axis.title.x=element_blank(), axis.text.x=element_text(size=14, angle = 45, hjust = 1), axis.text.y=element_text(size=18), axis.title.y = element_text(margin = margin(t = 0, r = -10, b = 0, l = 0)))+
  ggtitle("(A)")

main_fig <- ggarrange(main_fig2, main_fig1, align="hv")
ggsave(main_fig, filename = "RR_F_net.png", height=6, width=14)
ggsave(main_fig, filename = "RR_F_net.pdf", height=6, width=14)


# Compute averages for table
comp %>% filter(subsistence_mode == "horticulturalist") %>% group_by(population, sex) %>% summarise_if(is.numeric, mean, na.rm=T) %>% group_by(sex) %>% summarise(mu= mean(Ea..kcal., na.rm=T), mu_ef= mean(Ef..kcal., na.rm=T), mu_tf= mean(Tf..hours., na.rm=T), muF= mean(F, na.rm=T), mu_Rg = mean(Rg..kcal.hr., na.rm=T))

comp %>% filter(subsistence_mode == "hunter-gatherer") %>% group_by(population, sex) %>% summarise_if(is.numeric, mean, na.rm=T) %>% group_by(sex) %>% summarise(mu= mean(Ea..kcal., na.rm=T), mu_ef= mean(Ef..kcal., na.rm=T), mu_tf= mean(Tf..hours., na.rm=T), muF= mean(F, na.rm=T), mu_Rg = mean(Rg..kcal.hr., na.rm=T))

comp %>% group_by(subsistence_mode) %>% summarise(n_pops=length(unique(population)))



##### Fig. 7: Lifecourse of efficiency and return rate in Hadza and Tsimane #####
tplot <- ggplot(filter(tsimane, age >= 15), aes(x=age, y=EE, fill=factor(male), color=factor(male))) + geom_line(lwd=1.6) +
  geom_ribbon(aes(ymin=EE_low, ymax=EE_hi), alpha=0.2, col=NA) + 
  labs(x= "Age (years)", y="Foraging efficiency (F)") +
  lims(x=c(15,75), y=c(0,22.5)) +
  theme_classic(base_size=20) +
  scale_fill_manual(name="", labels=c("Female", "Male"), values = cpal) +
  scale_colour_manual(name="", labels=c("Female", "Male"), values = cpal)+
  theme(legend.position = "hidden") +
  ggtitle("(B) Tsimane")


hplot <- ggplot(had, aes(x=age, y=EE, fill=factor(male), color=factor(male))) + geom_line(lwd=1.6) +
  geom_ribbon(aes(ymin=EE_low, ymax=EE_hi), alpha=0.2, col=NA) + 
  labs(x= "Age (years)", y="Foraging efficiency (F)") +
  lims(x=c(15,75), y=c(0,22.5)) +
  theme_classic(base_size=20) +
  scale_fill_manual(name="", labels=c("Female", "Male"), values = cpal) +
  scale_colour_manual(name="", labels=c("Female", "Male"), values = cpal)+
  theme(legend.position = c(0.6, 0.8)) +
  ggtitle("(A) Hadza")


tsimane$male <- factor(tsimane$male)
tplot_RR <- ggplot(filter(tsimane, age >= 15), aes(x=age, y=Rn, fill=factor(male), color=factor(male))) +
  geom_line(lwd=1.6) +
  labs(x= "Age (years)", y=expression("Net return rate (R"["n"]*"), kcal/hr")) +
  geom_ribbon(aes(ymin=Rn_low, ymax=Rn_hi), alpha=0.2, col=NA) + 
  lims(x=c(15,75), y=c(-105, 2050)) +
  theme_classic(base_size=20) +
  scale_fill_manual(name="", labels=c("Female", "Male"), values = cpal) +
  scale_colour_manual(name="", labels=c("Female", "Male"), values = cpal)+
  theme(legend.position = "hidden") +
  ggtitle("(D) Tsimane")

had$male <- factor(had$male)
hplot_RR <- ggplot(filter(had, age > 15), aes(x=age, y=Rn, fill=male, color=male)) + geom_line(lwd=1.6) +
  labs(x= "Age (years)", y=expression("Net return rate (R"["n"]*"), kcal/hr")) +
  geom_ribbon(aes(ymin=Rn_low, ymax=Rn_hi), alpha=0.2, col=NA) + 
  lims(x=c(15,75), y=c(-105, 2050)) +
  theme_classic(base_size=20) +
  scale_fill_manual(name="", labels=c("Female", "Male"), values = cpal) +
  scale_colour_manual(name="", labels=c("Female", "Male"), values = cpal)+
  theme(legend.position = "hidden") +
  ggtitle("(C) Hadza")



png(file="RR_F_lifecourse.png", height=800, width=800)
cowplot::plot_grid(hplot, tplot, hplot_RR, tplot_RR,  align="hv", ncol=2)
dev.off()

pdf(file="RR_F_lifecourse_2021_12_09.pdf", height=10, width=10)
cowplot::plot_grid(hplot, tplot, hplot_RR, tplot_RR,  align="hv", ncol=2)
dev.off()




##### Figure 8 ####
# Note that this must be run in sequence following the figures above (where some processing is done)
comp4 <- comp
comp4$subsistence_mode[which(comp4$subsistence_mode == "horticulturalist")] <- "horticulture"
comp4$Species <- paste("Human (", comp4$subsistence_mode, ")", sep="")
all4 <- all
all4$label <- c("Hadza (male)", "Hadza (female)", "Tsimane (male)", "Tsimane (female)",  rep(NA, 6))
all4$Species <- as.character(all4$Species)
all4$Species[1:2] <- "Human (hunter-gatherer)"
all4$Species[3:4] <- "Human (horticulture)"
all4$Species[5:10] <- "Great ape"

comp4$sex[which(comp4$sex == "M")] <- "male"
comp4$sex[which(comp4$sex == "F")] <- "female"
comp4$sex[which(comp4$sex == "Both")] <- "combined"
comp4$label <- ""

comp4 <- select(comp4, Species, population, sex, Ea..kcal., Tf..hours., Rg..kcal.hr., label)
all4 <- select(all4, Species, Sex, prod, time_forage_and_moving, Rg, label)

# collapse Tsimane estimates to be a combined value like other horticulturalists
newrow <- all4 %>% filter(Species == "Human (horticulture)") %>% summarize(Species = "Human (horticulture)", Sex = "combined", prod = mean(prod), time_forage_and_moving = mean(time_forage_and_moving), Rg = mean(Rg), label="Tsimane")
all4 <- rbind(all4, newrow)
all4[which(all4$label %in% c("Tsimane (male)", "Tsimane (female)")), c("prod", "Rg")] <- NA  # remove non-combined Tsimane estimates for Ea and Rg (we only can compare to combined sexes for these variables)

names(comp4) <- c("species", "population", "sex", "Ea", "Tf", "Rg", "label")
names(all4) <- c("species", "sex", "Ea", "Tf", "Rg","label")

ts <- 2.5
d <- bind_rows(comp4, all4)
d$species <- factor(d$species, levels=c("Great ape", "Human (hunter-gatherer)", "Human (horticulture)"))
d$sex <- factor(d$sex, levels=c("female", "male", "combined"))
d2 <- filter(d, !is.na(Ea), !is.na(species), !is.na(sex)) %>% select(-population, -Tf, -Rg)
ea1 <- ggplot(d2, aes(x=species, y=Ea, fill=sex, col=sex)) + 
        geom_boxplot(width=0.6, position = position_dodge2(width = 0.6, preserve="total"), color="black", alpha=0.2) + 
        geom_point(position = position_dodge(0.6, preserve = 'total')) +
  labs(x="", y=expression("Production (E"[a]*"), kcal/d")) +
  scale_y_continuous(trans='log10') +
  scale_x_discrete(labels=c("Great apes", "Human\n(hunter-gatherer)", "Human\n(horticulture)"))+
  scale_fill_manual(values=c(cpal, "#D55E00")) +
    scale_colour_manual(values=c(cpal, "#D55E00")) +
  # geom_text_repel(data=d2, aes(label = label),
  #                                    size = 3,
  #                                    box.padding = unit(1, "lines"),
  #                                    point.padding = unit(1, "lines")) +
  annotate(
    geom = "segment", x = 2.6, y = 2100, xend = 2.99, yend = 5735.589) +
  annotate(geom = "text", x = 2.5, y = 2000, label = "Tsimane", hjust = "left", size=ts) +
  
  annotate(
    geom = "segment", x = 2.2, y = 6900, xend = 2.01, yend = 5050.883) +
  annotate(geom = "text", x = 2.1, y = 7400, label = "Hadza (male)", hjust = "left", size=ts) +
  
  annotate(
    geom = "segment", x = 2.1, y = 1350, xend = 1.81, yend = 2128.160) +
  annotate(geom = "text", x = 2.01, y = 1300, label = "Hadza (female)", hjust = "left", size=ts) +
  
  theme_classic(base_size = 16)+
  theme(legend.position="none")+
  ggtitle("(A)")



tf1 <- ggplot(d, aes(x=species, y=Tf, fill=sex, col=sex)) + 
  geom_boxplot(width=0.6, position = position_dodge2(width = 0.6, preserve="total"), color="black", alpha=0.2) + 
  geom_point(position = position_dodge(0.6, preserve = 'total')) +
      labs(x="Species/subsistence mode", y=expression("Time spent on subsistence (T"[f]*"), hrs")) +
      scale_x_discrete(labels=c("Great apes", "Human\n(hunter-gatherer)", "Human\n(horticulture)"))+
  scale_fill_manual(values=c(cpal, "#D55E00")) +
  scale_colour_manual(values=c(cpal, "#D55E00")) +
      theme_classic(base_size = 16) +
      theme(legend.position=c(0.4, .90), legend.title=element_blank()) +
  annotate(
    geom = "segment", x = 2.5, y = 2.5, xend = 2.79, yend = 3.938894) +
  annotate(geom = "text", x = 2.2, y = 2.35, label = "Tsimane (female)", hjust = "left", size=ts) +
  
  annotate(
    geom = "segment", x = 3.12, y = 5.8, xend = 3.01, yend = 4.418514) +
  annotate(geom = "text", x = 3.08, y = 6, label = "Tsimane (male)", hjust = "left", size=ts) +
  
  annotate(
    geom = "segment", x = 2.3, y = 8, xend = 2.01, yend = 7.089433) +
  annotate(geom = "text", x = 2.2, y = 8.2, label = "Hadza (male)", hjust = "left", size=ts) +
  
  annotate(
    geom = "segment", x = 1.6, y = 7, xend = 1.79, yend = 5.783077) +
  annotate(geom = "text", x = 1.4, y = 7.2, label = "Hadza (female)", hjust = "left", size=ts) +
  
  theme_classic(base_size = 16)+
  theme(legend.position = "none")+
  ggtitle("(B)")

rg1 <- ggplot(filter(d, !is.na(Rg), Rg < 5000), aes(x=species, y=Rg, fill=sex, col=sex)) + 
  geom_boxplot(width=0.6, position = position_dodge2(width = 0.6, preserve="total"), color="black", alpha=0.2) + 
  geom_point(position = position_dodge(0.6, preserve = 'total')) +
  labs(x="", y=expression("Return rate (R"["g"]*"), kcal/hr")) +
  scale_x_discrete(labels=c("Great apes", "Human\n(hunter-gatherer)", "Human\n(horticulture)"))+
  scale_fill_manual(values=c(cpal, "#D55E00")) +
  scale_colour_manual(values=c(cpal, "#D55E00")) +
  theme_classic(base_size = 16) +
  theme(legend.position=c(.2, .85), legend.title=element_blank()) +
annotate(
    geom = "segment", x = 2.6, y = 1900, xend = 2.99, yend = 1352.5622) +
  annotate(geom = "text", x = 2.35, y = 2050, label = "Tsimane", hjust = "left", size=ts) +
  
  annotate(
    geom = "segment", x = 2.19, y = 250, xend = 2.01, yend = 714.0516) +
  annotate(geom = "text", x = 2.2, y = 200, label = "Hadza (male)", hjust = "left", size=ts) +
  
  annotate(
    geom = "segment", x = 1.6, y = 135, xend = 1.79, yend = 369.9871) +
  annotate(geom = "text", x = 1.25, y = 45, label = "Hadza (female)", hjust = "left", size=ts) +
  ggtitle("(C)")

f4 <- ggarrange(ea1, tf1, rg1, ncol=3)
ggsave(f4, file="fig4_ea_tf_rg_v2.pdf", height=5, width=16)
ggsave(f4, file="fig4_ea_tf_rg_v2.png", height=5, width=16)


##### Figure 6: Percent TEE spent on subsistence ######
# Percentage of TEE spent on subsistence
TEE_bar <- ggplot(all, aes(fill=Sex, y=(cost/TEE)*100, x=Species)) +
  geom_bar(position="dodge", stat="identity") +
  labs(y = "% TEE spent on subsistence", fill = "", x="") +
  theme_classic(base_size = 16) +
  scale_fill_manual(values=cpal, labels=c("Female", "Male")) +
  scale_x_discrete(labels= c("chimpanzee" = "Chimpanzee", "gorilla" = "Gorilla", "orangutan" = "Orangutan", "Human (Hadza)" = "Hadza", "Human (Tsimane)" = "Tsimane")) +
  theme(legend.position = c(0.2, .8))
ggsave(TEE_bar, file="percent_TEE.png", height = 6, width=8)
ggsave(TEE_bar, file="percent_TEE.pdf", height = 6, width=8)


# Figure 4: Stacked bar plot of costs. To load the subplots that get compiled here, see code that generates them within "3_hadza_combine_data.R" and "3_tsimane_combine_data.R"
stacks <- ggpubr::ggarrange(male_stack_hadza, female_stack_hadza, male_cost, female_cost, ncol=2, nrow=2)
ggsave(stacks, file="all_stack_costs.pdf", height=10, width=16)
ggsave(stacks, file="all_stack_costs.png", height=10, width=16)




###### Figure S2: Supplementary figure for energetics of cross-cultural sample ####

comp$sex <- factor(comp$sex, levels=c("M", "F", "Both"))
F_plot <- ggplot(filter(comp, !is.na(EE)), aes(x=population, y=EE, fill=sex)) +
  geom_col(position=position_dodge2(width=0.1, preserve = "single", padding=0), width=.9) +
  labs(y = expression("Foraging efficiency (F)"), fill="") +
  theme_classic(base_size = 18) +
  scale_fill_manual(values=c("#D55E00"), labels=c("Male", "Female", "Combined")) +
  theme(strip.text.x = element_text(size = 15), legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_text(size=16, angle = 90), axis.text.y=element_text(size=18), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

Ef_plot <- ggplot(filter(comp, !is.na(Ef..kcal.)), aes(x=population, y=Ef..kcal., fill=sex)) +
  geom_col(position=position_dodge2(width=0.1, preserve = "single", padding=0), width=.9) +
  labs(y = expression("Subsistence cost (E"[f]*"), kcal/d"), fill="") +
  scale_fill_manual(values=c(cpal, "#D55E00"), labels=c("Male", "Female", "Combined")) +
  theme_classic(base_size = 18) +
  theme(strip.text.x = element_text(size = 15), legend.position=c(.2,.85), axis.title.x=element_blank(), axis.text.x=element_text(size=16, angle = 90), axis.text.y=element_text(size=18), axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))



#legend_Ef <- get_legend(Ef_plot + theme(legend.position="right"))
p1 <- plot_grid(Ef_plot, F_plot,  labels="AUTO", align = "h") 
ggsave("comparative_Ef_F.png", width=15, height = 8)
ggsave("comparative_Ef_F.pdf", width=15, height = 8)



####### Figure S5: Hadza foraging plot #########
library(brms)
load("hadza_lognormal_returns_model.RData")
hadza_age_mean <- 30.94   # these are needed to back transform between age z-scores used in the model and natural scale 
hadza_age_sd <- 18.79
p <- conditional_effects(hadza_lognormal_returns_model, effects=c("age_z:sex"))  # includes the hurdle component
p2 <- conditional_effects(hadza_lognormal_returns_model, effects=c("age_z:sex"), dpar="hu") # to see the hurdle component (Pr zero-days)

h1 <- plot(p)[[1]] + scale_x_continuous(breaks = c(-1, 0, 1, 2), labels= round(c(-1, 0, 1, 2)*hadza_age_sd +hadza_age_mean)) + 
  ggplot2::labs(x="Age", y="Foraging returns (kcal/day)") +
  ggplot2::lims(y=c(0,6500)) + theme_classic() +
  theme_classic(base_size=14) +
  theme(legend.position = c(0.8, 0.8), legend.title = element_blank()) +
  scale_color_manual(values = cpal, labels=c("Female", "Male")) +
  scale_fill_manual(values=cpal, labels=c("Female", "Male")) +
  ggtitle("A")

h2 <- plot(p2)[[1]] + ggplot2::lims(y=c(0,1)) + 
  labs(x="Age", y="Probability of zero day\n(hurdle component of model)") + theme_classic()+
  scale_x_continuous(breaks = c(-1, 0, 1, 2), labels= round(c(-1, 0, 1, 2)*hadza_age_sd + hadza_age_mean)) +
  theme_classic(base_size=14) +
  theme(legend.position = "none") +
  scale_color_manual(values = cpal, labels=c("Female", "Male")) +
  scale_fill_manual(values=cpal, labels=c("Female", "Male"))  +
  ggtitle("B")

h_both <- ggarrange(h1, h2)
ggsave(h_both, filename = "hadza_foraging.png", width=15, height = 8)
ggsave(h_both, filename = "hadza_foraging.pdf", width=15, height = 8)




#####################################################
########## Statistical comparisons in text #########
#####################################################
## Ea
p1 <- comp2 %>% filter( sex == "Both") %>% group_by(population, subsistence_mode) %>%
  summarise(Ea = mean(Ea..kcal., na.rm=T)) %>% filter(!is.na(Ea))

p2 <- comp2 %>% filter( sex != "Both") %>% group_by(population, subsistence_mode) %>%
  summarise(Ea = mean(Ea..kcal., na.rm=T)) %>% filter(!is.na(Ea))

# take preference to use combined if it already existed
p2 <- p2[-which(p2$population %in% p1$population), ]

p3 <- bind_rows(p1, p2)

t.test(p3$Ea[which(p3$subsistence_mode == "horticulturalist")], p3$Ea[which(p3$subsistence_mode == "hunter-gatherer")])

t.test(
  comp2$Ea..kcal.[which(comp2$subsistence_mode == "hunter-gatherer" & comp2$sex == "M")],
  comp2$Ea..kcal.[which(comp2$subsistence_mode == "hunter-gatherer" & comp2$sex == "F")]
)


## Tf
p1 <- comp2 %>% filter( sex == "Both") %>% group_by(population, subsistence_mode) %>%
  summarise(Tf = mean(Tf..hours., na.rm=T)) %>% filter(!is.na(Tf))

p2 <- comp2 %>% filter( sex != "Both") %>% group_by(population, subsistence_mode) %>%
  summarise(Tf = mean(Tf..hours., na.rm=T)) %>% filter(!is.na(Tf))

# take preference to use combined if it already existed
p2 <- p2[-which(p2$population %in% p1$population), ]

p3 <- bind_rows(p1, p2)

# range of times spent
p3 %>% group_by(subsistence_mode) %>% summarise(mean(Tf), range(Tf)[1], range(Tf)[2])

# between subsistence modes
t.test(p3$Tf[which(p3$subsistence_mode == "horticulturalist")], p3$Tf[which(p3$subsistence_mode == "hunter-gatherer")])

# between sexes
t.test(
  comp2$Tf..hours.[which(comp2$subsistence_mode == "hunter-gatherer" & comp2$sex == "M")],
  comp2$Tf..hours.[which(comp2$subsistence_mode == "hunter-gatherer" & comp2$sex == "F")]
)

t.test(
  comp2$Tf..hours.[which(comp2$subsistence_mode == "horticulturalist" & comp2$sex == "M")],
  comp2$Tf..hours.[which(comp2$subsistence_mode == "horticulturalist" & comp2$sex == "F")]
)



## Rg
p1 <- comp2 %>% filter( sex == "Both") %>% group_by(population, subsistence_mode) %>%
  summarise(Rg = mean(Rg..kcal.hr., na.rm=T)) %>% filter(!is.na(Rg))

p2 <- comp2 %>% filter( sex != "Both") %>% group_by(population, subsistence_mode) %>%
  summarise(Rg = mean(Rg..kcal.hr., na.rm=T)) %>% filter(!is.na(Rg))

# take preference to use combined if it already existed
p2 <- p2[-which(p2$population %in% p1$population), ]

p3 <- bind_rows(p1, p2)

# between subsistence modes
t.test(p3$Rg[which(p3$subsistence_mode == "horticulturalist")], p3$Rg[which(p3$subsistence_mode == "hunter-gatherer")])


# between sexes
t.test(
  comp2$Rg..kcal.hr.[which(comp2$subsistence_mode == "hunter-gatherer" & comp2$sex == "M")],
  comp2$Rg..kcal.hr.[which(comp2$subsistence_mode == "hunter-gatherer" & comp2$sex == "F")]
)
