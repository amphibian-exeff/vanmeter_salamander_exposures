dim(rvm_gsh)
summary(rvm_gsh)
colnames(rvm_gsh)

#treatments
unique(rvm_gsh$treatment)
#[1] C   24D CPS
#Levels: 24D C CPS

p1 <- rvm_gsh %>%
  mutate(treatment = fct_relevel(treatment, 
                            "C", "24D", "CPS")) %>%
    ggplot(aes(x=treatment, y=weight_g, fill=treatment)) +
      geom_boxplot() +
      geom_jitter(shape=16, position=position_jitter(0.2), color="darkorchid1") +
      xlab("Treatment") + ylab("body weight (g)") + ggtitle("Salamander Body Weight") +
      theme_bw()
p1

p2 <- rvm_gsh %>%
  mutate(treatment = fct_relevel(treatment, 
                                 "C", "24D", "CPS")) %>%
    ggplot(aes(x=treatment, y=svl_mm, fill=treatment)) +
      geom_boxplot(show.legend = FALSE) +
      geom_jitter(shape=16, position=position_jitter(0.2), show.legend = FALSE, color="darkorchid1") +
      xlab("Treatment") + ylab("snout-vent-length (mm)") + ggtitle("Salamander Snout-Vent-Length") +
      theme_bw()
p2

#combined figure jpg of body weights and svl boxplots
metrics_combined <- ggarrange(p1, p2, heights = c(4, 4), widths=c(3.7,3.3),
                              labels = c("A", "B"),
                              ncol = 2, nrow = 1)
metrics_combined

metrics_figure_filename <- paste(rvm_graphics,"/rvm_salamander_gsh_metrics_figure.jpg",sep="")
jpeg(metrics_figure_filename, width = 7, height = 4, units = "in",res=600)
  metrics_combined
dev.off()

# gsh alternatives for the dilutions, 1:5, 1:8 and average of teh 2 (we use the average)
# average of 1:5 and 1:8 for gsh
p3 <- rvm_gsh %>%
  mutate(treatment = fct_relevel(treatment, 
                                 "C", "24D", "CPS")) %>%
  ggplot(aes(x=treatment, y=gsh_nM_mL, fill=treatment)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), color="darkorchid1") +
  xlab("Treatment") + ylab("Glutathione (nM_mL)") + ggtitle("Glutathione") +
  ylim(0,38) +
  theme_bw() + theme(legend.position = "none")
p3

#1:5 gsh
p4 <- rvm_gsh %>%
  mutate(treatment = fct_relevel(treatment, 
                                 "C", "24D", "CPS")) %>%
    ggplot(aes(x=treatment, y=gsh_1_5_dilution_nM_mL, fill=treatment)) +
      geom_boxplot() +
      geom_jitter(shape=16, position=position_jitter(0.2), color="darkorchid1") +
      xlab("Treatment") + ylab("Glutathione (1:5) (nM_mL)") + ggtitle("Glutathione (1:5)") +
      ylim(0,38) +
      theme_bw() + theme(legend.position = "none")
p4

#1:8 gsh
p5 <- rvm_gsh %>%
  mutate(treatment = fct_relevel(treatment, 
                                 "C", "24D", "CPS")) %>%
    ggplot(aes(x=treatment, y=gsh_1_8_dilution_nM_mL, fill=treatment)) +
      geom_boxplot() +
      geom_jitter(shape=16, position=position_jitter(0.2), color="darkorchid1") +
      xlab("Treatment") + ylab("gsh_1_8_dilution_nM_mL") + ggtitle("Glutathion (1:8)") +
      ylim(0,38) +
      theme_bw()
p5

#combined gsh figure jpg
gsh_combined <- ggarrange(p3, p4, p5, heights = c(4,4,4), widths=c(3.3,3.3,4.8),
                              labels = c("A", "B", "C"),
                              ncol = 3, nrow = 1)
gsh_combined

gsh_figure_filename <- paste(rvm_graphics,"/rvm_salamander_gsh_figure.jpg",sep="")
jpeg(gsh_figure_filename, width = 7, height = 4, units = "in",res=600)
  gsh_combined
dev.off()

#acetylcholinesterase
p6 <- rvm_gsh %>%
  mutate(treatment = fct_relevel(treatment, 
                                 "C", "24D", "CPS")) %>%
  ggplot(aes(x=treatment, y=ache_ug_min_mg, fill=treatment)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), color="darkorchid1") +
  xlab("Treatment") + ylab("Acetylcholinesterase (ug_min_mg)") + ggtitle("Acetylcholinesterase") +
  theme_bw()
p6

# glutathione swabs
colnames(rvm_gsh)
colnames(rvm_gsh_swabs)
unique(rvm_gsh_swabs$treatment)
max(rvm_gsh_swabs$total_GSH)
p7 <- rvm_gsh_swabs %>%
  mutate(treatment = fct_relevel(treatment, 
                                 "CON", "D", "CHL")) %>%
  ggplot(aes(x=treatment, y=total_GSH, fill=treatment)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), color="darkorchid1") +
  xlab("Treatment") + ylab("Glutathione (nM_mL) Swabs") + ggtitle("Glutathione") +
  ylim(0,0.3) +
  theme_bw() + theme(legend.position = "none")
p7

#combined figure jpg
responses_combined <- ggarrange(p6, p3, heights = c(4, 4), widths=c(3.7,2.7),
                          labels = c("A", "B"),
                          ncol = 2, nrow = 1)
responses_combined

# combined figure
response_figure_filename <- paste(rvm_graphics,"/rvm_salamander_response_figure.jpg",sep="")
jpeg(response_figure_filename, width = 7, height = 4, units = "in",res=600)
 responses_combined
dev.off()
