install.packages("Matrix")
install.packages("lmerTest")
install.packages("lme4")
install.packages("emmeans")
install.packages("readxl")
install.packages("kableExtra")
install.packages("sjPlot")
install.packages("sjmisc")
install.packages("sjlabelled")
install.packages("dplyr")
install.packages("stringr")
install.packages("ggplot2")
install.packages("readxl")
install.packages("tidyr")

library(Matrix)
library(lmerTest)
library(lme4)
library(emmeans)
library(readxl)
library(kableExtra)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(dplyr)
library(stringr)
library(ggplot2)
library(readxl)
library(tidyr)

Summary_table_L1NDNF <- read_excel("//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_data/L1NDNF_mice/Summary_table_L1NDNF.xlsx")

####################################################
# Sélectionner Session Training
###################################################

matrix <- Summary_table_L1NDNF
matrix <- matrix [matrix $session_type != "Probe", ]
matrix <- matrix [matrix $session_type != "Habituation", ]


# Latence !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
matrix$session <- as.character(matrix$session)
matrix$trial <- as.character(matrix$trial)
model <-lmer(latency_to_reward_s~session + (1|mice) + (1|trial), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~session, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)


# Distance !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(distance_to_reward_cm~session + (1|mice) + (1|trial), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~session, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)


# Speed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(average_speed_cm_s~session + (1|mice), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~session, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)


######### Latency & Distance (withtout trial random factor)
model <-lmer(latency_to_reward_s~session + (1|mice), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~session, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)


model <-lmer(distance_to_reward_cm~session + (1|mice), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~session, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)
#############


####################################################
# Sélectionner Session Prob/Hab
###################################################
matrix2 <- Summary_table_L1NDNF
matrix2 <- matrix2[matrix2$session_type != "Training", ]

# Crossing m  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(crossings_per_m~session_type + (1|mice), data=matrix2)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~session, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Time sec  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(time_spent_in_reward_zone_s~session_type + (1|mice), data=matrix2)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~session, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Time %  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(time_spent_in_reward_zone_percent~session_type + (1|mice), data=matrix2)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~session, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)



############ Graph (SCRIPT TO CORRECT)
emm_df <- as.data.frame(emm$emmeans)
p<-ggplot(emm_df, aes(x = Drug, y = emmean), fill = Substate) +
  geom_point(size = 3, color = "#2c3e50") +
  #geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "#34495e") +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.3, color = "#98695e") +
  facet_wrap(~Substate, nrow=1) +
  labs(
    title = cellpop,
    x = "",
    y = "Calcium activity 0 - 1 s"
  ) +
  #ylim(0,45) +  # Set y-axis limits now
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 14),
    axis.title.x = element_text( size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(5, "pt"),
    panel.border = element_blank(),  # No full border
    panel.grid = element_blank(),    # Clean look, no grid
    axis.line.x.bottom = element_line(color = "black", linewidth = 1),
    axis.line.y.left   = element_line(color = "black", linewidth = 1)
  )+ geom_hline(yintercept = 0, linetype = "dashed", color = "black")


contrast_df <- as.data.frame(emm$contrasts)
contrast_df$sig_label <- cut(contrast_df$p.value,
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", "ns"))
sig_contrasts <- contrast_df[contrast_df$sig_label != "ns", ]
#sig_contrasts <- contrast_df
pairs_list <- str_split(sig_contrasts$contrast, " - ")
pairs_list <- lapply(pairs_list, unlist)
if (length(pairs_list) > 0 && length(sig_contrasts$sig_label) == length(pairs_list)) {
  p_final <- p + ggsignif::geom_signif(
    comparisons = pairs_list,
    annotations = sig_contrasts$sig_label,
    y_position = max(emm_df$upper.CL) + seq(0.1, by = 0.3, length.out = length(pairs_list)),
    tip_length = 0.02,
    textsize = 5,
    vjust = 0.5
  )
} else {
  # No significant contrasts or mismatch, just base plot
  p_final <- p
}
print(p_final)

time_stamp <- format(Sys.time(), "%H%M%S")
filename <- paste0("emm_plot_VigStatesDRUG_", cellpop, "_", time_stamp, ".svg")
ggsave(filename, plot = p_final, width = 2.2, height = 3, units = "in")

