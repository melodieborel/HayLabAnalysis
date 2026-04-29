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

GLM_L1NDNF_mice_All_VigSt_Global <- read_excel("//10.69.168.1/crnldata/forgetting/Théa/MiniscopeOE_analysis/Exploration_task/1_VigSt_2026-04-15_11_17_10_firstversion/VigStates_Global_noted.xlsx")

####################################################
# Vigilance states per celltype 
###################################################

matrix <- GLM_L1NDNF_mice_All_VigSt_Global
matrix$UnitNumber <- as.character(matrix$UnitNumber)
matrix$UnitValue <- as.character(matrix$UnitValue)
matrix$UnitNumber <- as.character(matrix$Day)
matrix$UnitValue <- as.character(matrix$Trial)
matrix <- matrix [matrix $Substate != "undefined", ]
matrix <- matrix [matrix $Drug == "baseline", ]
#matrix <- matrix [matrix $Substate != "IS", ]



##### 
matrix <- matrix %>% #only cells that appeared the 4 vigilance states
  group_by(Unit_ID) %>%
  filter(n_distinct(Substate) > 3) %>%
  ungroup()
#####



####################################################
# F(Substate)
###################################################
# Calcium activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(CalciumActivity~Substate + (1|Mice) + (1|Day)+ (1|Trial) + (1|Substate_ID) +(1|Unit_ID) + (1|ExpeType), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~Substate, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Normalized AUC calcium !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(NormalizedAUC_calcium~Substate + (1|Mice) + (1|Day)+ (1|Trial) +(1|Substate_ID) +(1|Unit_ID) + (1|ExpeType), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~Substate, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Deconv Spike Mean Activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(DeconvSpikeMeanActivity~Substate + (1|Mice) + (1|Day)+ (1|Trial) +(1|Substate_ID) +(1|Unit_ID) + (1|ExpeType), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~Substate, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Spike activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(SpikeActivityHz~Substate + (1|Mice) + (1|Day)+ (1|Trial)+ (1|Substate_ID) +(1|Unit_ID) + (1|ExpeType), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~Substate, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)


####################################################
# F(Substate & Day)
###################################################

matrix$Day <- as.factor(matrix$Day)
matrix$Substate <- as.factor(matrix$Substate)

# Calcium activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(CalciumActivity~Day*Substate + (1|Mice) + (1|Substate_ID) +(1|Unit_ID) + (1|ExpeType), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Normalized AUC calcium !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(NormalizedAUC_calcium~Day*Substate + (1|Mice) + (1|Trial) + (1|Substate_ID) +(1|Unit_ID) + (1|ExpeType), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Deconv Spike Mean Activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(DeconvSpikeMeanActivity~Day*Substate + (1|Mice) + (1|Trial) + (1|Substate_ID) +(1|Unit_ID) + (1|ExpeType), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Spike activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(SpikeActivityHz~Day*Substate + (1|Mice) + (1|Trial) + (1|Substate_ID) +(1|Unit_ID) + (1|ExpeType), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)



####################################################
# F(Substate & ExpeType) Without cheesboard session
###################################################

matrix2 <- GLM_L1NDNF_mice_All_VigSt_Global
matrix2 <- matrix2 [matrix2 $ExpeType != "Cheeseboard", ]
matrix2 <- matrix2 [matrix2 $Substate != "undefined", ]
#matrix <- matrix [matrix $Substate != "IS", ]
matrix2 <- matrix2 [matrix2 $Drug == "baseline", ]
matrix$UnitNumber <- as.character(matrix$UnitNumber)
matrix$UnitValue <- as.character(matrix$UnitValue)
matrix$UnitNumber <- as.character(matrix$Day)
matrix$UnitValue <- as.character(matrix$Trial)
matrix$Day <- as.factor(matrix$Day)
matrix$Substate <- as.factor(matrix$Substate)

# Calcium activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(CalciumActivity~ExpeType*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) + (1|Day), data=matrix2)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|ExpeType, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Normalized AUC calcium !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(NormalizedAUC_calcium~ExpeType*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) + (1|Day), data=matrix2)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|ExpeType, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Deconv Spike Mean Activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(DeconvSpikeMeanActivity~ExpeType*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) + (1|Day), data=matrix2)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|ExpeType, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Spike activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(SpikeActivityHz~ExpeType*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) + (1|Day), data=matrix2)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|ExpeType, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)


####################################################
# F(sub & day & expe type) 
###################################################

matrix3 <- GLM_L1NDNF_mice_All_VigSt_Global
matrix3 <- matrix3 [matrix3 $ExpeType == "SleepAfter", ]
matrix3 <- matrix3 [matrix3 $Substate != "undefined", ]
matrix3 <- matrix3 [matrix3 $Substate != "IS", ]
matrix3 <- matrix3 [matrix3 $Drug == "baseline", ]
matrix3$UnitNumber <- as.character(matrix3$UnitNumber)
matrix3$UnitValue <- as.character(matrix3$UnitValue)
matrix3$Day <- as.factor(matrix3$Day)
matrix3$Substate <- as.factor(matrix3$Substate)


# Normalized AUC calcium !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(NormalizedAUC_calcium~Day*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) , data=matrix3)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

####################################################
# F(sub & day) With only sleep after
###################################################

matrix3 <- GLM_L1NDNF_mice_All_VigSt_Global
matri3 <- matrix3 [matrix3 $ExpeType != "Cheeseboard", ]
matri3 <- matrix3 [matrix3 $ExpeType != "SleepBefore", ]
matrix3 <- matrix3 [matrix3 $Substate != "undefined", ]
#matrix <- matrix [matrix $Substate != "IS", ]
matrix3 <- matrix3 [matrix3 $Drug == "baseline", ]
matrix3$UnitNumber <- as.character(matrix3$UnitNumber)
matrix3$UnitValue <- as.character(matrix3$UnitValue)
matrix3$Day <- as.factor(matrix3$Day)
matrix3$Substate <- as.factor(matrix3$Substate)

# Calcium activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(CalciumActivity~Day*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID), data=matrix3)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)


# Normalized AUC calcium !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(NormalizedAUC_calcium~Day*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) + (1|Substate), data=matrix3)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Deconv Spike Mean Activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(DeconvSpikeMeanActivity~Day*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) + (1|Substate), data=matrix3)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Spike activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(SpikeActivityHz~Day*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) + (1|Substate), data=matrix3)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)


####################################################
# F(session day& substate) With only sleep before
###################################################

matrix4 <- GLM_L1NDNF_mice_All_VigSt_Global
matri4 <- matrix4 [matrix4 $ExpeType != "Cheeseboard", ]
matri4 <- matrix4 [matrix4 $ExpeType != "SleepAfter", ]
matrix4 <- matrix4 [matrix4 $Substate != "undefined", ]
#matrix <- matrix [matrix $Substate != "IS", ]
matrix4 <- matrix4 [matrix4 $Drug == "baseline", ]
matrix4$UnitNumber <- as.character(matrix4$UnitNumber)
matrix4$UnitValue <- as.character(matrix4$UnitValue)
matrix4$Day <- as.factor(matrix4$Day)
matrix4$Substate <- as.factor(matrix4$Substate)

# Calcium activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(CalciumActivity~Day*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) , data=matrix4)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Normalized AUC calcium !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(NormalizedAUC_calcium~Day*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) , data=matrix4)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Deconv Spike Mean Activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(DeconvSpikeMeanActivity~Day*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) , data=matrix4)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# Spike activity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(SpikeActivityHz~Day*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) , data=matrix4)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Substate|Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)


matrix$UnitNumber <- as.character(matrix$UnitNumber)
matrix$UnitValue <- as.character(matrix$UnitValue)
model <-lmer(NormalizedAUC_calcium~Substate + (1|Mice) + (1|Session_ID)+ (1|Substate_ID) +(1|Unit_ID), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~Substate, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

matrix <- GLM_L1NDNF_mice_All_VigSt_Global
matrix <- matrix [matrix $Substate != "undefined", ]
matrix <- matrix [matrix $Substate != "IS", ]
matrix <- matrix [matrix $Drug == "baseline", ]
matrix$ClusterHDBSCAN <- as.character(matrix$ClusterHDBSCAN)

matrix <- matrix %>% #only cells that appeared the 4 vigilance states
  group_by(Unit_ID) %>%
  filter(n_distinct(Substate) > 3) %>%
  ungroup()

model <-lmer(NormalizedAUC_calcium~Substate + (1|Mice) + (1|Session_ID)+ (1|Substate_ID) +(1|Unit_ID), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~Substate, lmer.df = "satterthwaite")
#kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)



# NormalizedAUC_calcium !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

matrix <- GLM_L1NDNF_mice_All_VigSt_Global
matrix <- matrix [matrix $Substate != "undefined", ]
matrix <- matrix [matrix $Substate != "IS", ]
matrix <- matrix [matrix $Drug == "baseline", ]
matrix$ClusterHDBSCAN <- as.character(matrix$ClusterHDBSCAN)
model <-lmer(NormalizedAUC_calcium~Substate*ClusterHDBSCAN + (1|Mice) + (1|Session_ID)+ (1|Substate_ID) +(1|Unit_ID), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~ClusterHDBSCAN, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# DeconvSpikeMeanActivity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

matrix <- GLM_L1NDNF_mice_All_VigSt_Global
matrix <- matrix [matrix $Substate != "undefined", ]
matrix <- matrix [matrix $Substate != "IS", ]
matrix <- matrix [matrix $Drug == "baseline", ]
matrix$ClusterHDBSCAN <- as.character(matrix$ClusterHDBSCAN)
model <-lmer(DeconvSpikeMeanActivity~Substate*ClusterHDBSCAN + (1|Mice) + (1|Session_ID)+ (1|Substate_ID) +(1|Unit_ID), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~Substate|ClusterHDBSCAN, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

# EstimatedFiringRate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

matrix <- GLM_L1NDNF_mice_All_VigSt_Global
matrix <- matrix [matrix $Substate != "undefined", ]
matrix <- matrix [matrix $Substate != "IS", ]
matrix <- matrix [matrix $Drug == "baseline", ]
matrix$ClusterHDBSCAN <- as.character(matrix$ClusterHDBSCAN)
model <-lmer(SpikeActivityHz~Substate*ClusterHDBSCAN + (1|Mice) + (1|Session_ID)+ (1|Substate_ID) +(1|Unit_ID), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~Substate|ClusterHDBSCAN, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)


####################################################
# Vigilance states DRUG
###################################################

cellpop='L1NDNF_cells' #L1NDNF_cells L2_3_mice
matrix <- GLM_L1NDNF_mice_All_VigSt_Global 
matrix <- matrix [matrix $Substate != "undefined", ]
matrix <- matrix [matrix $Substate != "IS", ]
matrix$ClusterHDBSCAN <- as.character(matrix$ClusterHDBSCAN)

#matrix <- matrix %>% #only cells that appeared in baseline & CGP
#  group_by(Unit_ID) %>%
#  filter(n_distinct(Drug) > 1) %>%
#  ungroup()
model <-lmer(NormalizedAUC_calcium~Substate*Drug + (1|Mice) + (1|Session_ID) +(1|Unit_ID) +(1|Substate_ID), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~Drug|Substate, lmer.df = "satterthwaite")
#kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)

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

