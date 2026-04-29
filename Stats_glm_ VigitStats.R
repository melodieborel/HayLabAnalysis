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



install.packages(c(
  "Matrix", "lme4", "lmerTest", "emmeans",
  "readxl", "kableExtra", "sjPlot", "sjmisc", "sjlabelled",
  "dplyr", "stringr", "ggplot2", "tidyr"
))
install.packages("ggsignif")


library(Matrix)
library(ggplot2)      
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(kableExtra)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(ggsignif)

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
matrix <- matrix [matrix $Substate != "IS", ]
matrix3$Day <- as.factor(matrix3$Day)
matrix3$Substate <- as.factor(matrix3$Substate)



##### 
matrix <- matrix %>% #only cells that appeared the 4 vigilance states
  group_by(Unit_ID) %>%
  filter(n_distinct(Substate) > 3) %>%
  ungroup()
#####



####################################################
# F(Substate)
###################################################


# Normalized AUC calcium !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(NormalizedAUC_calcium~Substate + (1|Mice) + (1|Day) +(1|Substate_ID) +(1|Unit_ID) + (1|ExpeType), data=matrix)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model, pairwise~Substate, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)




####################################################
# F(Substate & Day) ExpeType to def 
###################################################

matrix3 <- GLM_L1NDNF_mice_All_VigSt_Global
matrix3 <- matrix3 [matrix3 $ExpeType == "SleepBefore", ]
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
emm <- emmeans(model,pairwise~Day|Substate, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)





####################################################
# F(Substate & SleepB and SleepA)
###################################################

matrix2 <- GLM_L1NDNF_mice_All_VigSt_Global


matrix2 <- matrix2[matrix2$ExpeType != "Cheeseboard", ]
matrix2 <- matrix2[matrix2$Substate != "undefined", ]
matrix2 <- matrix2[matrix2$Substate != "IS", ]
matrix2 <- matrix2[matrix2$Drug == "baseline", ]


matrix2$ExpeType <- factor(as.character(matrix2$ExpeType),
                           levels = c("SleepBefore", "SleepAfter"))

matrix2$Substate <- factor(as.character(matrix2$Substate),
                           levels = c("AW", "NREM", "QW", "REM"))

matrix2$UnitNumber <- as.character(matrix2$UnitNumber)
matrix2$UnitValue  <- as.character(matrix2$UnitValue)
matrix2$UnitNumber <- as.character(matrix2$Day)
matrix2$UnitValue  <- as.character(matrix2$Trial)
matrix2$Day        <- as.factor(matrix2$Day)

levels(matrix2$ExpeType)

# Normalized AUC calcium !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(NormalizedAUC_calcium~ExpeType*Substate + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) + (1|Day), data=matrix2)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~ExpeType|Substate, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)




####################################################
# F(sub & day & expe type) 
###################################################

matrix3 <- GLM_L1NDNF_mice_All_VigSt_Global
matrix3 <- matrix3 [matrix3 $ExpeType == "Cheeseboard", ]
matrix3 <- matrix3 [matrix3 $Substate != "undefined", ]
matrix3 <- matrix3 [matrix3 $Substate != "IS", ]
matrix3 <- matrix3 [matrix3 $Drug == "baseline", ]
matrix3$UnitNumber <- as.character(matrix3$UnitNumber)
matrix3$UnitValue <- as.character(matrix3$UnitValue)
matrix3$Day <- as.factor(matrix3$Day)
matrix3$Substate <- as.factor(matrix3$Substate)


# Normalized AUC calcium !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
model <-lmer(NormalizedAUC_calcium~Day + (1|Mice) + (1|Substate_ID) + (1|Unit_ID) , data=matrix3)
summary(model)
#tab_model(model)
anova(model)
emm <- emmeans(model,pairwise~Day, lmer.df = "satterthwaite")
kable(emm$emmeans)
kable(emm$contrasts)
plot(emm)






###################################
  # PLOT 
####################################


##### Cheeseboard plot 
emm_df <- as.data.frame(emm$emmeans)
contrast_df <- as.data.frame(emm$contrasts)
contrast_df$contrast <- as.character(contrast_df$contrast)

# significance labels
contrast_df$sig_label <- cut(
  contrast_df$p.value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns")
)

sig_contrasts <- contrast_df[contrast_df$sig_label != "ns", ]


p <- ggplot(emm_df, aes(x = Day, y = emmean)) +
  
  geom_point(
    size = 4.5,
    stroke = 1.2,
    color = "#2c3e50"
  ) +
  
  geom_errorbar(
    aes(ymin = emmean - SE,
        ymax = emmean + SE),
    width = 0.25,
    linewidth = 1.4,
    color = "#007C91"
  ) +
  
  
  labs(
    title = "Calcium activity during task across Days",
    x = "",
    y = "Normalized AUC calcium (AU)"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 12
    ),
    axis.text.y = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 12),
    
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(5, "pt"),
    
    panel.grid = element_blank(),
    
    axis.line.x.bottom = element_line(color = "black", linewidth = 1),
    axis.line.y.left   = element_line(color = "black", linewidth = 1)
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black"
  )

if(nrow(sig_contrasts) > 0){
  
  pairs_list <- strsplit(
    as.character(sig_contrasts$contrast),
    " - ",
    fixed = TRUE
  )
  
  sig_contrasts$group1 <- sapply(pairs_list, `[`, 1)
  sig_contrasts$group2 <- sapply(pairs_list, `[`, 2)
  
  # max Y per panel
  ymax_panel <- emm_df %>%
    group_by(Substate) %>%
    summarise(ymax = max(emmean + SE))
  
  sig_contrasts <- left_join(sig_contrasts, ymax_panel, by = "Substate")
  
  sig_contrasts <- sig_contrasts %>%
    group_by(Substate) %>%
    mutate(y_position = ymax + seq(0.15, by = 0.35, length.out = n()))
  
  p_final <- p +
    geom_signif(
      data = sig_contrasts,
      aes(
        xmin = group1,
        xmax = group2,
        annotations = sig_label,
        y_position = y_position
      ),
      manual = TRUE,
      tip_length = 0.02,
      textsize = 5,
      vjust = 0.6
    )
  
} else {
  
  p_final <- p
}

time_stamp <- format(Sys.time(), "%H%M%S")

filename <- paste0(
  "emm_plot_Day_by_Substate_",
  time_stamp,
  ".svg"
)

ggsave(
  filename,
  plot = p_final,
  width = 8,
  height = 3,
  units = "in"
)
print(p_final)



#### Substate 
emm_df <- as.data.frame(emm$emmeans)

p<-ggplot(emm_df, aes(x = Substate, y = emmean)) +
  geom_point(size = 4.5, stroke = 1.2, color = "#2c3e50") +
  #geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "#34495e") +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.25, linewidth = 1.4, color = "#007C91") +
  labs(
    title = "Substates calcium activity",
    x = "",
    y = "Normalized AUC calcium (AU)"
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
  )+ 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")


contrast_df <- as.data.frame(emm$contrasts)
contrast_df$sig_label <- cut(contrast_df$p.value,
                             breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                             labels = c("***", "**", "*", "ns"))
sig_contrasts <- contrast_df[contrast_df$sig_label != "ns", ]


#sig_contrasts <- contrast_df
pairs_list <- strsplit(sig_contrasts$contrast, " - ", fixed = TRUE)
pairs_list <- lapply(pairs_list, unlist)

if (length(pairs_list) > 0) {
  y_base <- max(emm_df$upper.CL, na.rm = TRUE)
  p_final <- p + ggsignif::geom_signif(
    comparisons = pairs_list,
    annotations = sig_contrasts$sig_label,
    y_position = y_base + seq(0.1, by = 0.45, length.out = length(pairs_list)),
    tip_length = 0.02,
    textsize = 5,
    vjust = 0.70
  )
  
} else {
  p_final <- p
}
print(p_final)

time_stamp <- format(Sys.time(), "%H%M%S")
filename <- paste0("emm_plot_VigStates_", time_stamp, ".svg")
ggsave(filename, plot = p_final, width = 2.2, height = 3, units = "in")


##### Substate per day 

library(ggplot2)
library(ggsignif)
library(dplyr)
library(grid)


emm_grid <- emmeans(emm, ~ Day | Substate)
emm_df <- as.data.frame(emm_grid)
emm_contrasts_df <- as.data.frame(contrast(emm_grid, method = "pairwise"))


sig_contrasts <- emm_contrasts_df %>%
  mutate(
    sig_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sig_label))


p <- ggplot(emm_df, aes(x = Day, y = emmean)) +
  
  geom_point(
    size = 4.5,
    stroke = 1.2,
    color = "#2c3e50"
  ) +
  
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.25,
    linewidth = 1.4,
    color = "#007C91"
  ) +
  
  facet_wrap(~Substate) +
  
  labs(
    title = "SleepAfter substates calcium activity ",
    x = "",
    y = "Normalized AUC calcium"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(5, "pt"),
    panel.grid = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 1),
    axis.line.y.left = element_line(color = "black", linewidth = 1)
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30)))


day_levels <- sort(unique(emm_df$Day))

sig_contrasts <- sig_contrasts %>%
  mutate(
    contrast_clean = gsub("[()]", "", contrast),
    split = strsplit(trimws(contrast_clean), "\\s+-\\s+"),
    
    xmin_name = trimws(sapply(split, `[`, 1)),
    xmax_name = trimws(sapply(split, `[`, 2)),
    
    xmin_name_clean = gsub("^Day", "", xmin_name),
    xmax_name_clean = gsub("^Day", "", xmax_name),
    
    xmin = match(xmin_name_clean, day_levels),
    xmax = match(xmax_name_clean, day_levels)
  )


y_range <- diff(range(emm_df$emmean + emm_df$SE, na.rm = TRUE))

y_max_df <- emm_df %>%
  group_by(Substate) %>%
  summarise(
    y_max = max(emmean + SE, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# POSITION DES BARRES
# -----------------------------
sig_contrasts <- sig_contrasts %>%
  group_by(Substate) %>%
  mutate(
    y_position = y_max_df$y_max[match(Substate, y_max_df$Substate)] +
      0.10 * y_range +
      (row_number() - 1) * 0.08 * y_range
  ) %>%
  ungroup()
  
  # Tirets verticaux adaptés à l'échelle
  tick_drop <- 0.02 * y_range
  
  p_final <- p +
    geom_segment(
      data = sig_contrasts,
      aes(x = xmin, xend = xmax,
          y = y_position, yend = y_position),
      linewidth = 0.8
    ) +
    geom_segment(
      data = sig_contrasts,
      aes(x = xmin, xend = xmin,
          y = y_position, yend = y_position - tick_drop),
      linewidth = 0.8
    ) +
    geom_segment(
      data = sig_contrasts,
      aes(x = xmax, xend = xmax,
          y = y_position, yend = y_position - tick_drop),
      linewidth = 0.8
    ) +
    geom_text(
      data = sig_contrasts,
      aes(x = (xmin + xmax) / 2,
          y = y_position + 0.02 * y_range,
          label = sig_label),
      size = 5
    )




time_stamp <- format(Sys.time(), "%H%M%S")

filename <- paste0(
  "emm_plot_VigStates_",
  time_stamp,
  ".svg"
)

ggsave(
  filename,
  plot = p_final,
  width = 2.2,
  height = 3,
  units = "in"
)

print(p_final)




###### SleepB vs SleepA 
emm_grid <- emmeans(model, ~ ExpeType | Substate)
emm_df <- as.data.frame(emm_grid)
emm_contrasts_df <- contrast(
  emm_grid,
  method = "pairwise",
  by = "Substate"
) %>%
  as.data.frame()
emm_contrasts_df <- emm_contrasts_df %>%
  filter(contrast %in% c("sleepafter - sleepbefore",
                         "sleepbefore - sleepafter"))

sig_contrasts <- emm_contrasts_df %>%
  mutate(
    sig_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sig_label))


# -----------------------------
# ORDRE DES SUBSTATES
# -----------------------------
emm_df$Substate <- factor(emm_df$Substate,
                          levels = c("AW", "NREM", "QW", "REM"))

state_levels <- levels(emm_df$Substate)


# -----------------------------
# PLOT PRINCIPAL
# -----------------------------
p <- ggplot(emm_df, aes(x = Substate, y = emmean, color = ExpeType, group = ExpeType)) +
  
  geom_point(
    size = 4.5,
    stroke = 1.2,
    position = position_dodge(width = 0.35),
    color = "#2c3e50"
  ) +
  
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2,
    linewidth = 1.2,
    position = position_dodge(width = 0.35),
    color = "#007C91"
  ) +
  
  labs(
    title = "Substates calcium activity by sleep condition",
    x = "",
    y = "Normalized AUC calcium (AU)"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(5, "pt"),
    panel.grid = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 1),
    axis.line.y.left = element_line(color = "black", linewidth = 1),
    legend.position = "none"
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30)))


# -----------------------------
# CONTRAST CLEANING (optionnel si tu gardes les stats)
# -----------------------------
sig_contrasts <- sig_contrasts %>%
  mutate(
    contrast_clean = gsub("[()]", "", contrast),
    split = strsplit(trimws(contrast_clean), "\\s+-\\s+"),
    
    xmin_name = trimws(sapply(split, `[`, 1)),
    xmax_name = trimws(sapply(split, `[`, 2)),
    
    xmin_name_clean = gsub("^ExpeType", "", xmin_name),
    xmax_name_clean = gsub("^ExpeType", "", xmax_name),
    
    xmin = match(xmin_name_clean, levels(emm_df$ExpeType)),
    xmax = match(xmax_name_clean, levels(emm_df$ExpeType))
  )


# -----------------------------
# Y POSITIONS POUR LES BARRES DE SIGNIFICATIVITÉ
# -----------------------------
y_range <- diff(range(emm_df$emmean + emm_df$SE, na.rm = TRUE))

y_max_df <- emm_df %>%
  group_by(Substate) %>%
  summarise(
    y_max = max(emmean + SE, na.rm = TRUE),
    .groups = "drop"
  )

sig_contrasts <- sig_contrasts %>%
  group_by(Substate) %>%
  mutate(
    y_position = y_max_df$y_max[match(Substate, y_max_df$Substate)] +
      0.10 * y_range +
      (row_number() - 1) * 0.08 * y_range
  ) %>%
  ungroup()

tick_drop <- 0.02 * y_range


# -----------------------------
# AJOUT SIGNIFICATIVITÉ
# -----------------------------
p_final <- p +
  
  geom_segment(
    data = sig_contrasts,
    inherit.aes = FALSE,
    aes(x = xmin, xend = xmax,
        y = y_position, yend = y_position),
    linewidth = 0.8
  ) +
  
  geom_segment(
    data = sig_contrasts,
    inherit.aes = FALSE,
    aes(x = xmin, xend = xmin,
        y = y_position, yend = y_position - tick_drop),
    linewidth = 0.8
  ) +
  
  geom_segment(
    data = sig_contrasts,
    inherit.aes = FALSE,
    aes(x = xmax, xend = xmax,
        y = y_position, yend = y_position - tick_drop),
    linewidth = 0.8
  ) +
  
  geom_text(
    data = sig_contrasts,
    inherit.aes = FALSE,
    aes(x = (xmin + xmax) / 2,
        y = y_position + 0.02 * y_range,
        label = sig_label),
    size = 5
  )

# -----------------------------
# SAUVEGARDE
# -----------------------------
time_stamp <- format(Sys.time(), "%H%M%S")

filename <- paste0(
  "emm_plot_VigStates_",
  time_stamp,
  ".svg"
)

ggsave(
  filename,
  plot = p_final,
  width = 3.2,
  height = 3,
  units = "in"
)

print(p_final)





####### plot between sleep Av et Ap


emm_grid <- emmeans(model, ~ ExpeType | Substate)
emm_df <- as.data.frame(emm_grid)
emm_contrasts_df <- contrast(
  emm_grid,
  method = "pairwise",
  by = "Substate"
) %>%
  as.data.frame()
emm_contrasts_df <- emm_contrasts_df %>%
  filter(contrast %in% c(
    "SleepBefore - SleepAfter",
    "SleepAfter - SleepBefore"
  ))

sig_contrasts <- emm_contrasts_df %>%
  mutate(
    sig_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sig_label))


# -----------------------------
# ORDRE DES SUBSTATES
# -----------------------------
emm_df$Substate <- factor(emm_df$Substate,
                          levels = c("AW", "NREM", "QW", "REM"))


# -----------------------------
# PLOT PRINCIPAL
# -----------------------------
p <- ggplot(emm_df, aes(x = Substate, y = emmean, color = ExpeType, group = ExpeType)) +
  
  geom_point(
    size = 4.5,
    stroke = 1.2,
    position = position_dodge(width = 0.35),
    color = "#2c3e50"
  ) +
  
  geom_errorbar(
    aes(ymin = emmean - SE, ymax = emmean + SE),
    width = 0.2,
    linewidth = 1.2,
    position = position_dodge(width = 0.35),
    color = "#007C91"
  ) +
  
  labs(
    title = "Substates calcium activity by sleep condition",
    x = "",
    y = "Normalized AUC calcium (AU)"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(5, "pt"),
    panel.grid = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 1),
    axis.line.y.left = element_line(color = "black", linewidth = 1),
    legend.position = "none"
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30)))


# -----------------------------
# Y POSITIONS POUR LES BARRES DE SIGNIFICATIVITÉ
# -----------------------------
y_range <- diff(range(emm_df$emmean + emm_df$SE, na.rm = TRUE))

y_max_df <- emm_df %>%
  group_by(Substate) %>%
  summarise(
    y_max = max(emmean + SE, na.rm = TRUE),
    .groups = "drop"
  )
substate_levels <- levels(emm_df$Substate)  # "AW" "NREM" "QW" "REM"
dodge_width <- 0.35

sig_contrasts <- sig_contrasts %>%
  mutate(
    Substate = factor(Substate, levels = substate_levels),
    # Position x centrale du Substate
    x_center = as.numeric(Substate),
    # Les deux points dodgés sont à ± dodge_width/2 autour du centre
    xmin = x_center - dodge_width / 2,
    xmax = x_center + dodge_width / 2,
    y_position = y_max_df$y_max[match(Substate, y_max_df$Substate)] +
      0.12 * y_range
  )

tick_drop <- 0.02 * y_range


# -----------------------------
# AJOUT SIGNIFICATIVITÉ
# -----------------------------
p_final <- p +
  
  geom_segment(
    data = sig_contrasts,
    aes(x = xmin, xend = xmax,
        y = y_position, yend = y_position),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  
  geom_segment(
    data = sig_contrasts,
    aes(x = xmin, xend = xmin,
        y = y_position, yend = y_position - tick_drop),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  
  geom_segment(
    data = sig_contrasts,
    aes(x = xmax, xend = xmax,
        y = y_position, yend = y_position - tick_drop),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  
  geom_text(
    data = sig_contrasts,
    aes(x = (xmin + xmax) / 2,
        y = y_position + 0.02 * y_range,
        label = sig_label),
    inherit.aes = FALSE,
    size = 5
  )

# -----------------------------
# SAUVEGARDE
# -----------------------------
time_stamp <- format(Sys.time(), "%H%M%S")

filename <- paste0(
  "emm_plot_VigStates_",
  time_stamp,
  ".svg"
)

ggsave(
  filename,
  plot = p_final,
  width = 3.2,
  height = 3,
  units = "in"
)

print(p_final)


