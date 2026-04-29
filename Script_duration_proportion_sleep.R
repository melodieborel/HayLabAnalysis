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
install.packages("janitor")
install.packages("magick")
install.packages("patchwork")


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
library(janitor)
library (patchwork)

GLM_L1NDNF_mice_All_VigSt_Global <- read_excel("//10.69.168.1/crnldata/forgetting/Théa/MiniscopeOE_analysis/Exploration_task/1_VigSt_2026-04-15_11_17_10_firstversion/VigStates_Global_noted.xlsx")




####################################
# Duration & proportion & mean duration of Substats per trials 
####################################
matrix <- rbind(GLM_L1NDNF_mice_All_VigSt_Global)
matrix <- matrix [matrix $Substate != "undefined", ]
matrix <- matrix [matrix $Substate != "IS", ]
matrix <- matrix [matrix $ExpeType != "Cheeseboard", ]
matrix <- matrix [matrix $Drug == "baseline", ]
matrix <- matrix[!duplicated(matrix$Substate_ID), ]

result_prop <- matrix %>%
  mutate(DurationSubstate = as.numeric(DurationSubstate)) %>%
  
  group_by(Mice, Day,  ExpeType, Substate) %>%
  summarise(
    Total_Duration = sum(DurationSubstate, na.rm = TRUE),
    Mean_Duration = mean(DurationSubstate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  group_by(Mice, Day, ExpeType) %>%
  mutate(
    Trial_Total = sum(Total_Duration),
    Proportion = Total_Duration / Trial_Total
  ) %>%
  ungroup()

print(result_prop)


##########  Total duration 
#Plot substates duration per mice depending on the day !!!!!!!!!!!!!!!!!!!!!!!!!!!!
plot_data <- result_prop %>%
  mutate(
    Day_num = as.numeric(as.character(Day)),
    
    ExpeType = factor(
      ExpeType,
      levels = c("SleepBefore", "SleepAfter")
    ),
    
    Substate = factor(
      Substate,
      levels = c("AW", "QW", "NREM", "REM")
    ),
    
    # ordre voulu : 1 Before, 1 After, 2 Before, 2 After...
    X_order = factor(
      paste(Day_num, ExpeType, sep = "_"),
      levels = as.vector(rbind(
        paste(sort(unique(Day_num)), "SleepBefore", sep = "_"),
        paste(sort(unique(Day_num)), "SleepAfter",  sep = "_")
      ))
    )
  )

plot_data <- result_prop %>%
  mutate(
    Day_num = as.numeric(as.character(Day)),
    
    ExpeType = factor(
      ExpeType,
      levels = c("SleepBefore", "SleepAfter")
    ),
    
    Substate = factor(
      Substate,
      levels = c("AW", "QW", "NREM", "REM")
    )
  )

levels_order <- c()

for(d in sort(unique(plot_data$Day_num))){
  levels_order <- c(
    levels_order,
    paste(d, "SleepBefore", sep = "_"),
    paste(d, "SleepAfter", sep = "_")
  )
}


plot_data$X_order <- factor(
  paste(plot_data$Day_num, plot_data$ExpeType, sep = "_"),
  levels = levels_order
)

ggplot(plot_data, aes(x = X_order,
                      y = Total_Duration,
                      fill = Substate)) +
  
  geom_bar(stat = "identity", width = 0.8) +
  
  facet_wrap(~Mice, scales = "free_x") +
  
  scale_fill_manual(values = c(
    "AW"   = "orange2",
    "QW"   = "tomato3",
    "NREM" = "royalblue",
    "REM"  = "orchid4"
  )) +
  
  scale_x_discrete(labels = function(x) {
    x <- gsub("_SleepBefore", "\nBefore", x)
    x <- gsub("_SleepAfter", "\nAfter", x)
    x
  }) +
  
  labs(
    x = "Day / Sleep condition",
    y = "Total duration(s)",
    fill = "Substate",
    title = "Substates Duration"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 10),
    strip.text = element_text(face = "bold"),
    strip.background = element_blank()
  )



# Stats des substates duration en fonction des jours !!!!!!!!!!!!!!!!!!!!!!!!!!
data_glm <- result_prop %>%
  mutate(
    Mice = as.factor(Mice),
    Day = as.factor(Day),
    ExpeType = as.factor(ExpeType),
    Substate = as.factor(Substate),
    Total_Duration = as.numeric(Total_Duration)
  )

model <- lmer(Total_Duration ~ Day*Substate + (1|Mice),data = data_glm)
summary(model)
anova(model)
emmeans(model, pairwise ~ Day | Substate)


# Box plot duration substate per day 
result_prop <- result_prop %>%
  mutate(
    Day = factor(Day),
    Substates = factor(Substate),
    ExpeType = factor(ExpeType),
    Mice = factor(Mice)
  )

mouse_colors <- c(
  "NB" = "cornflowerblue",
  "NV" = "darkgoldenrod2",
  "OW" = "yellowgreen",
  "PW" = "indianred3"
)

make_plot <- function(data, expetype_label){
  
  ggplot(data %>% filter(ExpeType == expetype_label),
         aes(x = Day, y = Total_Duration)) +
    
    geom_boxplot(fill = "grey85", color = "grey40",
                 outlier.shape = NA) +
    
    geom_point(aes(color = Mice),
               position = position_jitter(width = 0.15),
               size = 2, alpha = 0.9) +
    
    scale_color_manual(values = mouse_colors) +
    
    facet_wrap(~Substates, ncol = 2) +
    
    labs(
      title = expetype_label,
      x = "Days",
      y = "Substates duration (secondes)",
      color = "Mice"
    ) +
    
    theme_bw() +
    
    theme(
      strip.background = element_blank(),  # ❌ enlève le fond gris
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

p_before <- make_plot(result_prop, "SleepBefore")
p_after  <- make_plot(result_prop, "SleepAfter")

p_before + p_after





#### Mean duration part 

# plot 
plot_data <- result_prop %>%
  mutate(
    Day_num = as.numeric(as.character(Day)),
    
    ExpeType = factor(
      ExpeType,
      levels = c("SleepBefore", "SleepAfter")
    ),
    
    Substate = factor(
      Substate,
      levels = c("AW", "QW", "NREM", "REM")
    ),
    
    # ordre voulu : 1 Before, 1 After, 2 Before, 2 After...
    X_order = factor(
      paste(Day_num, ExpeType, sep = "_"),
      levels = as.vector(rbind(
        paste(sort(unique(Day_num)), "SleepBefore", sep = "_"),
        paste(sort(unique(Day_num)), "SleepAfter",  sep = "_")
      ))
    )
  )

plot_data <- result_prop %>%
  mutate(
    Day_num = as.numeric(as.character(Day)),
    
    ExpeType = factor(
      ExpeType,
      levels = c("SleepBefore", "SleepAfter")
    ),
    
    Substate = factor(
      Substate,
      levels = c("AW", "QW", "NREM", "REM")
    )
  )

levels_order <- c()

for(d in sort(unique(plot_data$Day_num))){
  levels_order <- c(
    levels_order,
    paste(d, "SleepBefore", sep = "_"),
    paste(d, "SleepAfter", sep = "_")
  )
}


plot_data$X_order <- factor(
  paste(plot_data$Day_num, plot_data$ExpeType, sep = "_"),
  levels = levels_order
)

ggplot(plot_data, aes(x = X_order,
                      y = Mean_Duration,
                      fill = Substate)) +
  
  geom_bar(stat = "identity", width = 0.8) +
  
  facet_wrap(~Mice, scales = "free_x") +
  
  scale_fill_manual(values = c(
    "AW"   = "orange2",
    "QW"   = "tomato3",
    "NREM" = "royalblue",
    "REM"  = "orchid4"
  )) +
  
  scale_x_discrete(labels = function(x) {
    x <- gsub("_SleepBefore", "\nBefore", x)
    x <- gsub("_SleepAfter", "\nAfter", x)
    x
  }) +
  
  labs(
    x = "Day / Sleep condition",
    y = "Mean duration(s)",
    fill = "Substate",
    title = "Substates Duration"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 10),
    strip.text = element_text(face = "bold"),
    strip.background = element_blank()
  )

# Stats 
data_glm <- result_prop %>%
  mutate(
    Mice = as.factor(Mice),
    Day = as.factor(Day),
    ExpeType = as.factor(ExpeType),
    Substate = as.factor(Substate),
    Mean_Duration = as.numeric(Mean_Duration)
  )

model <- lmer(Mean_Duration ~ Day*Substate + (1|Mice),data = data_glm)
summary(model)
anova(model)
emmeans(model, pairwise ~ Day | Substate)


# Boxplot 
result_prop <- result_prop %>%
  mutate(
    Day = factor(Day),
    Substates = factor(Substate),
    ExpeType = factor(ExpeType),
    Mice = factor(Mice)
  )

mouse_colors <- c(
  "NB" = "cornflowerblue",
  "NV" = "darkgoldenrod2",
  "OW" = "yellowgreen",
  "PW" = "indianred3"
)

make_plot <- function(data, expetype_label){
  
  ggplot(data %>% filter(ExpeType == expetype_label),
         aes(x = Day, y = Mean_Duration)) +
    
    geom_boxplot(fill = "grey85", color = "grey40",
                 outlier.shape = NA) +
    
    geom_point(aes(color = Mice),
               position = position_jitter(width = 0.15),
               size = 2, alpha = 0.9) +
    
    scale_color_manual(values = mouse_colors) +
    
    facet_wrap(~Substates, ncol = 2) +
    
    labs(
      title = expetype_label,
      x = "Days",
      y = "Mean duration (secondes)",
      color = "Mice"
    ) +
    
    theme_bw() +
    
    theme(
      strip.background = element_blank(),  # ❌ enlève le fond gris
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

p_before <- make_plot(result_prop, "SleepBefore")
p_after  <- make_plot(result_prop, "SleepAfter")

p_before + p_after

# other format box plot
mouse_colors <- c(
  "NB" = "cornflowerblue",
  "NV" = "darkgoldenrod2",
  "OW" = "yellowgreen",
  "PW" = "indianred3"
)

fill_colors <- c("SleepBefore" = "#d0e4f7", "SleepAfter" = "#fde8d0")

jours_ordonnes <- c("1", "2", "3", "4")

result_prop2$DayExpe <- factor(
  result_prop$DayExpe,
  levels = c(
    "1\nSleepBefore", "1\nSleepAfter",
    "2\nSleepBefore", "2\nSleepAfter",
    "3\nSleepBefore", "3\nSleepAfter",
    "4\nSleepBefore", "4\nSleepAfter"
  )
)

make_panel <- function(data, substate_label) {
  
  data_sub <- data %>% filter(Substate == substate_label)
  
  ggplot(data_sub, aes(x = DayExpe, y = Mean_Duration)) +
    
    geom_boxplot(aes(fill = ExpeType),
                 color         = "grey40",
                 outlier.shape = NA,
                 alpha         = 0.7,
                 width         = 0.6) +
    
    geom_point(aes(color = Mice),
               position = position_jitter(width = 0.12),
               size = 2.2, alpha = 0.9) +
    
    scale_fill_manual(values = fill_colors,
                      labels = c("SleepBefore", "SleepAfter"),
                      name   = "Condition") +
    scale_color_manual(values = mouse_colors, name = "Mice") +
    
    # Lignes verticales pour séparer les jours
    geom_vline(xintercept = seq(2.5, by = 2, length.out = length(jours_ordonnes) - 1),
               linetype = "dashed", color = "grey70", linewidth = 0.4) +
    
    labs(
      title = substate_label,
      x     = NULL,
      y     = "Mean duration (s)"
    ) +
    
    theme_bw(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", hjust = 0.5, size = 13),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold"),
      legend.title       = element_text(face = "bold"),
      axis.text.x        = element_text(size = 8),
      panel.grid.major.x = element_blank()
    )
}

p_REM  <- make_panel(result_prop2, "REM")
p_NREM <- make_panel(result_prop2, "NREM")
p_AW   <- make_panel(result_prop2, "AW")
p_QW   <- make_panel(result_prop2, "QW")

(p_REM + p_NREM) / (p_AW + p_QW) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Sleep states mean durations — Before vs After",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
  )



# Vérification
levels(result_prop2$DayExpe)


############## Proportion 

# plot per mice 
plot_data <- result_prop %>%
  mutate(
    Day_num = as.numeric(as.character(Day)),
    
    ExpeType = factor(
      ExpeType,
      levels = c("SleepBefore", "SleepAfter")
    ),
    
    Substate = factor(
      Substate,
      levels = c("AW", "QW", "NREM", "REM")
    ),
    
    # ordre voulu : 1 Before, 1 After, 2 Before, 2 After...
    X_order = factor(
      paste(Day_num, ExpeType, sep = "_"),
      levels = as.vector(rbind(
        paste(sort(unique(Day_num)), "SleepBefore", sep = "_"),
        paste(sort(unique(Day_num)), "SleepAfter",  sep = "_")
      ))
    )
  )

plot_data <- result_prop %>%
  mutate(
    Day_num = as.numeric(as.character(Day)),
    
    ExpeType = factor(
      ExpeType,
      levels = c("SleepBefore", "SleepAfter")
    ),
    
    Substate = factor(
      Substate,
      levels = c("AW", "QW", "NREM", "REM")
    )
  )

levels_order <- c()

for(d in sort(unique(plot_data$Day_num))){
  levels_order <- c(
    levels_order,
    paste(d, "SleepBefore", sep = "_"),
    paste(d, "SleepAfter", sep = "_")
  )
}


plot_data$X_order <- factor(
  paste(plot_data$Day_num, plot_data$ExpeType, sep = "_"),
  levels = levels_order
)

ggplot(plot_data, aes(x = X_order,
                      y = Proportion,
                      fill = Substate)) +
  
  geom_bar(stat = "identity", width = 0.8) +
  
  facet_wrap(~Mice, scales = "free_x") +
  
  scale_fill_manual(values = c(
    "AW"   = "orange2",
    "QW"   = "tomato3",
    "NREM" = "royalblue",
    "REM"  = "orchid4"
  )) +
  
  scale_x_discrete(labels = function(x) {
    x <- gsub("_SleepBefore", "\nBefore", x)
    x <- gsub("_SleepAfter", "\nAfter", x)
    x
  }) +
  
  labs(
    x = "Days / Sleep condition",
    y = "Proportion (ratio)",
    fill = "Substate",
    title = "Substates proportion"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 10),
    strip.text = element_text(face = "bold"),
    strip.background = element_blank()
  )

# all mices
ggplot(plot_data, aes(x = X_order,
                      y = Proportion,
                      fill = Substate)) +
  
  geom_bar(stat = "identity", width = 0.8) +
  
  scale_fill_manual(values = c(
    "AW"   = "orange2",
    "QW"   = "tomato3",
    "NREM" = "royalblue",
    "REM"  = "orchid4"
  )) +
  
  scale_x_discrete(labels = function(x) {
    x <- gsub("_SleepBefore", "\nBefore", x)
    x <- gsub("_SleepAfter", "\nAfter", x)
    x
  }) +
  
  labs(
    x = "Days / Sleep condition",
    y = "Proportion (ratio)",
    fill = "Substate",
    title = "Substates proportion (all mice)"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 10),
    strip.text = element_blank(),
    strip.background = element_blank()
  )

####################################################
# Vigilance states duration
###################################################

matrix <- rbind(GLM_L1NDNF_mice_All_VigSt_Global)
matrix <- matrix [matrix $Substate != "undefined", ]
matrix <- matrix [matrix $Substate != "IS", ]
matrix <- matrix [matrix $ExpeType != "Cheeseboard", ]
matrix <- matrix [matrix $Drug == "baseline", ]
matrix <- matrix[!duplicated(matrix$Substate_ID), ]

matrix$Day <- as.factor(matrix$Day)
matrix$Substate <- as.factor(matrix$Substate)


model<-lmer(DurationSubstate~Substate*Day + (1|Mice)+ (1|ExpeType), data=matrix)
summary(model)
anova(model)
emm <- emmeans(model, pairwise~Substate |Day, lmer.df = "satterthwaite")
emm$contrasts
plot(emm)




####################################################
# Vigilance states proportion
###################################################

matrix <- rbind(GLM_L1NDNF_mice_All_VigSt_Global)
matrix <- matrix [matrix $Substate != "undefined", ]
matrix <- matrix [matrix $Substate != "IS", ]
matrix <- matrix [matrix $ExpeType != "Cheeseboard", ]
matrix <- matrix [matrix $Drug == "baseline", ]
matrix <- matrix[!duplicated(matrix$Substate_ID), ]
matrix$Day <- as.factor(matrix$Day)
matrix$Substate <- as.factor(matrix$Substate)

result2 <- aggregate(DurationSubstate ~ Day + Substate, data = matrix, sum)
result <- aggregate(DurationSubstate ~ Day, data = matrix, sum)
result2$DurTot <- result$DurationSubstate[match(result2$Day, result$Day)]
result2$Mice <- matrix$Mice[match(result2$Day, matrix$Day)]
result2$Prop  <-result2$DurationSubstate/ result2$DurTot *100

model<-lmer(Prop~Substate + (1|Mice) + (1|Day), data=result2)
anova(model)
emm <- emmeans(model, pairwise~Substate, lmer.df = "satterthwaite")
emm$contrasts
plot(emm)

################################
# PLOT
###############################

# dif substates duration 

emm_df <- as.data.frame(emm$emmeans)

p<-ggplot(emm_df, aes(x = Substate, y = emmean)) +
  geom_point(size = 4.5, stroke = 1.2, color = "#2c3e50") +
  #geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, color = "#34495e") +
  geom_errorbar(aes(ymin = emmean-SE, ymax = emmean+SE), width = 0.25, linewidth = 1.4, color = "#007C91") +
  labs(
    title = "secondes",
    x = "",
    y = "Duration of Substates"
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
    y_position = y_base + seq(0.2, by = 2.8, length.out = length(pairs_list)),
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


# dif substate duration depending on day


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
    title = "Seconde",
    x = "",
    y = "Substates duration"
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