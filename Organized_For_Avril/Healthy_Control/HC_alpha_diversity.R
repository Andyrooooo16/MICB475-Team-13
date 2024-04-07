### Loading Packages ###

library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(picante)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)



#### Load data ####

load("ibd_filt.RData")
load("ibd_rare.RData")


#### Setting Theme for Plots ####
custom_theme <- theme_minimal() + 
  theme(legend.position = "none", 
        text = element_text(family = "Helvetica"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 14, face = "bold"))

# Set the custom theme as the default theme for all ggplot2 plots
theme_set(custom_theme)







#### Alpha diversity ####
plot_richness(ibd_rare) 

plot_richness(ibd_rare, measures = c("Shannon","Chao1")) 

# Calculate diversity indices

diversity_data <- estimate_richness(ibd_rare, measures = c("Shannon", "Chao1"))
# diversity_data$inflammation_with_surgery[is.na(diversity_data$inflammation_with_surgery)] <- 'Healthy Control'
# ^ James idk what this is here for






#### Inflammation with Surgery Alpha Diversity ####
diversity_data$inflammation_with_surgery <- sample_data(ibd_rare)$inflammation_with_surgery

diversity_data$inflammation_with_surgery <- factor(diversity_data$inflammation_with_surgery, 
                                                   levels = c("Healthy Control", "no_inflammation_no_surgery", "no_inflammation_with_surgery", 
                                                              "inflammation_no_surgery", "inflammation_with_surgery"), ordered = TRUE)

diversity_long_infsurg <- pivot_longer(diversity_data, c("Shannon", "Chao1"), 
                                       names_to = "Measure", values_to = "Value")

p_iws <- ggplot(diversity_long_infsurg, aes(x = inflammation_with_surgery, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#56B4E9", "Chao1" = "#E69F00")) +
  labs(x = "Inflammation with Surgery",
       y = "Diversity Index") +
  facet_wrap(~Measure, scales = "free_y", ncol = 2)

p_iws

compare_means(Shannon ~ inflammation_with_surgery,  data = diversity_data)
compare_means(Chao1 ~ inflammation_with_surgery,  data = diversity_data)

a_my_comparisons <- list( c("Healthy Control", "no_inflammation_with_surgery"), c("no_inflammation_no_surgery", "no_inflammation_with_surgery"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

ggsave("iwf.png", p_iws)



#### Disease Severity Alpha diversity ####

# diversity_data$disease_severity[is.na(diversity_data$disease_severity)] <- 'Healthy Control'
# ^ Again james not sure whats going with this

diversity_data$disease_severity <- factor(sample_data(ibd_rare)$disease_severity, 
                                          levels = c("Healthy Control", "Low", "Medium", "High"), ordered = TRUE)

diversity_long_severity <- pivot_longer(diversity_data, c("Shannon", "Chao1"), 
                                        names_to = "Measure", values_to = "Value")

diversity_long_severity$disease_severity[is.na(diversity_long_severity$disease_severity)] <- 'Healthy Control'

diversity_long_severity$disease_severity <- factor(diversity_long_severity$disease_severity, 
                                                   levels = c("Healthy Control","Low", "Medium", "High"), ordered = TRUE)




p_ds <- ggplot(diversity_long_severity, aes(x = disease_severity, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#56B4E9", "Chao1" = "#E69F00")) +
  labs(x = "Disease Severity",
       y = "Diversity Index") +
  facet_wrap(~Measure, scales = "free_y", ncol = 2)

p_ds

# p_ds +
  # geom_pwc(
  #   aes(group = disease_severity), tip.length = 0,
  #   method = "t_test", label = "p.adj.format",
  #   bracket.nudge.y = 0.1
  # )

ggsave("disease_severity_hc.png", p_ds)




#### CD Location alpha diversity (Not sure if using for project) ####

diversity_data$cd_location <- sample_data(ibd_rare)$cd_location

filtered_diversity_data <- diversity_data %>%
  filter(cd_location != "Ileocolonic and Upper-GI (L3+L4)")

diversity_long_location_filtered <- pivot_longer(filtered_diversity_data, c("Shannon", "Chao1"), names_to = "Measure", values_to = "Value")

p_cdloc <- ggplot(diversity_long_location_filtered, aes(x = cd_location, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#56B4E9", "Chao1" = "#E69F00")) +
  labs(x = "CD Location",
       y = "Diversity Index") +
  facet_wrap(~Measure, scales = "free_y", ncol = 2) 

p_cdloc + 
  geom_pwc(
    aes(group = cd_location), tip.length = 0,
    method = "t_test", label = "p.adj.format",
    bracket.nudge.y = 0.1
  )




# Disease Severity with Inflammation Alpha Diversity
diversity_data$inflammation <- factor(sample_data(ibd_rare)$inflammation) 

diversity_data$DiseaseSev_with_inflammation <- sample_data(ibd_rare)$DiseaseSev_with_inflammation

diversity_long_ds_inf <- pivot_longer(diversity_data, c("Shannon", "Chao1"), 
                                      names_to = "Measure", values_to = "Value")

true_inflammation <- filter(diversity_long_ds_inf, inflammation == TRUE)
false_inflammation <- filter(diversity_long_ds_inf, inflammation == FALSE)

data_Shannon <- subset(true_inflammation, Measure == "Shannon")
data_Chao1 <- subset(true_inflammation, Measure == "Chao1")

data_Shannon$DiseaseSev_with_inflammation <- factor(data_Shannon$DiseaseSev_with_inflammation, 
                                                    levels = c("Healthy Control","Low with inflammation", "Medium with inflammation", "High with inflammation"), ordered = TRUE)
data_Chao1$DiseaseSev_with_inflammation <- factor(data_Chao1$DiseaseSev_with_inflammation, 
                                                  levels = c("Healthy Control","Low with inflammation", "Medium with inflammation", "High with inflammation"), ordered = TRUE)


p_Shannon <- ggplot(data_Shannon, aes(x = DiseaseSev_with_inflammation, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#AED6F1")) +
  labs(x = "Disease Severity and Inflammation", y = "") +
  facet_wrap(~Measure, scales = "free_y", ncol = 2) +
  coord_cartesian(ylim = c(0, 5))

p_Chao1 <- ggplot(data_Chao1, aes(x = DiseaseSev_with_inflammation, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Chao1" = "#F9E79F")) +
  labs(x = "", y = "Diversity Index") +
  facet_wrap(~Measure, scales = "free_y", ncol = 2) +
  coord_cartesian(ylim = c(0, 500))

p_dsinf_true <-  p_Chao1 | p_Shannon

false_inflammation$DiseaseSev_with_inflammation <- factor(false_inflammation$DiseaseSev_with_inflammation, 
                                                          levels = c("Healthy Control","Low no inflammation", "Medium no inflammation", "High no inflammation"), ordered = TRUE)

p_dsinf_false <- ggplot(false_inflammation, aes(x = DiseaseSev_with_inflammation, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#56B4E9", "Chao1" = "#E69F00")) +
  labs(x = "Disease Severity and Inflammation",
       y = "Diversity Index") +
  facet_wrap(~Measure, scales = "free_y", ncol = 4)

p_dsinf_true
p_dsinf_false

ggsave("dsinf_true.png", p_dsinf_true)
ggsave("dsinf_false.png", p_dsinf_false)

# p_dsinf +
#   geom_pwc(
#     aes(group = DiseaseSev_with_inflammation), tip.length = 0,
#     method = "t_test", label = "p.adj.format",
#     bracket.nudge.y = 0.1
#   )
