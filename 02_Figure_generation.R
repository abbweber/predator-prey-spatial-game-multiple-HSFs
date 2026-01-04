######---------------------------------------------------------------------#####
### 02: Figure generation ###
######---------------------------------------------------------------------#####

# Creator: Abigail M. Weber
# Correspondence: amw8076@psu.edu

##################
# load libraries #
##################

library(ggplot2)
library(stringr)
library(rcartocolor)

################################################################################
# Step 1: Plot the coefficients from top models (100m buffer) at 90% level #
################################################################################

# Read in the coefficient output
data <- read.csv("output/all_model_coefs_100m_buffer.csv")
data <- data[,-c(1)]

# Drop the intercept values
data <- subset(data, Covariate != "Intercept")

# Change name of models to match the different 'perspectives'
data$Model <- ifelse(data$Model == "Fawn Kill Site - Predator Use (LSDF)", "Predator perspective",
                     ifelse(data$Model == "Fawn Kill Site - Doe Use (LSDF)", "Prey perspective",
                            ifelse(data$Model == "Fawn Kill Site (RSF)", "Landscape perspective", "Predator selection")))

# Subset to be significant covariates for landscape and prey perspectives
data <- subset(data, Covariate == "Exurban" | Covariate == "Other" |
                 Covariate == "Crop" | Covariate == "Distance to Road" |
                 Covariate == "Edge" | Covariate == "Canopy Height")

# Create levels to the models
data$Model <- factor(data$Model, levels = c("Predator selection", "Predator perspective",
                                            "Prey perspective","Landscape perspective"))

# Set a covariate label wrap
data$Covar_wrap <- str_wrap(data$Covariate, width = 10)

# Coefficient plot for significance at 90% level
coefplot <- ggplot(data = subset(data)) +
  geom_point(mapping = aes(x = Covar_wrap, y = Estimate, color = Model, shape = Model),
             position = position_dodge(width = 0.9), size = 5) +
  scale_fill_carto_d("Model", palette = "SunsetDark") +
  scale_shape_manual(name = "Model",
                     breaks = c("Landscape perspective", "Prey perspective",
                                "Predator perspective", "Predator selection"),
                     values = c("Landscape perspective" = 24,
                                "Prey perspective" = 23,
                                "Predator perspective" = 22,
                                "Predator selection" = 21),
                     guide = guide_legend(nrow=2, byrow = T)) +
  geom_errorbar(mapping = aes(x = Covar_wrap, y = Estimate, ymin = Lower_90_CI, ymax = Upper_90_CI,
                              color = Model), width = 0.5,
                position = position_dodge(width = 0.9)) +
  scale_color_manual(name = "Model",
                     breaks = c("Landscape perspective", "Prey perspective",
                                "Predator perspective", "Predator selection"),
                     values = c("Landscape perspective" = "#7C1D6F",
                                "Prey perspective" = "#DC3977",
                                "Predator perspective" = "#F0746E",
                                "Predator selection" = "goldenrod1"),
                     guide = guide_legend(nrow=2, byrow = T)) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  labs(y = "Coefficient Estimates",
       x = NULL) +
  theme_classic() +
  facet_wrap(. ~ Covariate, scale = "free", ncol = 2, nrow = 3,
             labeller = labeller(Covariate = c("Canopy Height" = "(a)",
                                               "Crop" = "(b)",
                                               "Distance to Road" = "(c)",
                                               "Edge" = "(d)",
                                               "Exurban" = "(e)",
                                               "Other" = "(f)"))) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(angle = 90, hjust = 0.5, size = 16, color = "black"),
        axis.title.x = element_text(size = 18, vjust = -1, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 18), legend.text = element_text(size = 18),
        panel.spacing = unit(1, "lines"), strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 14, color = "black"),legend.position = "top") +
  scale_y_continuous(labels = function(x) ifelse(x == 0, "0", scales::label_number()(x))) + # Ensure "0" is consistent
  geom_blank(aes(y = -Estimate, ymin = -Lower_90_CI, ymax = -Upper_90_CI))

plot(coefplot)

# Save the 90% CI plot
ggsave(filename = "figures/Fig_2.jpg",
       plot = coefplot, width = 8.38, height = 7, units = "in")
