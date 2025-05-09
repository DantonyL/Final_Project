library(tidyverse)
base_d <- read_csv("data/revised_data.csv")

# Base Viz

base_viz <- ggplot(data = base_d,
                   aes(x = Mutant,
                       y = `L/P`,
                       color = Mutant)) +
  scale_color_viridis_b() +
  geom_boxplot()
base_viz


# Base Model 
library(performance)
library(broom)
library(car)

base_d <- base_d |>
  mutate(Mutant = factor(Mutant, levels = c("WT", setdiff(unique(Mutant), "WT"))))


lm_model <-  lm(`L/P` ~ Mutant, data = base_d)

check_model(lm_model)

tidy(lm_model)
  
# Does model explain variability? - Type 1 SS
anova(lm_model)

# Tukey's Test
Anova_mod <- aov(lm_model)

Tukey_mod <- TukeyHSD(x = Anova_mod, "Mutant", conf.level = 0.95)

my_colors <- c("red", "blue", "green", "orange", "purple", "pink", "cyan", "gold", "turquoise", "coral")

Tukey_plot <- plot(Tukey_mod, las=1, col=my_colors)
