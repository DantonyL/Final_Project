library(tidyverse)
base_d <- read_csv("data/revised_data.csv")

# Base Viz

base_viz <- ggplot(data = base_d,
                   aes(x = Mutant,
                       y = `L/P`,
                       color = Mutant)) +
  scale_color_viridis_d() +
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
par(mar=c(5,6,4,1)+.1)

Tukey_plot <- plot(Tukey_mod, las=1, col=my_colors)

# AIC Model Average
library(AICcmodavg)
library(MuMIn)

phyto_lm <- lm(`L/P` ~ Mutant, data = base_d)
phyto_mean <- lm(`L/P` ~ 1, data = base_d)

aictab(list(phyto_lm, phyto_mean), 
       c("linear model", "mean only"),
       second.ord = FALSE)

phyto_lm_2 <- lm(`L/P` ~ poly(Mutant,2), data = base_d)
phyto_lm_3 <- lm(`L/P` ~ poly(Mutant,3), data = base_d)
phyto_lm_4 <- lm(`L/P` ~ poly(Mutant,4), data = base_d)

mod_list <- list(phyto_mean,
                 phyto_lm,
                 phyto_lm_2,
                 phyto_lm_3,
                 phyto_lm_4)

mod_names <- 0:4

aictab(mod_list, mod_names, second.ord = FALSE)
summary(phyto_lm_4)

avg_model <- model.avg(mod_list)
summary(avg_model)

new_data <- data.frame(Mutant = unique(base_d$Mutant))
preds <- predict(avg_model, newdata = new_data, se.fit = TRUE)

pred_df <- data.frame(
  Mutant = new_data$Mutant,
  fit = preds$fit,
  lower = preds$fit - 1.96 * preds$se.fit,
  upper = preds$fit + 1.96 * preds$se.fit)

ggplot(pred_df, aes(x = Mutant, y = fit)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  labs(y = "Predicted L/P", x = "Mutant") +
  theme_minimal()

# Using Lima for Gene Expression Analysis
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
library(limma)

exprs <- read.csv("data/expression_matrix.csv", header=TRUE, row.names=1)
exprs <- exprs[, c("RPKM_PhyB", "RPKM_WT")]
colnames(exprs) <- c("PhyB", "WT")
exprs <- as.matrix(exprs)

metadata <- data.frame(
  Sample = c("PhyB", "WT"),
  Genotype = c("PhyB", "WT"),
  Ratio = c(NA, NA))  # no phenotype data
rownames(metadata) <- metadata$Sample
metadata$Genotype <- factor(metadata$Genotype)
exprs <- log2(exprs + 1)

design <- model.matrix(~ 0 + Genotype, data=metadata)
colnames(design) <- levels(metadata$Genotype)

fit <- lmFit(exprs, design)
contrast.matrix <- makeContrasts(PhyB_vs_WT = PhyB - WT, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)

results <- fit2$coefficients

log2fc <- log2(exprs[, "PhyB"]) - log2(exprs[, "WT"])
results <- data.frame(
  Gene = rownames(exprs),
  log2FoldChange = log2fc
)
head(results[order(-abs(results$log2FoldChange)), ])

library(limma)
A <- (log2(exprs[, "PhyB"]) + log2(exprs[, "WT"])) / 2
M <- log2(exprs[, "PhyB"]) - log2(exprs[, "WT"])
plot(A, M, pch=16, cex=0.6, main="MA Plot of PhyB constrasted to WT", 
     xlab="Average log2 Expression (A)", ylab="log2 Fold Change (M)",
     col=ifelse(abs(M) > 1, "gray", "black"))
abline(h=0, col="red")
highlight_genes <- abs(M) > 2
points(A[highlight_genes], M[highlight_genes], col="skyblue", pchisq=16, cex=1.2)

