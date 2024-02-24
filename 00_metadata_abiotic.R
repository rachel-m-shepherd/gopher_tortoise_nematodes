### testing soil variables 

library(mctoolsr)
library(tidyverse)
library(emmeans)
library(lme4)
library(lmerTest)
library(effectsize)
library(corrplot)

# input data
tax_table_nem = 'nem_asv_gt.txt'
map_fp = 'metadata_fall_gt.txt'

input_nem = load_taxa_table(tax_table_nem, map_fp)

input_nem$map_loaded$sampleID <- row.names(input_nem$map_loaded)
input_nem = filter_data(input_nem, filter_cat = 'sampleID',
                           filter_vals = c("F.10.O_F"))

map = input_nem$map_loaded

# pH
pH.lm = lmer(pH ~ Habitat * Ecosystem +(1|Site), data = map)
anova(pH.lm)
emmeans(pH.lm, pairwise~Habitat*Ecosystem)

# plot
map$TRT = factor(map$TRT, 
                    levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(map, aes(x = TRT, y = pH)) +
  geom_boxplot(aes(fill = TRT), outlier.shape = NA,
               color = c("black")) +    
  geom_jitter(aes(fill = TRT),
              alpha = 0.8,
              shape = 21,
              size = 3,
              set.seed(666),
              color = "black") +
  scale_fill_manual(values = c("green4","khaki", "green4","khaki")) +
  labs(x = NULL,
       y = "pH",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# OM
OM.lm = lmer(OM ~ Habitat * Ecosystem +(1|Site), data = map)
anova(OM.lm)
emmeans(OM.lm, pairwise~Habitat*Ecosystem)

# plot
map$TRT = factor(map$TRT, 
                 levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(map, aes(x = TRT, y = OM)) +
  geom_boxplot(aes(fill = TRT), outlier.shape = NA,
               color = c("black")) +    
  geom_jitter(aes(fill = TRT),
              alpha = 0.8,
              shape = 21,
              size = 3,
              set.seed(666),
              color = "black") +
  scale_fill_manual(values = c("green4","khaki", "green4","khaki")) +
  labs(x = NULL,
       y = "Organic Matter",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# oC
oC.lm = lmer(oC ~ Habitat * Ecosystem +(1|Site), data = map)
anova(oC.lm)
emmeans(oC.lm, pairwise~Habitat*Ecosystem)

# plot
map$TRT = factor(map$TRT, 
                 levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(map, aes(x = TRT, y = oC)) +
  geom_boxplot(aes(fill = TRT), outlier.shape = NA,
               color = c("black")) +    
  geom_jitter(aes(fill = TRT),
              alpha = 0.8,
              shape = 21,
              size = 3,
              set.seed(666),
              color = "black") +
  scale_fill_manual(values = c("green4","khaki", "green4","khaki")) +
  labs(x = NULL,
       y = "Organic Carbon",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# NO3
NO3.lm = lmer(NO3 ~ Habitat * Ecosystem +(1|Site), data = map)
anova(NO3.lm)
emmeans(NO3.lm, pairwise~Habitat*Ecosystem)

# plot
map$TRT = factor(map$TRT, 
                 levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(map, aes(x = TRT, y = NO3)) +
  geom_boxplot(aes(fill = TRT), outlier.shape = NA,
               color = c("black")) +    
  geom_jitter(aes(fill = TRT),
              alpha = 0.8,
              shape = 21,
              size = 3,
              set.seed(666),
              color = "black") +
  scale_fill_manual(values = c("green4","khaki", "green4","khaki")) +
  labs(x = NULL,
       y = "NO3",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# NH4
NH4.lm = lmer(NH4 ~ Habitat * Ecosystem +(1|Site), data = map)
anova(NH4.lm)
emmeans(NH4.lm, pairwise~Habitat*Ecosystem)

# plot
map$TRT = factor(map$TRT, 
                 levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(map, aes(x = TRT, y = NH4)) +
  geom_boxplot(aes(fill = TRT), outlier.shape = NA,
               color = c("black")) +    
  geom_jitter(aes(fill = TRT),
              alpha = 0.8,
              shape = 21,
              size = 3,
              set.seed(666),
              color = "black") +
  scale_fill_manual(values = c("green4","khaki", "green4","khaki")) +
  labs(x = NULL,
       y = "NH4",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# tP
tP.lm = lmer(tP ~ Habitat * Ecosystem +(1|Site), data = map)
anova(tP.lm)
emmeans(tP.lm, pairwise~Habitat*Ecosystem)

# plot
map$TRT = factor(map$TRT, 
                 levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(map, aes(x = TRT, y = tP)) +
  geom_boxplot(aes(fill = TRT), outlier.shape = NA,
               color = c("black")) +    
  geom_jitter(aes(fill = TRT),
              alpha = 0.8,
              shape = 21,
              size = 3,
              set.seed(666),
              color = "black") +
  scale_fill_manual(values = c("green4","khaki", "green4","khaki")) +
  labs(x = NULL,
       y = "Total Phosphorus",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# H2O
H2O.lm = lmer(H2O ~ Habitat * Ecosystem +(1|Site), data = map)
anova(H2O.lm)
emmeans(H2O.lm, pairwise~Habitat*Ecosystem)

# plot
map$TRT = factor(map$TRT, 
                 levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(map, aes(x = TRT, y = H2O)) +
  geom_boxplot(aes(fill = TRT), outlier.shape = NA,
               color = c("black")) +    
  geom_jitter(aes(fill = TRT),
              alpha = 0.8,
              shape = 21,
              size = 3,
              set.seed(666),
              color = "black") +
  scale_fill_manual(values = c("green4","khaki", "green4","khaki")) +
  labs(x = NULL,
       y = "Soil Moisture",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

## correlations of all metadata

full_met <- read.csv("full_soil_variables.csv")
full_met_filt <- full_met %>%
  filter(Burrow.feature != "Mound") %>%
  filter(Burrow.feature != "Trail")

cormap<-full_met_filt[,6:40]
M=cor(cormap,use="pairwise.complete.obs")
corrplot(M, method = 'number', type = 'upper',number.cex = 0.5,tl.cex = 0.75)

write.csv(M, "full_correlations.csv")
