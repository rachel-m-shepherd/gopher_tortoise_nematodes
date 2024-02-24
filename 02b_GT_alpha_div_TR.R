# abundance of Nematode families

library(mctoolsr)
library(tidyverse)
library(emmeans)
library(lme4)
library(lmerTest)
library(effectsize)

# input data
tax_table_nem = 'nem_asv_gt.txt'
map_fp = 'metadata_fall_gt.txt'

input_nem = load_taxa_table(tax_table_nem, map_fp)

input_nem_ra = convert_to_relative_abundances(input_nem)
input_nem_ra$map_loaded$sampleID <- row.names(input_nem_ra$map_loaded)
input_nem_ra = filter_data(input_nem_ra, filter_cat = 'sampleID',
                           filter_vals = c("F.10.O_F"))

# calc nem richness in each sample for major groups
# BF
BF_alpha = filter_taxa_from_input(input_nem_ra, at_spec_level = 3,
                                  taxa_to_keep = "BF")
BF_rich = calc_diversity(BF_alpha$data_loaded, metric = 'richness')
BF_div = calc_diversity(BF_alpha$data_loaded, metric = 'shannon')

sort(BF_rich)
sort(BF_div)


# join to map
BF_map = BF_alpha$map_loaded
BF_map$BF_rich = BF_rich
BF_map$BF_div = BF_div

BF.rich.lm = lmer(BF_rich ~ Habitat * Ecosystem +(1|Site), data = BF_map)
anova(BF.rich.lm)
emmeans(BF.rich.lm, pairwise~Habitat*Ecosystem)


BF.div.lm = lmer(BF_div ~ Habitat * Ecosystem +(1|Site), data = BF_map)
anova(BF.div.lm)
emmeans(BF.div.lm, pairwise~Habitat*Ecosystem)

# effect size
BF_map$Ecosystem = factor(BF_map$Ecosystem, 
                           levels = c("Native", "Degraded"))

cohens_d(BF_rich ~ Ecosystem, data = BF_map)
cohens_d(BF_div ~ Ecosystem, data = BF_map)

# plot
BF_map$TRT = factor(BF_map$TRT, 
                     levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(BF_map, aes(x = TRT, y = BF_rich)) +
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
       y = "Richness",
       title = "Bacterial-Feeding") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(BF_map, aes(x = TRT, y = BF_div)) +
  geom_boxplot(aes(fill = TRT), outlier.shape = NA,
               color = c("black")) +    
  geom_jitter(aes(fill = TRT),
              alpha = 0.8,
              shape = 21,
              size = 3,
              set.seed(666),
              color = "black") +
  scale_fill_manual(values = c("green4","green4", "khaki","khaki")) +
  labs(x = NULL,
       y = "Shannon Diversity",
       title = "Bacterial-Feeding") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# FF
FF_alpha = filter_taxa_from_input(input_nem_ra, at_spec_level = 3,
                                  taxa_to_keep = "FF")
FF_rich = calc_diversity(FF_alpha$data_loaded, metric = 'richness')
FF_div = calc_diversity(FF_alpha$data_loaded, metric = 'shannon')

sort(FF_rich)
sort(FF_div)


# join to map
FF_map = FF_alpha$map_loaded
FF_map$FF_rich = FF_rich
FF_map$FF_div = FF_div

FF.rich.lm = lmer(FF_rich ~ Habitat * Ecosystem +(1|Site), data = FF_map)
anova(FF.rich.lm)
emmeans(FF.rich.lm, pairwise~Habitat*Ecosystem)


FF.div.lm = lmer(FF_div ~ Habitat * Ecosystem +(1|Site), data = FF_map)
anova(FF.div.lm)
emmeans(FF.div.lm, pairwise~Habitat*Ecosystem)

# effect size
FF_map$Ecosystem = factor(FF_map$Ecosystem, 
                          levels = c("Native", "Degraded"))

cohens_d(FF_rich ~ Ecosystem, data = FF_map)
cohens_d(FF_div ~ Ecosystem, data = FF_map)

# plot
FF_map$TRT = factor(FF_map$TRT, 
                    levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(FF_map, aes(x = TRT, y = FF_rich)) +
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
       y = "Richness",
       title = "Fungal-Feeding") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(FF_map, aes(x = TRT, y = FF_div)) +
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
       y = "Shannon Diversity",
       title = "Fungal-Feeding") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# PP
PP_alpha = filter_taxa_from_input(input_nem_ra, at_spec_level = 3,
                                  taxa_to_keep = "PP")
PP_rich = calc_diversity(PP_alpha$data_loaded, metric = 'richness')
PP_div = calc_diversity(PP_alpha$data_loaded, metric = 'shannon')

sort(PP_rich)
sort(PP_div)


# join to map
PP_map = PP_alpha$map_loaded
PP_map$PP_rich = PP_rich
PP_map$PP_div = PP_div

PP.rich.lm = lmer(PP_rich ~ Habitat * Ecosystem +(1|Site), data = PP_map)
anova(PP.rich.lm)
emmeans(PP.rich.lm, pairwise~Habitat*Ecosystem)


PP.div.lm = lmer(PP_div ~ Habitat * Ecosystem +(1|Site), data = PP_map)
anova(PP.div.lm)
emmeans(PP.div.lm, pairwise~Habitat*Ecosystem)

# effect size
PP_map$Ecosystem = factor(PP_map$Ecosystem, 
                          levels = c("Native", "Degraded"))

cohens_d(PP_rich ~ Ecosystem, data = PP_map)
cohens_d(PP_div ~ Ecosystem, data = PP_map)

# plot
PP_map$TRT = factor(PP_map$TRT, 
                    levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(PP_map, aes(x = TRT, y = PP_rich)) +
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
       y = "Richness",
       title = "Plant-Parasitic") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(PP_map, aes(x = TRT, y = PP_div)) +
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
       y = "Shannon Diversity",
       title = "Plant-Parasitic") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# OM
OM_alpha = filter_taxa_from_input(input_nem_ra, at_spec_level = 3,
                                  taxa_to_keep = "OM")
OM_rich = calc_diversity(OM_alpha$data_loaded, metric = 'richness')
OM_div = calc_diversity(OM_alpha$data_loaded, metric = 'shannon')

sort(OM_rich)
sort(OM_div)


# join to map
OM_map = OM_alpha$map_loaded
OM_map$OM_rich = OM_rich
OM_map$OM_div = OM_div

OM.rich.lm = lmer(OM_rich ~ Habitat * Ecosystem +(1|Site), data = OM_map)
anova(OM.rich.lm)
emmeans(OM.rich.lm, pairwise~Habitat*Ecosystem)


OM.div.lm = lmer(OM_div ~ Habitat * Ecosystem +(1|Site), data = OM_map)
anova(OM.div.lm)
emmeans(OM.div.lm, pairwise~Habitat*Ecosystem)

# effect size
OM_map$Ecosystem = factor(OM_map$Ecosystem, 
                          levels = c("Native", "Degraded"))

cohens_d(OM_rich ~ Ecosystem, data = OM_map)
cohens_d(OM_div ~ Ecosystem, data = OM_map)

# plot
OM_map$TRT = factor(OM_map$TRT, 
                    levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(OM_map, aes(x = TRT, y = OM_rich)) +
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
       y = "Richness",
       title = "Omnivorous") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(OM_map, aes(x = TRT, y = OM_div)) +
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
       y = "Shannon Diversity",
       title = "Omnivorous") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# RA
RA_alpha = filter_taxa_from_input(input_nem_ra, at_spec_level = 3,
                                  taxa_to_keep = "RA")
RA_rich = calc_diversity(RA_alpha$data_loaded, metric = 'richness')
RA_div = calc_diversity(RA_alpha$data_loaded, metric = 'shannon')

sort(RA_rich)
sort(RA_div)


# join to map
RA_map = RA_alpha$map_loaded
RA_map$RA_rich = RA_rich
RA_map$RA_div = RA_div

RA.rich.lm = lmer(RA_rich ~ Habitat * Ecosystem +(1|Site), data = RA_map)
anova(RA.rich.lm)
emmeans(RA.rich.lm, pairwise~Habitat*Ecosystem)


RA.div.lm = lmer(RA_div ~ Habitat * Ecosystem +(1|Site), data = RA_map)
anova(RA.div.lm)
emmeans(RA.div.lm, pairwise~Habitat*Ecosystem)

# effect size
RA_map$Ecosystem = factor(RA_map$Ecosystem, 
                          levels = c("Native", "Degraded"))

cohens_d(RA_rich ~ Ecosystem, data = RA_map)
cohens_d(RA_div ~ Ecosystem, data = RA_map)

# plot
RA_map$TRT = factor(RA_map$TRT, 
                    levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(RA_map, aes(x = TRT, y = RA_rich)) +
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
       y = "Richness",
       title = "Root-Associated") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))



ggplot(RA_map, aes(x = TRT, y = RA_div)) +
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
       y = "Shannon Diversity",
       title = "Root-Associated") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))


###### PREDICTORS
### abiotic/biotic predictors
## BF
deg_ap = BF_map %>%
  filter(TRT == "DegradedApron")

# relate to all abiotic variables
deg_ap.lm = lm(BF_rich ~ pH + OM + oC + NO3 + NH4 + tP + H2O + BactRich, data = deg_ap)
summary(deg_ap.lm)

deg_ar = BF_map %>%
  filter(TRT == "DegradedAround")

# relate to all abiotic variables
deg_ar.lm = lm(BF_rich ~ pH + OM + oC + NO3 + NH4 + tP + H2O + BactRich, data = deg_ar)
summary(deg_ar.lm)

### abiotic/biotic predictors
nat_ap = BF_map %>%
  filter(TRT == "NativeApron")

# relate to all abiotic variables
nat_ap.lm = lm(BF_rich ~ pH + OM + oC + NO3 + NH4 + tP + H2O + BactRich, data = nat_ap)
summary(nat_ap.lm)

nat_ar = BF_map %>%
  filter(TRT == "NativeAround")

# relate to all abiotic variables
nat_ar.lm = lm(BF_rich ~ pH + OM + oC + NO3 + NH4 + tP + H2O + BactRich, data = nat_ar)
summary(nat_ar.lm)

## FF 
deg_ap = FF_map %>%
  filter(TRT == "DegradedApron")

# relate to all abiotic variables
deg_ap.lm = lm(FF_rich ~ pH + OM + oC + NO3 + NH4 + tP + H2O + FungiRich, data = deg_ap)
summary(deg_ap.lm)

deg_ar = FF_map %>%
  filter(TRT == "DegradedAround")

# relate to all abiotic variables
deg_ar.lm = lm(FF_rich ~ pH + OM + oC + NO3 + NH4 + tP + H2O + FungiRich, data = deg_ar)
summary(deg_ar.lm)

nat_ap = FF_map %>%
  filter(TRT == "NativeApron")

# relate to all abiotic variables
nat_ap.lm = lm(FF_rich ~ pH + OM + oC + NO3 + NH4 + tP + H2O + FungiRich, data = nat_ap)
summary(nat_ap.lm)

nat_ar = FF_map %>%
  filter(TRT == "NativeAround")

# relate to all abiotic variables
nat_ar.lm = lm(FF_rich ~ pH + OM + oC + NO3 + NH4 + tP + H2O + FungiRich, data = nat_ar)
summary(nat_ar.lm)
