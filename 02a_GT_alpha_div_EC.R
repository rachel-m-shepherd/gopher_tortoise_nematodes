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
                           filter_vals = c("F.10.O_F", "F.7.A_F", "S.3.A_F","S.6.O_F"))

# first calc overall nem richness in each sample
nem_rich = calc_diversity(input_nem_ra$data_loaded, metric = 'richness')
nem_div = calc_diversity(input_nem_ra$data_loaded, metric = 'shannon')

sort(nem_rich)
sort(nem_div)


# join to map
nem_map = input_nem_ra$map_loaded
nem_map$nem_rich = nem_rich
nem_map$nem_div = nem_div

nem.rich.lm = lmer(nem_rich ~ Habitat * Ecosystem + (1|Site), data = nem_map)
anova(nem.rich.lm)
emmeans(nem.rich.lm, pairwise~Habitat*Ecosystem)


nem.div.lm = lmer(nem_div ~ Habitat * Ecosystem + (1|Site), data = nem_map)
anova(nem.div.lm)
emmeans(nem.div.lm, pairwise~Habitat*Ecosystem)

# effect size
nem_map$Ecosystem = factor(nem_map$Ecosystem, 
                           levels = c("Native", "Degraded"))

cohens_d(nem_rich ~ Ecosystem, data = nem_map)
cohens_d(nem_div ~ Ecosystem, data = nem_map)

# plot
nem_map$TRT = factor(nem_map$TRT, 
                     levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(nem_map, aes(x = TRT, y = nem_rich)) +
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
       title = "Entire Nematode Community") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(nem_map, aes(x = TRT, y = nem_div)) +
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
       title = "Entire Nematode Community") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

### abiotic/biotic predictors
deg_ap = nem_map %>%
  filter(TRT == "DegradedApron")

# relate to all abiotic variables
deg_ap.lm = lm(nem_rich ~ pH + OM + oC + NO3 + NH4 + tP + H2O + FungiRich + BactRich, data = deg_ap)
summary(deg_ap.lm)

deg_ar = nem_map %>%
  filter(TRT == "DegradedAround")

# relate to all abiotic variables
deg_ar.lm = lm(nem_rich ~ pH + OM + oC + NO3 + NH4 + tP + H2O + FungiRich + BactRich, data = deg_ar)
summary(deg_ar.lm)

### abiotic/biotic predictors
nat_ap = nem_map %>%
  filter(TRT == "NativeApron")

# relate to all abiotic variables
nat_ap.lm = lm(nem_rich ~ pH + OM + oC + NO3 + NH4 + tP + H2O + FungiRich + BactRich, data = nat_ap)
summary(nat_ap.lm)

nat_ar = nem_map %>%
  filter(TRT == "NativeAround")

# relate to all abiotic variables
nat_ar.lm = lm(nem_rich ~ pH + OM + oC + NO3 + NH4 + tP + H2O + FungiRich + BactRich, data = nat_ar)
summary(nat_ar.lm)
