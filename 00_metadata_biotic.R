### testing soil variables 

library(mctoolsr)
library(tidyverse)
library(emmeans)
library(lme4)
library(lmerTest)
library(effectsize)
library(vegan)

# input data
tax_table_nem = 'nem_asv_gt.txt'
map_fp = 'metadata_fall_gt.txt'

input_nem = load_taxa_table(tax_table_nem, map_fp)

input_nem$map_loaded$sampleID <- row.names(input_nem$map_loaded)
input_nem = filter_data(input_nem, filter_cat = 'sampleID',
                        filter_vals = c("F.10.O_F"))

map = input_nem$map_loaded

# Fungi Richness
fung.rich.lm = lmer(FungiRich ~ Habitat * Ecosystem +(1|Site), data = map)
anova(fung.rich.lm)
emmeans(fung.rich.lm, pairwise~Habitat*Ecosystem)

# plot
map$TRT = factor(map$TRT, 
                 levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(map, aes(x = TRT, y = FungiRich)) +
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
       y = "FungiRich",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# Fungi Diversity
fung.div.lm = lmer(FungDiv ~ Habitat * Ecosystem +(1|Site), data = map)
anova(fung.div.lm)
emmeans(fung.div.lm, pairwise~Habitat*Ecosystem)

# plot
map$TRT = factor(map$TRT, 
                 levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(map, aes(x = TRT, y = FungDiv)) +
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
       y = "Fungal Diversity",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# Bacti Diversity
Bact.div.lm = lmer(BactDiv ~ Habitat * Ecosystem +(1|Site), data = map)
anova(Bact.div.lm)
emmeans(Bact.div.lm, pairwise~Habitat*Ecosystem)

# plot
map$TRT = factor(map$TRT, 
                 levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(map, aes(x = TRT, y = BactDiv)) +
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
       y = "Bacterial Diversity",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

# Bacti Rich
Bact.rich.lm = lmer(BactRich ~ Habitat * Ecosystem +(1|Site), data = map)
anova(Bact.rich.lm)
emmeans(Bact.rich.lm, pairwise~Habitat*Ecosystem)

# plot
map$TRT = factor(map$TRT, 
                 levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))


ggplot(map, aes(x = TRT, y = BactRich)) +
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
       y = "Bacterial Richness",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))


###PLANTS
plant <- read.csv("plant_data_for.csv")
plant_nb_only <- plant %>%
  filter(B.NB == "NB")
plant_nb_only_no_t <- plant_nb_only %>%
  filter(Plot != "T")
plant_no_soil <- plant_nb_only_no_t %>%
  filter(FuncGrp != "S")
plant_no_litter <- plant_no_soil %>%
  filter(FuncGrp != "D")
plant_no_unknown <- plant_no_litter %>%
  filter(SppCode != "UNKNOWN")

plant_no_unknown$PercCov <- as.numeric(plant_no_unknown$PercCov)

plant_wide <- plant_no_unknown %>% 
  select(ID, Habitat, Site, Plot, SppCode, PercCov) %>%                      
  group_by(ID, Habitat, Site, Plot, SppCode) %>% 
  pivot_wider(names_from = SppCode, values_from = PercCov) 

plant_wide <- apply(plant_wide,2,as.character)
write.csv(plant_wide, "plant.wide.csv")

plant_wide_n<-read.csv("plant.wide_n.csv", row.names = 1)

str(plant_wide_n)
mapping <- plant_wide_n[,1:4]
div_calcs <- as.data.frame(t(plant_wide_n))
colnames(div_calcs) <- div_calcs[1,]
div_calcs_f <- div_calcs[5:88,]

write.csv(div_calcs_f, "div_calcs_f.csv")
div_calcs_n <- read.csv("div_calcs_f_n.csv", row.names = 1)

plant_div = calc_diversity(div_calcs_n, metric = "shannon")                   
plant_rich = calc_diversity(div_calcs_n, metric = "richness")

mapping$PlantDiv = plant_div
mapping$PlantRich = plant_rich

# Plant Diversity
Plant.div.lm = lmer(PlantDiv ~ Habitat + (1|Site), data = mapping)
anova(Plant.div.lm)

# Plant Rich
Plant.rich.lm = lmer(PlantRich ~ Habitat +(1|Site), data = mapping)
anova(Plant.rich.lm)

# plot
mapping$Habitat = factor(mapping$Habitat, 
                 levels = c("S", "F"))

ggplot(mapping, aes(x = Habitat, y = PlantDiv)) +
  geom_boxplot(aes(fill = Habitat), outlier.shape = NA,
               color = c("black")) +    
  geom_jitter(aes(fill = Habitat),
              alpha = 0.8,
              shape = 21,
              size = 3,
              set.seed(666),
              color = "black") +
  scale_fill_manual(values = c("green4","green4")) +
  labs(x = NULL,
       y = "Plant Diversity",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(mapping, aes(x = Habitat, y = PlantRich)) +
  geom_boxplot(aes(fill = Habitat), outlier.shape = NA,
               color = c("black")) +    
  geom_jitter(aes(fill = Habitat),
              alpha = 0.8,
              shape = 21,
              size = 3,
              set.seed(666),
              color = "black") +
  scale_fill_manual(values = c("green4","green4")) +
  labs(x = NULL,
       y = "Plant Richness",
       title = NULL) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))


ord_asv <- metaMDS(t(div_calcs_n), distance = "bray", k = 3,
                     trymax = 200)
nmds1 <- ord_asv$points[,1]
nmds2 <- ord_asv$points[,2]

ord_asv_plot <- as.data.frame(cbind(nmds1,nmds2))

ggplot(ord_asv_plot, aes(nmds1, nmds2)) +
  geom_jitter(aes(shape = mapping$Habitat),
              color = "green4", size = 4, fill="lightgrey", stroke=1) +
  scale_shape_manual(values = c(21, 24)) +
  theme_classic() + 
  labs(x="NMDS 1", y="NMDS 2") 

perm <- adonis2(ord_asv_plot ~ mapping$Habitat, strata=mapping$Site)
perm
