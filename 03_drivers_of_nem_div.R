# drivers of Nematode diversity

library(mctoolsr)
library(tidyverse)
library(emmeans)
library(lme4)
library(lmerTest)
library(effectsize)
library(ggpattern)

# input data
tax_table_nem = 'nem_asv_gt.txt'
map_fp = 'metadata_fall_gt.txt'

input_nem = load_taxa_table(tax_table_nem, map_fp)

input_nem_ra = convert_to_relative_abundances(input_nem)
input_nem_ra$map_loaded$sampleID <- row.names(input_nem_ra$map_loaded)
input_nem_ra = filter_data(input_nem_ra, filter_cat = 'sampleID',
                           filter_vals = c("F.10.O_F"))

EC_rich = calc_diversity(input_nem_ra$data_loaded, metric = 'richness')
EC_div = calc_diversity(input_nem_ra$data_loaded, metric = 'shannon')

EC_map = input_nem_ra$map_loaded
EC_map$EC_rich = EC_rich
EC_map$EC_div = EC_div

# make sub tables for each TRT 

EC_degradedapron = EC_map %>%
  filter(TRT == "DegradedApron")
EC_degradedaround  = EC_map %>%
  filter(TRT == "DegradedAround")
EC_nativeapron = EC_map %>%
  filter(TRT == "NativeApron")
EC_nativearound = EC_map %>%
  filter(TRT == "NativeAround")

mean(EC_nativeapron$EC_div)
sd(EC_nativeapron$EC_div)
mean(EC_nativearound$EC_div)
sd(EC_nativearound$EC_div)
mean(EC_degradedapron$EC_div)
sd(EC_degradedapron$EC_div)
mean(EC_degradedaround$EC_div)
sd(EC_degradedaround$EC_div)

mean(EC_nativeapron$EC_rich)
sd(EC_nativeapron$EC_rich)
mean(EC_nativearound$EC_rich)
sd(EC_nativearound$EC_rich)
mean(EC_degradedapron$EC_rich)
sd(EC_degradedapron$EC_rich)
mean(EC_degradedaround$EC_rich)
sd(EC_degradedaround$EC_rich)


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

# make sub tables for each TRT 

BF_degradedapron = BF_map %>%
  filter(TRT == "DegradedApron")
BF_degradedaround  = BF_map %>%
  filter(TRT == "DegradedAround")
BF_nativeapron = BF_map %>%
  filter(TRT == "NativeApron")
BF_nativearound = BF_map %>%
  filter(TRT == "NativeAround")

mean(BF_nativeapron$BF_div)
sd(BF_nativeapron$BF_div)
mean(BF_nativearound$BF_div)
sd(BF_nativearound$BF_div)
mean(BF_degradedapron$BF_div)
sd(BF_degradedapron$BF_div)
mean(BF_degradedaround$BF_div)
sd(BF_degradedaround$BF_div)

mean(BF_nativeapron$BF_rich)
sd(BF_nativeapron$BF_rich)
mean(BF_nativearound$BF_rich)
sd(BF_nativearound$BF_rich)
mean(BF_degradedapron$BF_rich)
sd(BF_degradedapron$BF_rich)
mean(BF_degradedaround$BF_rich)
sd(BF_degradedaround$BF_rich)

# Richness
# Run correlations 
# Degraded Apron
rich.BF.da_pH = cor.test(x = BF_degradedapron$BF_rich, y = BF_degradedapron$pH, method = "spearman", exact=FALSE)
rich.BF.da_OM = cor.test(x = BF_degradedapron$BF_rich, y = BF_degradedapron$OM, method = "spearman", exact=FALSE)
rich.BF.da_NO3 = cor.test(x = BF_degradedapron$BF_rich, y = BF_degradedapron$NO3, method = "spearman", exact=FALSE)
rich.BF.da_oC = cor.test(x = BF_degradedapron$BF_rich, y = BF_degradedapron$oC, method = "spearman", exact=FALSE)
rich.BF.da_tP = cor.test(x = BF_degradedapron$BF_rich, y = BF_degradedapron$tP, method = "spearman", exact=FALSE)
rich.BF.da_NH4 = cor.test(x = BF_degradedapron$BF_rich, y = BF_degradedapron$NH4, method = "spearman", exact=FALSE)
rich.BF.da_H2O = cor.test(x = BF_degradedapron$BF_rich, y = BF_degradedapron$H2O, method = "spearman", exact=FALSE)
rich.BF.da_BactDiv = cor.test(x = BF_degradedapron$BF_rich, y = BF_degradedapron$BactDiv, method = "spearman", exact=FALSE)
rich.BF.da_BactRich = cor.test(x = BF_degradedapron$BF_rich, y = BF_degradedapron$BactRich, method = "spearman", exact=FALSE)
rich.BF.da_quibit = cor.test(x = BF_degradedapron$BF_rich, y = BF_degradedapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Bacterial Diversity", "Bacterial Richness", 
             "Quibit")
signif = c(rich.BF.da_pH$p.value, rich.BF.da_OM$p.value, rich.BF.da_NO3$p.value,
           rich.BF.da_oC$p.value, rich.BF.da_tP$p.value, rich.BF.da_NH4$p.value,
           rich.BF.da_H2O$p.value, rich.BF.da_BactDiv$p.value, rich.BF.da_BactRich$p.value,
           rich.BF.da_quibit$p.value)

estimate = c(rich.BF.da_pH$estimate, rich.BF.da_OM$estimate, rich.BF.da_NO3$estimate,
             rich.BF.da_oC$estimate, rich.BF.da_tP$estimate, rich.BF.da_NH4$estimate,
             rich.BF.da_H2O$estimate, rich.BF.da_BactDiv$estimate, rich.BF.da_BactRich$estimate,
             rich.BF.da_quibit$estimate)

corr_BF.da_rich <- data.frame(variable, signif, estimate)
corr_BF.da_rich$group <- "BF Richness (Degraded Apron)"

# Run correlations 
# Degraded Around
rich.BF.do_pH = cor.test(x = BF_degradedaround$BF_rich, y = BF_degradedaround$pH, method = "spearman", exact=FALSE)
rich.BF.do_OM = cor.test(x = BF_degradedaround$BF_rich, y = BF_degradedaround$OM, method = "spearman", exact=FALSE)
rich.BF.do_NO3 = cor.test(x = BF_degradedaround$BF_rich, y = BF_degradedaround$NO3, method = "spearman", exact=FALSE)
rich.BF.do_oC = cor.test(x = BF_degradedaround$BF_rich, y = BF_degradedaround$oC, method = "spearman", exact=FALSE)
rich.BF.do_tP = cor.test(x = BF_degradedaround$BF_rich, y = BF_degradedaround$tP, method = "spearman", exact=FALSE)
rich.BF.do_NH4 = cor.test(x = BF_degradedaround$BF_rich, y = BF_degradedaround$NH4, method = "spearman", exact=FALSE)
rich.BF.do_H2O = cor.test(x = BF_degradedaround$BF_rich, y = BF_degradedaround$H2O, method = "spearman", exact=FALSE)
rich.BF.do_BactDiv = cor.test(x = BF_degradedaround$BF_rich, y = BF_degradedaround$BactDiv, method = "spearman", exact=FALSE)
rich.BF.do_BactRich = cor.test(x = BF_degradedaround$BF_rich, y = BF_degradedaround$BactRich, method = "spearman", exact=FALSE)
rich.BF.do_quibit = cor.test(x = BF_degradedaround$BF_rich, y = BF_degradedaround$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Bacterial Diversity", "Bacterial Richness", 
             "Quibit")
signif = c(rich.BF.do_pH$p.value, rich.BF.do_OM$p.value, rich.BF.do_NO3$p.value,
           rich.BF.do_oC$p.value, rich.BF.do_tP$p.value, rich.BF.do_NH4$p.value,
           rich.BF.do_H2O$p.value, rich.BF.do_BactDiv$p.value, rich.BF.do_BactRich$p.value,
           rich.BF.do_quibit$p.value)

estimate = c(rich.BF.do_pH$estimate, rich.BF.do_OM$estimate, rich.BF.do_NO3$estimate,
             rich.BF.do_oC$estimate, rich.BF.do_tP$estimate, rich.BF.do_NH4$estimate,
             rich.BF.do_H2O$estimate, rich.BF.do_BactDiv$estimate, rich.BF.do_BactRich$estimate,
             rich.BF.do_quibit$estimate)

corr_BF.do_rich <- data.frame(variable, signif, estimate)
corr_BF.do_rich$group <- "BF Richness (Degraded Around)"

# Run correlations 
# Native Apron
rich.BF.na_pH = cor.test(x = BF_nativeapron$BF_rich, y = BF_nativeapron$pH, method = "spearman", exact=FALSE)
rich.BF.na_OM = cor.test(x = BF_nativeapron$BF_rich, y = BF_nativeapron$OM, method = "spearman", exact=FALSE)
rich.BF.na_NO3 = cor.test(x = BF_nativeapron$BF_rich, y = BF_nativeapron$NO3, method = "spearman", exact=FALSE)
rich.BF.na_oC = cor.test(x = BF_nativeapron$BF_rich, y = BF_nativeapron$oC, method = "spearman", exact=FALSE)
rich.BF.na_tP = cor.test(x = BF_nativeapron$BF_rich, y = BF_nativeapron$tP, method = "spearman", exact=FALSE)
rich.BF.na_NH4 = cor.test(x = BF_nativeapron$BF_rich, y = BF_nativeapron$NH4, method = "spearman", exact=FALSE)
rich.BF.na_H2O = cor.test(x = BF_nativeapron$BF_rich, y = BF_nativeapron$H2O, method = "spearman", exact=FALSE)
rich.BF.na_BactDiv = cor.test(x = BF_nativeapron$BF_rich, y = BF_nativeapron$BactDiv, method = "spearman", exact=FALSE)
rich.BF.na_BactRich = cor.test(x = BF_nativeapron$BF_rich, y = BF_nativeapron$BactRich, method = "spearman", exact=FALSE)
rich.BF.na_quibit = cor.test(x = BF_nativeapron$BF_rich, y = BF_nativeapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Bacterial Diversity", "Bacterial Richness", 
             "Quibit")
signif = c(rich.BF.na_pH$p.value, rich.BF.na_OM$p.value, rich.BF.na_NO3$p.value,
           rich.BF.na_oC$p.value, rich.BF.na_tP$p.value, rich.BF.na_NH4$p.value,
           rich.BF.na_H2O$p.value, rich.BF.na_BactDiv$p.value, rich.BF.na_BactRich$p.value,
           rich.BF.na_quibit$p.value)

estimate = c(rich.BF.na_pH$estimate, rich.BF.na_OM$estimate, rich.BF.na_NO3$estimate,
             rich.BF.na_oC$estimate, rich.BF.na_tP$estimate, rich.BF.na_NH4$estimate,
             rich.BF.na_H2O$estimate, rich.BF.na_BactDiv$estimate, rich.BF.na_BactRich$estimate,
             rich.BF.na_quibit$estimate)

corr_BF.na_rich <- data.frame(variable, signif, estimate)
corr_BF.na_rich$group <- "BF Richness (Native Apron)"

# Run correlations 
# Native Around
rich.BF.no_pH = cor.test(x = BF_nativearound$BF_rich, y = BF_nativearound$pH, method = "spearman", exact=FALSE)
rich.BF.no_OM = cor.test(x = BF_nativearound$BF_rich, y = BF_nativearound$OM, method = "spearman", exact=FALSE)
rich.BF.no_NO3 = cor.test(x = BF_nativearound$BF_rich, y = BF_nativearound$NO3, method = "spearman", exact=FALSE)
rich.BF.no_oC = cor.test(x = BF_nativearound$BF_rich, y = BF_nativearound$oC, method = "spearman", exact=FALSE)
rich.BF.no_tP = cor.test(x = BF_nativearound$BF_rich, y = BF_nativearound$tP, method = "spearman", exact=FALSE)
rich.BF.no_NH4 = cor.test(x = BF_nativearound$BF_rich, y = BF_nativearound$NH4, method = "spearman", exact=FALSE)
rich.BF.no_H2O = cor.test(x = BF_nativearound$BF_rich, y = BF_nativearound$H2O, method = "spearman", exact=FALSE)
rich.BF.no_BactDiv = cor.test(x = BF_nativearound$BF_rich, y = BF_nativearound$BactDiv, method = "spearman", exact=FALSE)
rich.BF.no_BactRich = cor.test(x = BF_nativearound$BF_rich, y = BF_nativearound$BactRich, method = "spearman", exact=FALSE)
rich.BF.no_quibit = cor.test(x = BF_nativearound$BF_rich, y = BF_nativearound$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Bacterial Diversity", "Bacterial Richness", 
             "Quibit")
signif = c(rich.BF.no_pH$p.value, rich.BF.no_OM$p.value, rich.BF.no_NO3$p.value,
           rich.BF.no_oC$p.value, rich.BF.no_tP$p.value, rich.BF.no_NH4$p.value,
           rich.BF.no_H2O$p.value, rich.BF.no_BactDiv$p.value, rich.BF.no_BactRich$p.value,
           rich.BF.no_quibit$p.value)

estimate = c(rich.BF.no_pH$estimate, rich.BF.no_OM$estimate, rich.BF.no_NO3$estimate,
             rich.BF.no_oC$estimate, rich.BF.no_tP$estimate, rich.BF.no_NH4$estimate,
             rich.BF.no_H2O$estimate, rich.BF.no_BactDiv$estimate, rich.BF.no_BactRich$estimate,
             rich.BF.no_quibit$estimate)

corr_BF.no_rich <- data.frame(variable, signif, estimate)
corr_BF.no_rich$group <- "BF Richness (Native Around)"

# Diversity
# Run correlations 
# Degraded Apron
div.BF.da_pH = cor.test(x = BF_degradedapron$BF_div, y = BF_degradedapron$pH, method = "spearman", exact=FALSE)
div.BF.da_OM = cor.test(x = BF_degradedapron$BF_div, y = BF_degradedapron$OM, method = "spearman", exact=FALSE)
div.BF.da_NO3 = cor.test(x = BF_degradedapron$BF_div, y = BF_degradedapron$NO3, method = "spearman", exact=FALSE)
div.BF.da_oC = cor.test(x = BF_degradedapron$BF_div, y = BF_degradedapron$oC, method = "spearman", exact=FALSE)
div.BF.da_tP = cor.test(x = BF_degradedapron$BF_div, y = BF_degradedapron$tP, method = "spearman", exact=FALSE)
div.BF.da_NH4 = cor.test(x = BF_degradedapron$BF_div, y = BF_degradedapron$NH4, method = "spearman", exact=FALSE)
div.BF.da_H2O = cor.test(x = BF_degradedapron$BF_div, y = BF_degradedapron$H2O, method = "spearman", exact=FALSE)
div.BF.da_BactDiv = cor.test(x = BF_degradedapron$BF_div, y = BF_degradedapron$BactDiv, method = "spearman", exact=FALSE)
div.BF.da_BactRich = cor.test(x = BF_degradedapron$BF_div, y = BF_degradedapron$BactRich, method = "spearman", exact=FALSE)
div.BF.da_quibit = cor.test(x = BF_degradedapron$BF_div, y = BF_degradedapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Bacterial Diversity", "Bacterial Richness", 
             "Quibit")
signif = c(div.BF.da_pH$p.value, div.BF.da_OM$p.value, div.BF.da_NO3$p.value,
           div.BF.da_oC$p.value, div.BF.da_tP$p.value, div.BF.da_NH4$p.value,
           div.BF.da_H2O$p.value, div.BF.da_BactDiv$p.value, div.BF.da_BactRich$p.value,
           div.BF.da_quibit$p.value)

estimate = c(div.BF.da_pH$estimate, div.BF.da_OM$estimate, div.BF.da_NO3$estimate,
             div.BF.da_oC$estimate, div.BF.da_tP$estimate, div.BF.da_NH4$estimate,
             div.BF.da_H2O$estimate, div.BF.da_BactDiv$estimate, div.BF.da_BactRich$estimate,
             div.BF.da_quibit$estimate)

corr_BF.da_div <- data.frame(variable, signif, estimate)
corr_BF.da_div$group <- "BF Diversity (Degraded Apron)"

# Run correlations 
# Degraded Around
div.BF.do_pH = cor.test(x = BF_degradedaround$BF_div, y = BF_degradedaround$pH, method = "spearman", exact=FALSE)
div.BF.do_OM = cor.test(x = BF_degradedaround$BF_div, y = BF_degradedaround$OM, method = "spearman", exact=FALSE)
div.BF.do_NO3 = cor.test(x = BF_degradedaround$BF_div, y = BF_degradedaround$NO3, method = "spearman", exact=FALSE)
div.BF.do_oC = cor.test(x = BF_degradedaround$BF_div, y = BF_degradedaround$oC, method = "spearman", exact=FALSE)
div.BF.do_tP = cor.test(x = BF_degradedaround$BF_div, y = BF_degradedaround$tP, method = "spearman", exact=FALSE)
div.BF.do_NH4 = cor.test(x = BF_degradedaround$BF_div, y = BF_degradedaround$NH4, method = "spearman", exact=FALSE)
div.BF.do_H2O = cor.test(x = BF_degradedaround$BF_div, y = BF_degradedaround$H2O, method = "spearman", exact=FALSE)
div.BF.do_BactDiv = cor.test(x = BF_degradedaround$BF_div, y = BF_degradedaround$BactDiv, method = "spearman", exact=FALSE)
div.BF.do_BactRich = cor.test(x = BF_degradedaround$BF_div, y = BF_degradedaround$BactRich, method = "spearman", exact=FALSE)
div.BF.do_quibit = cor.test(x = BF_degradedaround$BF_div, y = BF_degradedaround$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Bacterial Diversity", "Bacterial Richness", 
             "Quibit")
signif = c(div.BF.do_pH$p.value, div.BF.do_OM$p.value, div.BF.do_NO3$p.value,
           div.BF.do_oC$p.value, div.BF.do_tP$p.value, div.BF.do_NH4$p.value,
           div.BF.do_H2O$p.value, div.BF.do_BactDiv$p.value, div.BF.do_BactRich$p.value,
           div.BF.do_quibit$p.value)

estimate = c(div.BF.do_pH$estimate, div.BF.do_OM$estimate, div.BF.do_NO3$estimate,
             div.BF.do_oC$estimate, div.BF.do_tP$estimate, div.BF.do_NH4$estimate,
             div.BF.do_H2O$estimate, div.BF.do_BactDiv$estimate, div.BF.do_BactRich$estimate,
             div.BF.do_quibit$estimate)

corr_BF.do_div <- data.frame(variable, signif, estimate)
corr_BF.do_div$group <- "BF Diversity (Degraded Around)"

# Run correlations 
# Native Apron
div.BF.na_pH = cor.test(x = BF_nativeapron$BF_div, y = BF_nativeapron$pH, method = "spearman", exact=FALSE)
div.BF.na_OM = cor.test(x = BF_nativeapron$BF_div, y = BF_nativeapron$OM, method = "spearman", exact=FALSE)
div.BF.na_NO3 = cor.test(x = BF_nativeapron$BF_div, y = BF_nativeapron$NO3, method = "spearman", exact=FALSE)
div.BF.na_oC = cor.test(x = BF_nativeapron$BF_div, y = BF_nativeapron$oC, method = "spearman", exact=FALSE)
div.BF.na_tP = cor.test(x = BF_nativeapron$BF_div, y = BF_nativeapron$tP, method = "spearman", exact=FALSE)
div.BF.na_NH4 = cor.test(x = BF_nativeapron$BF_div, y = BF_nativeapron$NH4, method = "spearman", exact=FALSE)
div.BF.na_H2O = cor.test(x = BF_nativeapron$BF_div, y = BF_nativeapron$H2O, method = "spearman", exact=FALSE)
div.BF.na_BactDiv = cor.test(x = BF_nativeapron$BF_div, y = BF_nativeapron$BactDiv, method = "spearman", exact=FALSE)
div.BF.na_BactRich = cor.test(x = BF_nativeapron$BF_div, y = BF_nativeapron$BactRich, method = "spearman", exact=FALSE)
div.BF.na_quibit = cor.test(x = BF_nativeapron$BF_div, y = BF_nativeapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Bacterial Diversity", "Bacterial Richness", 
             "Quibit")
signif = c(div.BF.na_pH$p.value, div.BF.na_OM$p.value, div.BF.na_NO3$p.value,
           div.BF.na_oC$p.value, div.BF.na_tP$p.value, div.BF.na_NH4$p.value,
           div.BF.na_H2O$p.value, div.BF.na_BactDiv$p.value, div.BF.na_BactRich$p.value,
           div.BF.na_quibit$p.value)

estimate = c(div.BF.na_pH$estimate, div.BF.na_OM$estimate, div.BF.na_NO3$estimate,
             div.BF.na_oC$estimate, div.BF.na_tP$estimate, div.BF.na_NH4$estimate,
             div.BF.na_H2O$estimate, div.BF.na_BactDiv$estimate, div.BF.na_BactRich$estimate,
             div.BF.na_quibit$estimate)

corr_BF.na_div <- data.frame(variable, signif, estimate)
corr_BF.na_div$group <- "BF Diversity (Native Apron)"

# Run correlations 
# Native Around
div.BF.no_pH = cor.test(x = BF_nativearound$BF_div, y = BF_nativearound$pH, method = "spearman", exact=FALSE)
div.BF.no_OM = cor.test(x = BF_nativearound$BF_div, y = BF_nativearound$OM, method = "spearman", exact=FALSE)
div.BF.no_NO3 = cor.test(x = BF_nativearound$BF_div, y = BF_nativearound$NO3, method = "spearman", exact=FALSE)
div.BF.no_oC = cor.test(x = BF_nativearound$BF_div, y = BF_nativearound$oC, method = "spearman", exact=FALSE)
div.BF.no_tP = cor.test(x = BF_nativearound$BF_div, y = BF_nativearound$tP, method = "spearman", exact=FALSE)
div.BF.no_NH4 = cor.test(x = BF_nativearound$BF_div, y = BF_nativearound$NH4, method = "spearman", exact=FALSE)
div.BF.no_H2O = cor.test(x = BF_nativearound$BF_div, y = BF_nativearound$H2O, method = "spearman", exact=FALSE)
div.BF.no_BactDiv = cor.test(x = BF_nativearound$BF_div, y = BF_nativearound$BactDiv, method = "spearman", exact=FALSE)
div.BF.no_BactRich = cor.test(x = BF_nativearound$BF_div, y = BF_nativearound$BactRich, method = "spearman", exact=FALSE)
div.BF.no_quibit = cor.test(x = BF_nativearound$BF_div, y = BF_nativearound$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Bacterial Diversity", "Bacterial Richness", 
             "Quibit")
signif = c(div.BF.no_pH$p.value, div.BF.no_OM$p.value, div.BF.no_NO3$p.value,
           div.BF.no_oC$p.value, div.BF.no_tP$p.value, div.BF.no_NH4$p.value,
           div.BF.no_H2O$p.value, div.BF.no_BactDiv$p.value, div.BF.no_BactRich$p.value,
           div.BF.no_quibit$p.value)

estimate = c(div.BF.no_pH$estimate, div.BF.no_OM$estimate, div.BF.no_NO3$estimate,
             div.BF.no_oC$estimate, div.BF.no_tP$estimate, div.BF.no_NH4$estimate,
             div.BF.no_H2O$estimate, div.BF.no_BactDiv$estimate, div.BF.no_BactRich$estimate,
             div.BF.no_quibit$estimate)

corr_BF.no_div <- data.frame(variable, signif, estimate)
corr_BF.no_div$group <- "BF Diversity (Native Around)"


# calc nem richness in each sample for major groups
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

# make sub tables for each TRT 

FF_degradedapron = FF_map %>%
  filter(TRT == "DegradedApron")
FF_degradedaround  = FF_map %>%
  filter(TRT == "DegradedAround")
FF_nativeapron = FF_map %>%
  filter(TRT == "NativeApron")
FF_nativearound = FF_map %>%
  filter(TRT == "NativeAround")

mean(FF_nativeapron$FF_div)
sd(FF_nativeapron$FF_div)
mean(FF_nativearound$FF_div)
sd(FF_nativearound$FF_div)
mean(FF_degradedapron$FF_div)
sd(FF_degradedapron$FF_div)
mean(FF_degradedaround$FF_div)
sd(FF_degradedaround$FF_div)

mean(FF_nativeapron$FF_rich)
sd(FF_nativeapron$FF_rich)
mean(FF_nativearound$FF_rich)
sd(FF_nativearound$FF_rich)
mean(FF_degradedapron$FF_rich)
sd(FF_degradedapron$FF_rich)
mean(FF_degradedaround$FF_rich)
sd(FF_degradedaround$FF_rich)

# Richness
# Run correlations 
# Degraded Apron
rich.FF.da_pH = cor.test(x = FF_degradedapron$FF_rich, y = FF_degradedapron$pH, method = "spearman", exact=FALSE)
rich.FF.da_OM = cor.test(x = FF_degradedapron$FF_rich, y = FF_degradedapron$OM, method = "spearman", exact=FALSE)
rich.FF.da_NO3 = cor.test(x = FF_degradedapron$FF_rich, y = FF_degradedapron$NO3, method = "spearman", exact=FALSE)
rich.FF.da_oC = cor.test(x = FF_degradedapron$FF_rich, y = FF_degradedapron$oC, method = "spearman", exact=FALSE)
rich.FF.da_tP = cor.test(x = FF_degradedapron$FF_rich, y = FF_degradedapron$tP, method = "spearman", exact=FALSE)
rich.FF.da_NH4 = cor.test(x = FF_degradedapron$FF_rich, y = FF_degradedapron$NH4, method = "spearman", exact=FALSE)
rich.FF.da_H2O = cor.test(x = FF_degradedapron$FF_rich, y = FF_degradedapron$H2O, method = "spearman", exact=FALSE)
rich.FF.da_FungDiv = cor.test(x = FF_degradedapron$FF_rich, y = FF_degradedapron$FungDiv, method = "spearman", exact=FALSE)
rich.FF.da_FungiRich = cor.test(x = FF_degradedapron$FF_rich, y = FF_degradedapron$FungiRich, method = "spearman", exact=FALSE)
rich.FF.da_quibit = cor.test(x = FF_degradedapron$FF_rich, y = FF_degradedapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit")
signif = c(rich.FF.da_pH$p.value, rich.FF.da_OM$p.value, rich.FF.da_NO3$p.value,
           rich.FF.da_oC$p.value, rich.FF.da_tP$p.value, rich.FF.da_NH4$p.value,
           rich.FF.da_H2O$p.value, rich.FF.da_FungDiv$p.value, rich.FF.da_FungiRich$p.value,
           rich.FF.da_quibit$p.value)

estimate = c(rich.FF.da_pH$estimate, rich.FF.da_OM$estimate, rich.FF.da_NO3$estimate,
             rich.FF.da_oC$estimate, rich.FF.da_tP$estimate, rich.FF.da_NH4$estimate,
             rich.FF.da_H2O$estimate, rich.FF.da_FungDiv$estimate, rich.FF.da_FungiRich$estimate,
             rich.FF.da_quibit$estimate)

corr_FF.da_rich <- data.frame(variable, signif, estimate)
corr_FF.da_rich$group <- "FF Richness (Degraded Apron)"

# Run correlations 
# Degraded Around
rich.FF.do_pH = cor.test(x = FF_degradedaround$FF_rich, y = FF_degradedaround$pH, method = "spearman", exact=FALSE)
rich.FF.do_OM = cor.test(x = FF_degradedaround$FF_rich, y = FF_degradedaround$OM, method = "spearman", exact=FALSE)
rich.FF.do_NO3 = cor.test(x = FF_degradedaround$FF_rich, y = FF_degradedaround$NO3, method = "spearman", exact=FALSE)
rich.FF.do_oC = cor.test(x = FF_degradedaround$FF_rich, y = FF_degradedaround$oC, method = "spearman", exact=FALSE)
rich.FF.do_tP = cor.test(x = FF_degradedaround$FF_rich, y = FF_degradedaround$tP, method = "spearman", exact=FALSE)
rich.FF.do_NH4 = cor.test(x = FF_degradedaround$FF_rich, y = FF_degradedaround$NH4, method = "spearman", exact=FALSE)
rich.FF.do_H2O = cor.test(x = FF_degradedaround$FF_rich, y = FF_degradedaround$H2O, method = "spearman", exact=FALSE)
rich.FF.do_FungDiv = cor.test(x = FF_degradedaround$FF_rich, y = FF_degradedaround$FungDiv, method = "spearman", exact=FALSE)
rich.FF.do_FungiRich = cor.test(x = FF_degradedaround$FF_rich, y = FF_degradedaround$FungiRich, method = "spearman", exact=FALSE)
rich.FF.do_quibit = cor.test(x = FF_degradedaround$FF_rich, y = FF_degradedaround$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit")
signif = c(rich.FF.do_pH$p.value, rich.FF.do_OM$p.value, rich.FF.do_NO3$p.value,
           rich.FF.do_oC$p.value, rich.FF.do_tP$p.value, rich.FF.do_NH4$p.value,
           rich.FF.do_H2O$p.value, rich.FF.do_FungDiv$p.value, rich.FF.do_FungiRich$p.value,
           rich.FF.do_quibit$p.value)

estimate = c(rich.FF.do_pH$estimate, rich.FF.do_OM$estimate, rich.FF.do_NO3$estimate,
             rich.FF.do_oC$estimate, rich.FF.do_tP$estimate, rich.FF.do_NH4$estimate,
             rich.FF.do_H2O$estimate, rich.FF.do_FungDiv$estimate, rich.FF.do_FungiRich$estimate,
             rich.FF.do_quibit$estimate)

corr_FF.do_rich <- data.frame(variable, signif, estimate)
corr_FF.do_rich$group <- "FF Richness (Degraded Around)"

# Run correlations 
# Native Apron
rich.FF.na_pH = cor.test(x = FF_nativeapron$FF_rich, y = FF_nativeapron$pH, method = "spearman", exact=FALSE)
rich.FF.na_OM = cor.test(x = FF_nativeapron$FF_rich, y = FF_nativeapron$OM, method = "spearman", exact=FALSE)
rich.FF.na_NO3 = cor.test(x = FF_nativeapron$FF_rich, y = FF_nativeapron$NO3, method = "spearman", exact=FALSE)
rich.FF.na_oC = cor.test(x = FF_nativeapron$FF_rich, y = FF_nativeapron$oC, method = "spearman", exact=FALSE)
rich.FF.na_tP = cor.test(x = FF_nativeapron$FF_rich, y = FF_nativeapron$tP, method = "spearman", exact=FALSE)
rich.FF.na_NH4 = cor.test(x = FF_nativeapron$FF_rich, y = FF_nativeapron$NH4, method = "spearman", exact=FALSE)
rich.FF.na_H2O = cor.test(x = FF_nativeapron$FF_rich, y = FF_nativeapron$H2O, method = "spearman", exact=FALSE)
rich.FF.na_FungDiv = cor.test(x = FF_nativeapron$FF_rich, y = FF_nativeapron$FungDiv, method = "spearman", exact=FALSE)
rich.FF.na_FungiRich = cor.test(x = FF_nativeapron$FF_rich, y = FF_nativeapron$FungiRich, method = "spearman", exact=FALSE)
rich.FF.na_quibit = cor.test(x = FF_nativeapron$FF_rich, y = FF_nativeapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit")
signif = c(rich.FF.na_pH$p.value, rich.FF.na_OM$p.value, rich.FF.na_NO3$p.value,
           rich.FF.na_oC$p.value, rich.FF.na_tP$p.value, rich.FF.na_NH4$p.value,
           rich.FF.na_H2O$p.value, rich.FF.na_FungDiv$p.value, rich.FF.na_FungiRich$p.value,
           rich.FF.na_quibit$p.value)

estimate = c(rich.FF.na_pH$estimate, rich.FF.na_OM$estimate, rich.FF.na_NO3$estimate,
             rich.FF.na_oC$estimate, rich.FF.na_tP$estimate, rich.FF.na_NH4$estimate,
             rich.FF.na_H2O$estimate, rich.FF.na_FungDiv$estimate, rich.FF.na_FungiRich$estimate,
             rich.FF.na_quibit$estimate)

corr_FF.na_rich <- data.frame(variable, signif, estimate)
corr_FF.na_rich$group <- "FF Richness (Native Apron)"

# Run correlations 
# Native Around
rich.FF.no_pH = cor.test(x = FF_nativearound$FF_rich, y = FF_nativearound$pH, method = "spearman", exact=FALSE)
rich.FF.no_OM = cor.test(x = FF_nativearound$FF_rich, y = FF_nativearound$OM, method = "spearman", exact=FALSE)
rich.FF.no_NO3 = cor.test(x = FF_nativearound$FF_rich, y = FF_nativearound$NO3, method = "spearman", exact=FALSE)
rich.FF.no_oC = cor.test(x = FF_nativearound$FF_rich, y = FF_nativearound$oC, method = "spearman", exact=FALSE)
rich.FF.no_tP = cor.test(x = FF_nativearound$FF_rich, y = FF_nativearound$tP, method = "spearman", exact=FALSE)
rich.FF.no_NH4 = cor.test(x = FF_nativearound$FF_rich, y = FF_nativearound$NH4, method = "spearman", exact=FALSE)
rich.FF.no_H2O = cor.test(x = FF_nativearound$FF_rich, y = FF_nativearound$H2O, method = "spearman", exact=FALSE)
rich.FF.no_FungDiv = cor.test(x = FF_nativearound$FF_rich, y = FF_nativearound$FungDiv, method = "spearman", exact=FALSE)
rich.FF.no_FungiRich = cor.test(x = FF_nativearound$FF_rich, y = FF_nativearound$FungiRich, method = "spearman", exact=FALSE)
rich.FF.no_quibit = cor.test(x = FF_nativearound$FF_rich, y = FF_nativearound$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit")
signif = c(rich.FF.no_pH$p.value, rich.FF.no_OM$p.value, rich.FF.no_NO3$p.value,
           rich.FF.no_oC$p.value, rich.FF.no_tP$p.value, rich.FF.no_NH4$p.value,
           rich.FF.no_H2O$p.value, rich.FF.no_FungDiv$p.value, rich.FF.no_FungiRich$p.value,
           rich.FF.no_quibit$p.value)

estimate = c(rich.FF.no_pH$estimate, rich.FF.no_OM$estimate, rich.FF.no_NO3$estimate,
             rich.FF.no_oC$estimate, rich.FF.no_tP$estimate, rich.FF.no_NH4$estimate,
             rich.FF.no_H2O$estimate, rich.FF.no_FungDiv$estimate, rich.FF.no_FungiRich$estimate,
             rich.FF.no_quibit$estimate)

corr_FF.no_rich <- data.frame(variable, signif, estimate)
corr_FF.no_rich$group <- "FF Richness (Native Around)"

# Diversity
# Run correlations 
# Degraded Apron
div.FF.da_pH = cor.test(x = FF_degradedapron$FF_div, y = FF_degradedapron$pH, method = "spearman", exact=FALSE)
div.FF.da_OM = cor.test(x = FF_degradedapron$FF_div, y = FF_degradedapron$OM, method = "spearman", exact=FALSE)
div.FF.da_NO3 = cor.test(x = FF_degradedapron$FF_div, y = FF_degradedapron$NO3, method = "spearman", exact=FALSE)
div.FF.da_oC = cor.test(x = FF_degradedapron$FF_div, y = FF_degradedapron$oC, method = "spearman", exact=FALSE)
div.FF.da_tP = cor.test(x = FF_degradedapron$FF_div, y = FF_degradedapron$tP, method = "spearman", exact=FALSE)
div.FF.da_NH4 = cor.test(x = FF_degradedapron$FF_div, y = FF_degradedapron$NH4, method = "spearman", exact=FALSE)
div.FF.da_H2O = cor.test(x = FF_degradedapron$FF_div, y = FF_degradedapron$H2O, method = "spearman", exact=FALSE)
div.FF.da_FungDiv = cor.test(x = FF_degradedapron$FF_div, y = FF_degradedapron$FungDiv, method = "spearman", exact=FALSE)
div.FF.da_FungiRich = cor.test(x = FF_degradedapron$FF_div, y = FF_degradedapron$FungiRich, method = "spearman", exact=FALSE)
div.FF.da_quibit = cor.test(x = FF_degradedapron$FF_div, y = FF_degradedapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit")
signif = c(div.FF.da_pH$p.value, div.FF.da_OM$p.value, div.FF.da_NO3$p.value,
           div.FF.da_oC$p.value, div.FF.da_tP$p.value, div.FF.da_NH4$p.value,
           div.FF.da_H2O$p.value, div.FF.da_FungDiv$p.value, div.FF.da_FungiRich$p.value,
           div.FF.da_quibit$p.value)

estimate = c(div.FF.da_pH$estimate, div.FF.da_OM$estimate, div.FF.da_NO3$estimate,
             div.FF.da_oC$estimate, div.FF.da_tP$estimate, div.FF.da_NH4$estimate,
             div.FF.da_H2O$estimate, div.FF.da_FungDiv$estimate, div.FF.da_FungiRich$estimate,
             div.FF.da_quibit$estimate)

corr_FF.da_div <- data.frame(variable, signif, estimate)
corr_FF.da_div$group <- "FF Diversity (Degraded Apron)"

# Run correlations 
# Degraded Around
div.FF.do_pH = cor.test(x = FF_degradedaround$FF_div, y = FF_degradedaround$pH, method = "spearman", exact=FALSE)
div.FF.do_OM = cor.test(x = FF_degradedaround$FF_div, y = FF_degradedaround$OM, method = "spearman", exact=FALSE)
div.FF.do_NO3 = cor.test(x = FF_degradedaround$FF_div, y = FF_degradedaround$NO3, method = "spearman", exact=FALSE)
div.FF.do_oC = cor.test(x = FF_degradedaround$FF_div, y = FF_degradedaround$oC, method = "spearman", exact=FALSE)
div.FF.do_tP = cor.test(x = FF_degradedaround$FF_div, y = FF_degradedaround$tP, method = "spearman", exact=FALSE)
div.FF.do_NH4 = cor.test(x = FF_degradedaround$FF_div, y = FF_degradedaround$NH4, method = "spearman", exact=FALSE)
div.FF.do_H2O = cor.test(x = FF_degradedaround$FF_div, y = FF_degradedaround$H2O, method = "spearman", exact=FALSE)
div.FF.do_FungDiv = cor.test(x = FF_degradedaround$FF_div, y = FF_degradedaround$FungDiv, method = "spearman", exact=FALSE)
div.FF.do_FungiRich = cor.test(x = FF_degradedaround$FF_div, y = FF_degradedaround$FungiRich, method = "spearman", exact=FALSE)
div.FF.do_quibit = cor.test(x = FF_degradedaround$FF_div, y = FF_degradedaround$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit")
signif = c(div.FF.do_pH$p.value, div.FF.do_OM$p.value, div.FF.do_NO3$p.value,
           div.FF.do_oC$p.value, div.FF.do_tP$p.value, div.FF.do_NH4$p.value,
           div.FF.do_H2O$p.value, div.FF.do_FungDiv$p.value, div.FF.do_FungiRich$p.value,
           div.FF.do_quibit$p.value)

estimate = c(div.FF.do_pH$estimate, div.FF.do_OM$estimate, div.FF.do_NO3$estimate,
             div.FF.do_oC$estimate, div.FF.do_tP$estimate, div.FF.do_NH4$estimate,
             div.FF.do_H2O$estimate, div.FF.do_FungDiv$estimate, div.FF.do_FungiRich$estimate,
             div.FF.do_quibit$estimate)

corr_FF.do_div <- data.frame(variable, signif, estimate)
corr_FF.do_div$group <- "FF Diversity (Degraded Around)"

# Run correlations 
# Native Apron
div.FF.na_pH = cor.test(x = FF_nativeapron$FF_div, y = FF_nativeapron$pH, method = "spearman", exact=FALSE)
div.FF.na_OM = cor.test(x = FF_nativeapron$FF_div, y = FF_nativeapron$OM, method = "spearman", exact=FALSE)
div.FF.na_NO3 = cor.test(x = FF_nativeapron$FF_div, y = FF_nativeapron$NO3, method = "spearman", exact=FALSE)
div.FF.na_oC = cor.test(x = FF_nativeapron$FF_div, y = FF_nativeapron$oC, method = "spearman", exact=FALSE)
div.FF.na_tP = cor.test(x = FF_nativeapron$FF_div, y = FF_nativeapron$tP, method = "spearman", exact=FALSE)
div.FF.na_NH4 = cor.test(x = FF_nativeapron$FF_div, y = FF_nativeapron$NH4, method = "spearman", exact=FALSE)
div.FF.na_H2O = cor.test(x = FF_nativeapron$FF_div, y = FF_nativeapron$H2O, method = "spearman", exact=FALSE)
div.FF.na_FungDiv = cor.test(x = FF_nativeapron$FF_div, y = FF_nativeapron$FungDiv, method = "spearman", exact=FALSE)
div.FF.na_FungiRich = cor.test(x = FF_nativeapron$FF_div, y = FF_nativeapron$FungiRich, method = "spearman", exact=FALSE)
div.FF.na_quibit = cor.test(x = FF_nativeapron$FF_div, y = FF_nativeapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit")
signif = c(div.FF.na_pH$p.value, div.FF.na_OM$p.value, div.FF.na_NO3$p.value,
           div.FF.na_oC$p.value, div.FF.na_tP$p.value, div.FF.na_NH4$p.value,
           div.FF.na_H2O$p.value, div.FF.na_FungDiv$p.value, div.FF.na_FungiRich$p.value,
           div.FF.na_quibit$p.value)

estimate = c(div.FF.na_pH$estimate, div.FF.na_OM$estimate, div.FF.na_NO3$estimate,
             div.FF.na_oC$estimate, div.FF.na_tP$estimate, div.FF.na_NH4$estimate,
             div.FF.na_H2O$estimate, div.FF.na_FungDiv$estimate, div.FF.na_FungiRich$estimate,
             div.FF.na_quibit$estimate)

corr_FF.na_div <- data.frame(variable, signif, estimate)
corr_FF.na_div$group <- "FF Diversity (Native Apron)"

# Run correlations 
# Native Around
div.FF.no_pH = cor.test(x = FF_nativearound$FF_div, y = FF_nativearound$pH, method = "spearman", exact=FALSE)
div.FF.no_OM = cor.test(x = FF_nativearound$FF_div, y = FF_nativearound$OM, method = "spearman", exact=FALSE)
div.FF.no_NO3 = cor.test(x = FF_nativearound$FF_div, y = FF_nativearound$NO3, method = "spearman", exact=FALSE)
div.FF.no_oC = cor.test(x = FF_nativearound$FF_div, y = FF_nativearound$oC, method = "spearman", exact=FALSE)
div.FF.no_tP = cor.test(x = FF_nativearound$FF_div, y = FF_nativearound$tP, method = "spearman", exact=FALSE)
div.FF.no_NH4 = cor.test(x = FF_nativearound$FF_div, y = FF_nativearound$NH4, method = "spearman", exact=FALSE)
div.FF.no_H2O = cor.test(x = FF_nativearound$FF_div, y = FF_nativearound$H2O, method = "spearman", exact=FALSE)
div.FF.no_FungDiv = cor.test(x = FF_nativearound$FF_div, y = FF_nativearound$FungDiv, method = "spearman", exact=FALSE)
div.FF.no_FungiRich = cor.test(x = FF_nativearound$FF_div, y = FF_nativearound$FungiRich, method = "spearman", exact=FALSE)
div.FF.no_quibit = cor.test(x = FF_nativearound$FF_div, y = FF_nativearound$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit")
signif = c(div.FF.no_pH$p.value, div.FF.no_OM$p.value, div.FF.no_NO3$p.value,
           div.FF.no_oC$p.value, div.FF.no_tP$p.value, div.FF.no_NH4$p.value,
           div.FF.no_H2O$p.value, div.FF.no_FungDiv$p.value, div.FF.no_FungiRich$p.value,
           div.FF.no_quibit$p.value)

estimate = c(div.FF.no_pH$estimate, div.FF.no_OM$estimate, div.FF.no_NO3$estimate,
             div.FF.no_oC$estimate, div.FF.no_tP$estimate, div.FF.no_NH4$estimate,
             div.FF.no_H2O$estimate, div.FF.no_FungDiv$estimate, div.FF.no_FungiRich$estimate,
             div.FF.no_quibit$estimate)

corr_FF.no_div <- data.frame(variable, signif, estimate)
corr_FF.no_div$group <- "FF Diversity (Native Around)"

# calc nem richness in each sample for major groups
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

# make sub tables for each TRT 

RA_degradedapron = RA_map %>%
  filter(TRT == "DegradedApron")
RA_degradedaround  = RA_map %>%
  filter(TRT == "DegradedAround")
RA_nativeapron = RA_map %>%
  filter(TRT == "NativeApron")
RA_nativearound = RA_map %>%
  filter(TRT == "NativeAround")

mean(RA_nativeapron$RA_div)
sd(RA_nativeapron$RA_div)
mean(RA_nativearound$RA_div)
sd(RA_nativearound$RA_div)
mean(RA_degradedapron$RA_div)
sd(RA_degradedapron$RA_div)
mean(RA_degradedaround$RA_div)
sd(RA_degradedaround$RA_div)

mean(RA_nativeapron$RA_rich)
sd(RA_nativeapron$RA_rich)
mean(RA_nativearound$RA_rich)
sd(RA_nativearound$RA_rich)
mean(RA_degradedapron$RA_rich)
sd(RA_degradedapron$RA_rich)
mean(RA_degradedaround$RA_rich)
sd(RA_degradedaround$RA_rich)

# Richness
# Run correlations 
# Degraded Apron
rich.RA.da_pH = cor.test(x = RA_degradedapron$RA_rich, y = RA_degradedapron$pH, method = "spearman", exact=FALSE)
rich.RA.da_OM = cor.test(x = RA_degradedapron$RA_rich, y = RA_degradedapron$OM, method = "spearman", exact=FALSE)
rich.RA.da_NO3 = cor.test(x = RA_degradedapron$RA_rich, y = RA_degradedapron$NO3, method = "spearman", exact=FALSE)
rich.RA.da_oC = cor.test(x = RA_degradedapron$RA_rich, y = RA_degradedapron$oC, method = "spearman", exact=FALSE)
rich.RA.da_tP = cor.test(x = RA_degradedapron$RA_rich, y = RA_degradedapron$tP, method = "spearman", exact=FALSE)
rich.RA.da_NH4 = cor.test(x = RA_degradedapron$RA_rich, y = RA_degradedapron$NH4, method = "spearman", exact=FALSE)
rich.RA.da_H2O = cor.test(x = RA_degradedapron$RA_rich, y = RA_degradedapron$H2O, method = "spearman", exact=FALSE)
rich.RA.da_FungDiv = cor.test(x = RA_degradedapron$RA_rich, y = RA_degradedapron$FungDiv, method = "spearman", exact=FALSE)
rich.RA.da_FungiRich = cor.test(x = RA_degradedapron$RA_rich, y = RA_degradedapron$FungiRich, method = "spearman", exact=FALSE)
rich.RA.da_quibit = cor.test(x = RA_degradedapron$RA_rich, y = RA_degradedapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit")
signif = c(rich.RA.da_pH$p.value, rich.RA.da_OM$p.value, rich.RA.da_NO3$p.value,
           rich.RA.da_oC$p.value, rich.RA.da_tP$p.value, rich.RA.da_NH4$p.value,
           rich.RA.da_H2O$p.value, rich.RA.da_FungDiv$p.value, rich.RA.da_FungiRich$p.value,
           rich.RA.da_quibit$p.value)

estimate = c(rich.RA.da_pH$estimate, rich.RA.da_OM$estimate, rich.RA.da_NO3$estimate,
             rich.RA.da_oC$estimate, rich.RA.da_tP$estimate, rich.RA.da_NH4$estimate,
             rich.RA.da_H2O$estimate, rich.RA.da_FungDiv$estimate, rich.RA.da_FungiRich$estimate,
             rich.RA.da_quibit$estimate)

corr_RA.da_rich <- data.frame(variable, signif, estimate)
corr_RA.da_rich$group <- "RA Richness (Degraded Apron)"

# Run correlations 
# Degraded Around
rich.RA.do_pH = cor.test(x = RA_degradedaround$RA_rich, y = RA_degradedaround$pH, method = "spearman", exact=FALSE)
rich.RA.do_OM = cor.test(x = RA_degradedaround$RA_rich, y = RA_degradedaround$OM, method = "spearman", exact=FALSE)
rich.RA.do_NO3 = cor.test(x = RA_degradedaround$RA_rich, y = RA_degradedaround$NO3, method = "spearman", exact=FALSE)
rich.RA.do_oC = cor.test(x = RA_degradedaround$RA_rich, y = RA_degradedaround$oC, method = "spearman", exact=FALSE)
rich.RA.do_tP = cor.test(x = RA_degradedaround$RA_rich, y = RA_degradedaround$tP, method = "spearman", exact=FALSE)
rich.RA.do_NH4 = cor.test(x = RA_degradedaround$RA_rich, y = RA_degradedaround$NH4, method = "spearman", exact=FALSE)
rich.RA.do_H2O = cor.test(x = RA_degradedaround$RA_rich, y = RA_degradedaround$H2O, method = "spearman", exact=FALSE)
rich.RA.do_FungDiv = cor.test(x = RA_degradedaround$RA_rich, y = RA_degradedaround$FungDiv, method = "spearman", exact=FALSE)
rich.RA.do_FungiRich = cor.test(x = RA_degradedaround$RA_rich, y = RA_degradedaround$FungiRich, method = "spearman", exact=FALSE)
rich.RA.do_quibit = cor.test(x = RA_degradedaround$RA_rich, y = RA_degradedaround$Qubit, method = "spearman", exact=FALSE)
rich.RA.do_PlantRich = cor.test(x = RA_degradedaround$RA_rich, y = RA_degradedaround$PlantRich, method = "spearman", exact=FALSE)
rich.RA.do_PlantDiv = cor.test(x = RA_degradedaround$RA_rich, y = RA_degradedaround$PlantDiv, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit", "Plant Diversity","Plant Richness")
signif = c(rich.RA.do_pH$p.value, rich.RA.do_OM$p.value, rich.RA.do_NO3$p.value,
           rich.RA.do_oC$p.value, rich.RA.do_tP$p.value, rich.RA.do_NH4$p.value,
           rich.RA.do_H2O$p.value, rich.RA.do_FungDiv$p.value, rich.RA.do_FungiRich$p.value,
           rich.RA.do_quibit$p.value, rich.RA.do_PlantDiv$p.value, rich.RA.do_PlantRich$p.value)

estimate = c(rich.RA.do_pH$estimate, rich.RA.do_OM$estimate, rich.RA.do_NO3$estimate,
             rich.RA.do_oC$estimate, rich.RA.do_tP$estimate, rich.RA.do_NH4$estimate,
             rich.RA.do_H2O$estimate, rich.RA.do_FungDiv$estimate, rich.RA.do_FungiRich$estimate,
             rich.RA.do_quibit$estimate, rich.RA.do_PlantDiv$estimate, rich.RA.do_PlantRich$estimate)

corr_RA.do_rich <- data.frame(variable, signif, estimate)
corr_RA.do_rich$group <- "RA Richness (Degraded Around)"

# Run correlations 
# Native Apron
rich.RA.na_pH = cor.test(x = RA_nativeapron$RA_rich, y = RA_nativeapron$pH, method = "spearman", exact=FALSE)
rich.RA.na_OM = cor.test(x = RA_nativeapron$RA_rich, y = RA_nativeapron$OM, method = "spearman", exact=FALSE)
rich.RA.na_NO3 = cor.test(x = RA_nativeapron$RA_rich, y = RA_nativeapron$NO3, method = "spearman", exact=FALSE)
rich.RA.na_oC = cor.test(x = RA_nativeapron$RA_rich, y = RA_nativeapron$oC, method = "spearman", exact=FALSE)
rich.RA.na_tP = cor.test(x = RA_nativeapron$RA_rich, y = RA_nativeapron$tP, method = "spearman", exact=FALSE)
rich.RA.na_NH4 = cor.test(x = RA_nativeapron$RA_rich, y = RA_nativeapron$NH4, method = "spearman", exact=FALSE)
rich.RA.na_H2O = cor.test(x = RA_nativeapron$RA_rich, y = RA_nativeapron$H2O, method = "spearman", exact=FALSE)
rich.RA.na_FungDiv = cor.test(x = RA_nativeapron$RA_rich, y = RA_nativeapron$FungDiv, method = "spearman", exact=FALSE)
rich.RA.na_FungiRich = cor.test(x = RA_nativeapron$RA_rich, y = RA_nativeapron$FungiRich, method = "spearman", exact=FALSE)
rich.RA.na_quibit = cor.test(x = RA_nativeapron$RA_rich, y = RA_nativeapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit")
signif = c(rich.RA.na_pH$p.value, rich.RA.na_OM$p.value, rich.RA.na_NO3$p.value,
           rich.RA.na_oC$p.value, rich.RA.na_tP$p.value, rich.RA.na_NH4$p.value,
           rich.RA.na_H2O$p.value, rich.RA.na_FungDiv$p.value, rich.RA.na_FungiRich$p.value,
           rich.RA.na_quibit$p.value)

estimate = c(rich.RA.na_pH$estimate, rich.RA.na_OM$estimate, rich.RA.na_NO3$estimate,
             rich.RA.na_oC$estimate, rich.RA.na_tP$estimate, rich.RA.na_NH4$estimate,
             rich.RA.na_H2O$estimate, rich.RA.na_FungDiv$estimate, rich.RA.na_FungiRich$estimate,
             rich.RA.na_quibit$estimate)

corr_RA.na_rich <- data.frame(variable, signif, estimate)
corr_RA.na_rich$group <- "RA Richness (Native Apron)"

# Run correlations 
# Native Around
rich.RA.no_pH = cor.test(x = RA_nativearound$RA_rich, y = RA_nativearound$pH, method = "spearman", exact=FALSE)
rich.RA.no_OM = cor.test(x = RA_nativearound$RA_rich, y = RA_nativearound$OM, method = "spearman", exact=FALSE)
rich.RA.no_NO3 = cor.test(x = RA_nativearound$RA_rich, y = RA_nativearound$NO3, method = "spearman", exact=FALSE)
rich.RA.no_oC = cor.test(x = RA_nativearound$RA_rich, y = RA_nativearound$oC, method = "spearman", exact=FALSE)
rich.RA.no_tP = cor.test(x = RA_nativearound$RA_rich, y = RA_nativearound$tP, method = "spearman", exact=FALSE)
rich.RA.no_NH4 = cor.test(x = RA_nativearound$RA_rich, y = RA_nativearound$NH4, method = "spearman", exact=FALSE)
rich.RA.no_H2O = cor.test(x = RA_nativearound$RA_rich, y = RA_nativearound$H2O, method = "spearman", exact=FALSE)
rich.RA.no_FungDiv = cor.test(x = RA_nativearound$RA_rich, y = RA_nativearound$FungDiv, method = "spearman", exact=FALSE)
rich.RA.no_FungiRich = cor.test(x = RA_nativearound$RA_rich, y = RA_nativearound$FungiRich, method = "spearman", exact=FALSE)
rich.RA.no_quibit = cor.test(x = RA_nativearound$RA_rich, y = RA_nativearound$Qubit, method = "spearman", exact=FALSE)
rich.RA.no_PlantRich = cor.test(x = RA_nativearound$RA_rich, y = RA_nativearound$PlantRich, method = "spearman", exact=FALSE)
rich.RA.no_PlantDiv = cor.test(x = RA_nativearound$RA_rich, y = RA_nativearound$PlantDiv, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit", "Plant Diversity","Plant Richness")
signif = c(rich.RA.no_pH$p.value, rich.RA.no_OM$p.value, rich.RA.no_NO3$p.value,
           rich.RA.no_oC$p.value, rich.RA.no_tP$p.value, rich.RA.no_NH4$p.value,
           rich.RA.no_H2O$p.value, rich.RA.no_FungDiv$p.value, rich.RA.no_FungiRich$p.value,
           rich.RA.no_quibit$p.value, rich.RA.no_PlantDiv$p.value, rich.RA.no_PlantRich$p.value)

estimate = c(rich.RA.no_pH$estimate, rich.RA.no_OM$estimate, rich.RA.no_NO3$estimate,
             rich.RA.no_oC$estimate, rich.RA.no_tP$estimate, rich.RA.no_NH4$estimate,
             rich.RA.no_H2O$estimate, rich.RA.no_FungDiv$estimate, rich.RA.no_FungiRich$estimate,
             rich.RA.no_quibit$estimate, rich.RA.no_PlantDiv$estimate, rich.RA.no_PlantRich$estimate)

corr_RA.no_rich <- data.frame(variable, signif, estimate)
corr_RA.no_rich$group <- "RA Richness (Native Around)"

# Diversity
# Run correlations 
# Degraded Apron
div.RA.da_pH = cor.test(x = RA_degradedapron$RA_div, y = RA_degradedapron$pH, method = "spearman", exact=FALSE)
div.RA.da_OM = cor.test(x = RA_degradedapron$RA_div, y = RA_degradedapron$OM, method = "spearman", exact=FALSE)
div.RA.da_NO3 = cor.test(x = RA_degradedapron$RA_div, y = RA_degradedapron$NO3, method = "spearman", exact=FALSE)
div.RA.da_oC = cor.test(x = RA_degradedapron$RA_div, y = RA_degradedapron$oC, method = "spearman", exact=FALSE)
div.RA.da_tP = cor.test(x = RA_degradedapron$RA_div, y = RA_degradedapron$tP, method = "spearman", exact=FALSE)
div.RA.da_NH4 = cor.test(x = RA_degradedapron$RA_div, y = RA_degradedapron$NH4, method = "spearman", exact=FALSE)
div.RA.da_H2O = cor.test(x = RA_degradedapron$RA_div, y = RA_degradedapron$H2O, method = "spearman", exact=FALSE)
div.RA.da_FungDiv = cor.test(x = RA_degradedapron$RA_div, y = RA_degradedapron$FungDiv, method = "spearman", exact=FALSE)
div.RA.da_FungiRich = cor.test(x = RA_degradedapron$RA_div, y = RA_degradedapron$FungiRich, method = "spearman", exact=FALSE)
div.RA.da_quibit = cor.test(x = RA_degradedapron$RA_div, y = RA_degradedapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit")
signif = c(div.RA.da_pH$p.value, div.RA.da_OM$p.value, div.RA.da_NO3$p.value,
           div.RA.da_oC$p.value, div.RA.da_tP$p.value, div.RA.da_NH4$p.value,
           div.RA.da_H2O$p.value, div.RA.da_FungDiv$p.value, div.RA.da_FungiRich$p.value,
           div.RA.da_quibit$p.value)

estimate = c(div.RA.da_pH$estimate, div.RA.da_OM$estimate, div.RA.da_NO3$estimate,
             div.RA.da_oC$estimate, div.RA.da_tP$estimate, div.RA.da_NH4$estimate,
             div.RA.da_H2O$estimate, div.RA.da_FungDiv$estimate, div.RA.da_FungiRich$estimate,
             div.RA.da_quibit$estimate)

corr_RA.da_div <- data.frame(variable, signif, estimate)
corr_RA.da_div$group <- "RA Diversity (Degraded Apron)"

# Run correlations 
# Degraded Around
div.RA.do_pH = cor.test(x = RA_degradedaround$RA_div, y = RA_degradedaround$pH, method = "spearman", exact=FALSE)
div.RA.do_OM = cor.test(x = RA_degradedaround$RA_div, y = RA_degradedaround$OM, method = "spearman", exact=FALSE)
div.RA.do_NO3 = cor.test(x = RA_degradedaround$RA_div, y = RA_degradedaround$NO3, method = "spearman", exact=FALSE)
div.RA.do_oC = cor.test(x = RA_degradedaround$RA_div, y = RA_degradedaround$oC, method = "spearman", exact=FALSE)
div.RA.do_tP = cor.test(x = RA_degradedaround$RA_div, y = RA_degradedaround$tP, method = "spearman", exact=FALSE)
div.RA.do_NH4 = cor.test(x = RA_degradedaround$RA_div, y = RA_degradedaround$NH4, method = "spearman", exact=FALSE)
div.RA.do_H2O = cor.test(x = RA_degradedaround$RA_div, y = RA_degradedaround$H2O, method = "spearman", exact=FALSE)
div.RA.do_FungDiv = cor.test(x = RA_degradedaround$RA_div, y = RA_degradedaround$FungDiv, method = "spearman", exact=FALSE)
div.RA.do_FungiRich = cor.test(x = RA_degradedaround$RA_div, y = RA_degradedaround$FungiRich, method = "spearman", exact=FALSE)
div.RA.do_quibit = cor.test(x = RA_degradedaround$RA_div, y = RA_degradedaround$Qubit, method = "spearman", exact=FALSE)
div.RA.do_PlantDiv = cor.test(x = RA_degradedaround$RA_div, y = RA_degradedaround$PlantDiv, method = "spearman", exact=FALSE)
div.RA.do_PlantRich = cor.test(x = RA_degradedaround$RA_div, y = RA_degradedaround$PlantRich, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit", "Plant Diversity", "Plant Richness")
signif = c(div.RA.do_pH$p.value, div.RA.do_OM$p.value, div.RA.do_NO3$p.value,
           div.RA.do_oC$p.value, div.RA.do_tP$p.value, div.RA.do_NH4$p.value,
           div.RA.do_H2O$p.value, div.RA.do_FungDiv$p.value, div.RA.do_FungiRich$p.value,
           div.RA.do_quibit$p.value, div.RA.do_PlantDiv$p.value, div.RA.do_PlantRich$p.value)

estimate = c(div.RA.do_pH$estimate, div.RA.do_OM$estimate, div.RA.do_NO3$estimate,
             div.RA.do_oC$estimate, div.RA.do_tP$estimate, div.RA.do_NH4$estimate,
             div.RA.do_H2O$estimate, div.RA.do_FungDiv$estimate, div.RA.do_FungiRich$estimate,
             div.RA.do_quibit$estimate, div.RA.do_PlantDiv$estimate, div.RA.do_PlantRich$estimate)

corr_RA.do_div <- data.frame(variable, signif, estimate)
corr_RA.do_div$group <- "RA Diversity (Degraded Around)"

# Run correlations 
# Native Apron
div.RA.na_pH = cor.test(x = RA_nativeapron$RA_div, y = RA_nativeapron$pH, method = "spearman", exact=FALSE)
div.RA.na_OM = cor.test(x = RA_nativeapron$RA_div, y = RA_nativeapron$OM, method = "spearman", exact=FALSE)
div.RA.na_NO3 = cor.test(x = RA_nativeapron$RA_div, y = RA_nativeapron$NO3, method = "spearman", exact=FALSE)
div.RA.na_oC = cor.test(x = RA_nativeapron$RA_div, y = RA_nativeapron$oC, method = "spearman", exact=FALSE)
div.RA.na_tP = cor.test(x = RA_nativeapron$RA_div, y = RA_nativeapron$tP, method = "spearman", exact=FALSE)
div.RA.na_NH4 = cor.test(x = RA_nativeapron$RA_div, y = RA_nativeapron$NH4, method = "spearman", exact=FALSE)
div.RA.na_H2O = cor.test(x = RA_nativeapron$RA_div, y = RA_nativeapron$H2O, method = "spearman", exact=FALSE)
div.RA.na_FungDiv = cor.test(x = RA_nativeapron$RA_div, y = RA_nativeapron$FungDiv, method = "spearman", exact=FALSE)
div.RA.na_FungiRich = cor.test(x = RA_nativeapron$RA_div, y = RA_nativeapron$FungiRich, method = "spearman", exact=FALSE)
div.RA.na_quibit = cor.test(x = RA_nativeapron$RA_div, y = RA_nativeapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit")
signif = c(div.RA.na_pH$p.value, div.RA.na_OM$p.value, div.RA.na_NO3$p.value,
           div.RA.na_oC$p.value, div.RA.na_tP$p.value, div.RA.na_NH4$p.value,
           div.RA.na_H2O$p.value, div.RA.na_FungDiv$p.value, div.RA.na_FungiRich$p.value,
           div.RA.na_quibit$p.value)

estimate = c(div.RA.na_pH$estimate, div.RA.na_OM$estimate, div.RA.na_NO3$estimate,
             div.RA.na_oC$estimate, div.RA.na_tP$estimate, div.RA.na_NH4$estimate,
             div.RA.na_H2O$estimate, div.RA.na_FungDiv$estimate, div.RA.na_FungiRich$estimate,
             div.RA.na_quibit$estimate)

corr_RA.na_div <- data.frame(variable, signif, estimate)
corr_RA.na_div$group <- "RA Diversity (Native Apron)"

# Run correlations 
# Native Around
div.RA.no_pH = cor.test(x = RA_nativearound$RA_div, y = RA_nativearound$pH, method = "spearman", exact=FALSE)
div.RA.no_OM = cor.test(x = RA_nativearound$RA_div, y = RA_nativearound$OM, method = "spearman", exact=FALSE)
div.RA.no_NO3 = cor.test(x = RA_nativearound$RA_div, y = RA_nativearound$NO3, method = "spearman", exact=FALSE)
div.RA.no_oC = cor.test(x = RA_nativearound$RA_div, y = RA_nativearound$oC, method = "spearman", exact=FALSE)
div.RA.no_tP = cor.test(x = RA_nativearound$RA_div, y = RA_nativearound$tP, method = "spearman", exact=FALSE)
div.RA.no_NH4 = cor.test(x = RA_nativearound$RA_div, y = RA_nativearound$NH4, method = "spearman", exact=FALSE)
div.RA.no_H2O = cor.test(x = RA_nativearound$RA_div, y = RA_nativearound$H2O, method = "spearman", exact=FALSE)
div.RA.no_FungDiv = cor.test(x = RA_nativearound$RA_div, y = RA_nativearound$FungDiv, method = "spearman", exact=FALSE)
div.RA.no_FungiRich = cor.test(x = RA_nativearound$RA_div, y = RA_nativearound$FungiRich, method = "spearman", exact=FALSE)
div.RA.no_quibit = cor.test(x = RA_nativearound$RA_div, y = RA_nativearound$Qubit, method = "spearman", exact=FALSE)
div.RA.no_PlantDiv = cor.test(x = RA_nativearound$RA_div, y = RA_nativearound$PlantDiv, method = "spearman", exact=FALSE)
div.RA.no_PlantRich = cor.test(x = RA_nativearound$RA_div, y = RA_nativearound$PlantRich, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", 
             "Quibit", "Plant Diversity", "Plant Richness")
signif = c(div.RA.no_pH$p.value, div.RA.no_OM$p.value, div.RA.no_NO3$p.value,
           div.RA.no_oC$p.value, div.RA.no_tP$p.value, div.RA.no_NH4$p.value,
           div.RA.no_H2O$p.value, div.RA.no_FungDiv$p.value, div.RA.no_FungiRich$p.value,
           div.RA.no_quibit$p.value, div.RA.no_PlantDiv$p.value, div.RA.no_PlantRich$p.value)

estimate = c(div.RA.no_pH$estimate, div.RA.no_OM$estimate, div.RA.no_NO3$estimate,
             div.RA.no_oC$estimate, div.RA.no_tP$estimate, div.RA.no_NH4$estimate,
             div.RA.no_H2O$estimate, div.RA.no_FungDiv$estimate, div.RA.no_FungiRich$estimate,
             div.RA.no_quibit$estimate, div.RA.no_PlantDiv$estimate, div.RA.no_PlantRich$estimate)

corr_RA.no_div <- data.frame(variable, signif, estimate)
corr_RA.no_div$group <- "RA Diversity (Native Around)"

# calc nem richness in each sample for major groups
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

# make sub tables for each TRT 

OM_degradedapron = OM_map %>%
  filter(TRT == "DegradedApron")
OM_degradedaround  = OM_map %>%
  filter(TRT == "DegradedAround")
OM_nativeapron = OM_map %>%
  filter(TRT == "NativeApron")
OM_nativearound = OM_map %>%
  filter(TRT == "NativeAround")

mean(OM_nativeapron$OM_div)
sd(OM_nativeapron$OM_div)
mean(OM_nativearound$OM_div)
sd(OM_nativearound$OM_div)
mean(OM_degradedapron$OM_div)
sd(OM_degradedapron$OM_div)
mean(OM_degradedaround$OM_div)
sd(OM_degradedaround$OM_div)

mean(OM_nativeapron$OM_rich)
sd(OM_nativeapron$OM_rich)
mean(OM_nativearound$OM_rich)
sd(OM_nativearound$OM_rich)
mean(OM_degradedapron$OM_rich)
sd(OM_degradedapron$OM_rich)
mean(OM_degradedaround$OM_rich)
sd(OM_degradedaround$OM_rich)

# Richness
# Run correlations 
# Degraded Apron
rich.OM.da_pH = cor.test(x = OM_degradedapron$OM_rich, y = OM_degradedapron$pH, method = "spearman", exact=FALSE)
rich.OM.da_OM = cor.test(x = OM_degradedapron$OM_rich, y = OM_degradedapron$OM, method = "spearman", exact=FALSE)
rich.OM.da_NO3 = cor.test(x = OM_degradedapron$OM_rich, y = OM_degradedapron$NO3, method = "spearman", exact=FALSE)
rich.OM.da_oC = cor.test(x = OM_degradedapron$OM_rich, y = OM_degradedapron$oC, method = "spearman", exact=FALSE)
rich.OM.da_tP = cor.test(x = OM_degradedapron$OM_rich, y = OM_degradedapron$tP, method = "spearman", exact=FALSE)
rich.OM.da_NH4 = cor.test(x = OM_degradedapron$OM_rich, y = OM_degradedapron$NH4, method = "spearman", exact=FALSE)
rich.OM.da_H2O = cor.test(x = OM_degradedapron$OM_rich, y = OM_degradedapron$H2O, method = "spearman", exact=FALSE)
rich.OM.da_FungDiv = cor.test(x = OM_degradedapron$OM_rich, y = OM_degradedapron$FungDiv, method = "spearman", exact=FALSE)
rich.OM.da_FungiRich = cor.test(x = OM_degradedapron$OM_rich, y = OM_degradedapron$FungiRich, method = "spearman", exact=FALSE)
rich.OM.da_BactDiv = cor.test(x = OM_degradedapron$OM_rich, y = OM_degradedapron$BactDiv, method = "spearman", exact=FALSE)
rich.OM.da_BactRich = cor.test(x = OM_degradedapron$OM_rich, y = OM_degradedapron$BactRich, method = "spearman", exact=FALSE)
rich.OM.da_quibit = cor.test(x = OM_degradedapron$OM_rich, y = OM_degradedapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for OM model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", "Bacterial Diversity", "Bacterial Richness",  
             "Quibit")
signif = c(rich.OM.da_pH$p.value, rich.OM.da_OM$p.value, rich.OM.da_NO3$p.value,
           rich.OM.da_oC$p.value, rich.OM.da_tP$p.value, rich.OM.da_NH4$p.value,
           rich.OM.da_H2O$p.value, rich.OM.da_FungDiv$p.value, rich.OM.da_FungiRich$p.value, rich.OM.da_BactDiv$p.value, rich.OM.da_BactRich$p.value,
           rich.OM.da_quibit$p.value)

estimate = c(rich.OM.da_pH$estimate, rich.OM.da_OM$estimate, rich.OM.da_NO3$estimate,
             rich.OM.da_oC$estimate, rich.OM.da_tP$estimate, rich.OM.da_NH4$estimate,
             rich.OM.da_H2O$estimate, rich.OM.da_FungDiv$estimate, rich.OM.da_FungiRich$estimate, rich.OM.da_BactDiv$estimate, rich.OM.da_BactRich$estimate,
             rich.OM.da_quibit$estimate)

corr_OM.da_rich <- data.frame(variable, signif, estimate)
corr_OM.da_rich$group <- "OM Richness (Degraded Apron)"

# Run correlations 
# Degraded Around
rich.OM.do_pH = cor.test(x = OM_degradedaround$OM_rich, y = OM_degradedaround$pH, method = "spearman", exact=FALSE)
rich.OM.do_OM = cor.test(x = OM_degradedaround$OM_rich, y = OM_degradedaround$OM, method = "spearman", exact=FALSE)
rich.OM.do_NO3 = cor.test(x = OM_degradedaround$OM_rich, y = OM_degradedaround$NO3, method = "spearman", exact=FALSE)
rich.OM.do_oC = cor.test(x = OM_degradedaround$OM_rich, y = OM_degradedaround$oC, method = "spearman", exact=FALSE)
rich.OM.do_tP = cor.test(x = OM_degradedaround$OM_rich, y = OM_degradedaround$tP, method = "spearman", exact=FALSE)
rich.OM.do_NH4 = cor.test(x = OM_degradedaround$OM_rich, y = OM_degradedaround$NH4, method = "spearman", exact=FALSE)
rich.OM.do_H2O = cor.test(x = OM_degradedaround$OM_rich, y = OM_degradedaround$H2O, method = "spearman", exact=FALSE)
rich.OM.do_FungDiv = cor.test(x = OM_degradedaround$OM_rich, y = OM_degradedaround$FungDiv, method = "spearman", exact=FALSE)
rich.OM.do_FungiRich = cor.test(x = OM_degradedaround$OM_rich, y = OM_degradedaround$FungiRich, method = "spearman", exact=FALSE)
rich.OM.do_BactDiv = cor.test(x = OM_degradedaround$OM_rich, y = OM_degradedaround$BactDiv, method = "spearman", exact=FALSE)
rich.OM.do_BactRich = cor.test(x = OM_degradedaround$OM_rich, y = OM_degradedaround$BactRich, method = "spearman", exact=FALSE)
rich.OM.do_quibit = cor.test(x = OM_degradedaround$OM_rich, y = OM_degradedaround$Qubit, method = "spearman", exact=FALSE)

# make dataframe for OM model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", "Bacterial Diversity", "Bacterial Richness",  
             "Quibit")
signif = c(rich.OM.do_pH$p.value, rich.OM.do_OM$p.value, rich.OM.do_NO3$p.value,
           rich.OM.do_oC$p.value, rich.OM.do_tP$p.value, rich.OM.do_NH4$p.value,
           rich.OM.do_H2O$p.value, rich.OM.do_FungDiv$p.value, rich.OM.do_FungiRich$p.value, rich.OM.do_BactDiv$p.value, rich.OM.do_BactRich$p.value,
           rich.OM.do_quibit$p.value)

estimate = c(rich.OM.do_pH$estimate, rich.OM.do_OM$estimate, rich.OM.do_NO3$estimate,
             rich.OM.do_oC$estimate, rich.OM.do_tP$estimate, rich.OM.do_NH4$estimate,
             rich.OM.do_H2O$estimate, rich.OM.do_FungDiv$estimate, rich.OM.do_FungiRich$estimate, rich.OM.do_BactDiv$estimate, rich.OM.do_BactRich$estimate,
             rich.OM.do_quibit$estimate)

corr_OM.do_rich <- data.frame(variable, signif, estimate)
corr_OM.do_rich$group <- "OM Richness (Degraded Around)"

# Run correlations 
# Native Apron
rich.OM.na_pH = cor.test(x = OM_nativeapron$OM_rich, y = OM_nativeapron$pH, method = "spearman", exact=FALSE)
rich.OM.na_OM = cor.test(x = OM_nativeapron$OM_rich, y = OM_nativeapron$OM, method = "spearman", exact=FALSE)
rich.OM.na_NO3 = cor.test(x = OM_nativeapron$OM_rich, y = OM_nativeapron$NO3, method = "spearman", exact=FALSE)
rich.OM.na_oC = cor.test(x = OM_nativeapron$OM_rich, y = OM_nativeapron$oC, method = "spearman", exact=FALSE)
rich.OM.na_tP = cor.test(x = OM_nativeapron$OM_rich, y = OM_nativeapron$tP, method = "spearman", exact=FALSE)
rich.OM.na_NH4 = cor.test(x = OM_nativeapron$OM_rich, y = OM_nativeapron$NH4, method = "spearman", exact=FALSE)
rich.OM.na_H2O = cor.test(x = OM_nativeapron$OM_rich, y = OM_nativeapron$H2O, method = "spearman", exact=FALSE)
rich.OM.na_FungDiv = cor.test(x = OM_nativeapron$OM_rich, y = OM_nativeapron$FungDiv, method = "spearman", exact=FALSE)
rich.OM.na_FungiRich = cor.test(x = OM_nativeapron$OM_rich, y = OM_nativeapron$FungiRich, method = "spearman", exact=FALSE)
rich.OM.na_BactDiv = cor.test(x = OM_nativeapron$OM_rich, y = OM_nativeapron$BactDiv, method = "spearman", exact=FALSE)
rich.OM.na_BactRich = cor.test(x = OM_nativeapron$OM_rich, y = OM_nativeapron$BactRich, method = "spearman", exact=FALSE)
rich.OM.na_quibit = cor.test(x = OM_nativeapron$OM_rich, y = OM_nativeapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for OM model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", "Bacterial Diversity", "Bacterial Richness",  
             "Quibit")
signif = c(rich.OM.na_pH$p.value, rich.OM.na_OM$p.value, rich.OM.na_NO3$p.value,
           rich.OM.na_oC$p.value, rich.OM.na_tP$p.value, rich.OM.na_NH4$p.value,
           rich.OM.na_H2O$p.value, rich.OM.na_FungDiv$p.value, rich.OM.na_FungiRich$p.value, rich.OM.na_BactDiv$p.value, rich.OM.na_BactRich$p.value,
           rich.OM.na_quibit$p.value)

estimate = c(rich.OM.na_pH$estimate, rich.OM.na_OM$estimate, rich.OM.na_NO3$estimate,
             rich.OM.na_oC$estimate, rich.OM.na_tP$estimate, rich.OM.na_NH4$estimate,
             rich.OM.na_H2O$estimate, rich.OM.na_FungDiv$estimate, rich.OM.na_FungiRich$estimate, rich.OM.na_BactDiv$estimate, rich.OM.na_BactRich$estimate,
             rich.OM.na_quibit$estimate)

corr_OM.na_rich <- data.frame(variable, signif, estimate)
corr_OM.na_rich$group <- "OM Richness (Native Apron)"

# Run correlations 
# Native Around
rich.OM.no_pH = cor.test(x = OM_nativearound$OM_rich, y = OM_nativearound$pH, method = "spearman", exact=FALSE)
rich.OM.no_OM = cor.test(x = OM_nativearound$OM_rich, y = OM_nativearound$OM, method = "spearman", exact=FALSE)
rich.OM.no_NO3 = cor.test(x = OM_nativearound$OM_rich, y = OM_nativearound$NO3, method = "spearman", exact=FALSE)
rich.OM.no_oC = cor.test(x = OM_nativearound$OM_rich, y = OM_nativearound$oC, method = "spearman", exact=FALSE)
rich.OM.no_tP = cor.test(x = OM_nativearound$OM_rich, y = OM_nativearound$tP, method = "spearman", exact=FALSE)
rich.OM.no_NH4 = cor.test(x = OM_nativearound$OM_rich, y = OM_nativearound$NH4, method = "spearman", exact=FALSE)
rich.OM.no_H2O = cor.test(x = OM_nativearound$OM_rich, y = OM_nativearound$H2O, method = "spearman", exact=FALSE)
rich.OM.no_FungDiv = cor.test(x = OM_nativearound$OM_rich, y = OM_nativearound$FungDiv, method = "spearman", exact=FALSE)
rich.OM.no_FungiRich = cor.test(x = OM_nativearound$OM_rich, y = OM_nativearound$FungiRich, method = "spearman", exact=FALSE)
rich.OM.no_BactDiv = cor.test(x = OM_nativearound$OM_rich, y = OM_nativearound$BactDiv, method = "spearman", exact=FALSE)
rich.OM.no_BactRich = cor.test(x = OM_nativearound$OM_rich, y = OM_nativearound$BactRich, method = "spearman", exact=FALSE)
rich.OM.no_quibit = cor.test(x = OM_nativearound$OM_rich, y = OM_nativearound$Qubit, method = "spearman", exact=FALSE)

# make dataframe for OM model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", "Bacterial Diversity", "Bacterial Richness",  
             "Quibit")
signif = c(rich.OM.no_pH$p.value, rich.OM.no_OM$p.value, rich.OM.no_NO3$p.value,
           rich.OM.no_oC$p.value, rich.OM.no_tP$p.value, rich.OM.no_NH4$p.value,
           rich.OM.no_H2O$p.value, rich.OM.no_FungDiv$p.value, rich.OM.no_FungiRich$p.value, rich.OM.no_BactDiv$p.value, rich.OM.no_BactRich$p.value,
           rich.OM.no_quibit$p.value)

estimate = c(rich.OM.no_pH$estimate, rich.OM.no_OM$estimate, rich.OM.no_NO3$estimate,
             rich.OM.no_oC$estimate, rich.OM.no_tP$estimate, rich.OM.no_NH4$estimate,
             rich.OM.no_H2O$estimate, rich.OM.no_FungDiv$estimate, rich.OM.no_FungiRich$estimate, rich.OM.no_BactDiv$estimate, rich.OM.no_BactRich$estimate,
             rich.OM.no_quibit$estimate)

corr_OM.no_rich <- data.frame(variable, signif, estimate)
corr_OM.no_rich$group <- "OM Richness (Native Around)"

# Diversity
# Run correlations 
# Degraded Apron
div.OM.da_pH = cor.test(x = OM_degradedapron$OM_div, y = OM_degradedapron$pH, method = "spearman", exact=FALSE)
div.OM.da_OM = cor.test(x = OM_degradedapron$OM_div, y = OM_degradedapron$OM, method = "spearman", exact=FALSE)
div.OM.da_NO3 = cor.test(x = OM_degradedapron$OM_div, y = OM_degradedapron$NO3, method = "spearman", exact=FALSE)
div.OM.da_oC = cor.test(x = OM_degradedapron$OM_div, y = OM_degradedapron$oC, method = "spearman", exact=FALSE)
div.OM.da_tP = cor.test(x = OM_degradedapron$OM_div, y = OM_degradedapron$tP, method = "spearman", exact=FALSE)
div.OM.da_NH4 = cor.test(x = OM_degradedapron$OM_div, y = OM_degradedapron$NH4, method = "spearman", exact=FALSE)
div.OM.da_H2O = cor.test(x = OM_degradedapron$OM_div, y = OM_degradedapron$H2O, method = "spearman", exact=FALSE)
div.OM.da_FungDiv = cor.test(x = OM_degradedapron$OM_div, y = OM_degradedapron$FungDiv, method = "spearman", exact=FALSE)
div.OM.da_FungiRich = cor.test(x = OM_degradedapron$OM_div, y = OM_degradedapron$FungiRich, method = "spearman", exact=FALSE)
div.OM.da_BactDiv = cor.test(x = OM_degradedapron$OM_div, y = OM_degradedapron$BactDiv, method = "spearman", exact=FALSE)
div.OM.da_BactRich = cor.test(x = OM_degradedapron$OM_div, y = OM_degradedapron$BactRich, method = "spearman", exact=FALSE)
div.OM.da_quibit = cor.test(x = OM_degradedapron$OM_div, y = OM_degradedapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for OM model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", "Bacterial Diversity", "Bacterial Richness",  
             "Quibit")
signif = c(div.OM.da_pH$p.value, div.OM.da_OM$p.value, div.OM.da_NO3$p.value,
           div.OM.da_oC$p.value, div.OM.da_tP$p.value, div.OM.da_NH4$p.value,
           div.OM.da_H2O$p.value, div.OM.da_FungDiv$p.value, div.OM.da_FungiRich$p.value, div.OM.da_BactDiv$p.value, div.OM.da_BactRich$p.value,
           div.OM.da_quibit$p.value)

estimate = c(div.OM.da_pH$estimate, div.OM.da_OM$estimate, div.OM.da_NO3$estimate,
             div.OM.da_oC$estimate, div.OM.da_tP$estimate, div.OM.da_NH4$estimate,
             div.OM.da_H2O$estimate, div.OM.da_FungDiv$estimate, div.OM.da_FungiRich$estimate, div.OM.da_BactDiv$estimate, div.OM.da_BactRich$estimate,
             div.OM.da_quibit$estimate)

corr_OM.da_div <- data.frame(variable, signif, estimate)
corr_OM.da_div$group <- "OM Diversity (Degraded Apron)"

# Run correlations 
# Degraded Around
div.OM.do_pH = cor.test(x = OM_degradedaround$OM_div, y = OM_degradedaround$pH, method = "spearman", exact=FALSE)
div.OM.do_OM = cor.test(x = OM_degradedaround$OM_div, y = OM_degradedaround$OM, method = "spearman", exact=FALSE)
div.OM.do_NO3 = cor.test(x = OM_degradedaround$OM_div, y = OM_degradedaround$NO3, method = "spearman", exact=FALSE)
div.OM.do_oC = cor.test(x = OM_degradedaround$OM_div, y = OM_degradedaround$oC, method = "spearman", exact=FALSE)
div.OM.do_tP = cor.test(x = OM_degradedaround$OM_div, y = OM_degradedaround$tP, method = "spearman", exact=FALSE)
div.OM.do_NH4 = cor.test(x = OM_degradedaround$OM_div, y = OM_degradedaround$NH4, method = "spearman", exact=FALSE)
div.OM.do_H2O = cor.test(x = OM_degradedaround$OM_div, y = OM_degradedaround$H2O, method = "spearman", exact=FALSE)
div.OM.do_FungDiv = cor.test(x = OM_degradedaround$OM_div, y = OM_degradedaround$FungDiv, method = "spearman", exact=FALSE)
div.OM.do_FungiRich = cor.test(x = OM_degradedaround$OM_div, y = OM_degradedaround$FungiRich, method = "spearman", exact=FALSE)
div.OM.do_BactDiv = cor.test(x = OM_degradedaround$OM_div, y = OM_degradedaround$BactDiv, method = "spearman", exact=FALSE)
div.OM.do_BactRich = cor.test(x = OM_degradedaround$OM_div, y = OM_degradedaround$BactRich, method = "spearman", exact=FALSE)
div.OM.do_quibit = cor.test(x = OM_degradedaround$OM_div, y = OM_degradedaround$Qubit, method = "spearman", exact=FALSE)

# make dataframe for OM model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", "Bacterial Diversity", "Bacterial Richness",  
             "Quibit")
signif = c(div.OM.do_pH$p.value, div.OM.do_OM$p.value, div.OM.do_NO3$p.value,
           div.OM.do_oC$p.value, div.OM.do_tP$p.value, div.OM.do_NH4$p.value,
           div.OM.do_H2O$p.value, div.OM.do_FungDiv$p.value, div.OM.do_FungiRich$p.value, div.OM.do_BactDiv$p.value, div.OM.do_BactRich$p.value,
           div.OM.do_quibit$p.value)

estimate = c(div.OM.do_pH$estimate, div.OM.do_OM$estimate, div.OM.do_NO3$estimate,
             div.OM.do_oC$estimate, div.OM.do_tP$estimate, div.OM.do_NH4$estimate,
             div.OM.do_H2O$estimate, div.OM.do_FungDiv$estimate, div.OM.do_FungiRich$estimate, div.OM.do_BactDiv$estimate, div.OM.do_BactRich$estimate,
             div.OM.do_quibit$estimate)

corr_OM.do_div <- data.frame(variable, signif, estimate)
corr_OM.do_div$group <- "OM Diversity (Degraded Around)"

# Run correlations 
# Native Apron
div.OM.na_pH = cor.test(x = OM_nativeapron$OM_div, y = OM_nativeapron$pH, method = "spearman", exact=FALSE)
div.OM.na_OM = cor.test(x = OM_nativeapron$OM_div, y = OM_nativeapron$OM, method = "spearman", exact=FALSE)
div.OM.na_NO3 = cor.test(x = OM_nativeapron$OM_div, y = OM_nativeapron$NO3, method = "spearman", exact=FALSE)
div.OM.na_oC = cor.test(x = OM_nativeapron$OM_div, y = OM_nativeapron$oC, method = "spearman", exact=FALSE)
div.OM.na_tP = cor.test(x = OM_nativeapron$OM_div, y = OM_nativeapron$tP, method = "spearman", exact=FALSE)
div.OM.na_NH4 = cor.test(x = OM_nativeapron$OM_div, y = OM_nativeapron$NH4, method = "spearman", exact=FALSE)
div.OM.na_H2O = cor.test(x = OM_nativeapron$OM_div, y = OM_nativeapron$H2O, method = "spearman", exact=FALSE)
div.OM.na_FungDiv = cor.test(x = OM_nativeapron$OM_div, y = OM_nativeapron$FungDiv, method = "spearman", exact=FALSE)
div.OM.na_FungiRich = cor.test(x = OM_nativeapron$OM_div, y = OM_nativeapron$FungiRich, method = "spearman", exact=FALSE)
div.OM.na_BactDiv = cor.test(x = OM_nativeapron$OM_div, y = OM_nativeapron$BactDiv, method = "spearman", exact=FALSE)
div.OM.na_BactRich = cor.test(x = OM_nativeapron$OM_div, y = OM_nativeapron$BactRich, method = "spearman", exact=FALSE)
div.OM.na_quibit = cor.test(x = OM_nativeapron$OM_div, y = OM_nativeapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for OM model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", "Bacterial Diversity", "Bacterial Richness",  
             "Quibit")
signif = c(div.OM.na_pH$p.value, div.OM.na_OM$p.value, div.OM.na_NO3$p.value,
           div.OM.na_oC$p.value, div.OM.na_tP$p.value, div.OM.na_NH4$p.value,
           div.OM.na_H2O$p.value, div.OM.na_FungDiv$p.value, div.OM.na_FungiRich$p.value, div.OM.na_BactDiv$p.value, div.OM.na_BactRich$p.value,
           div.OM.na_quibit$p.value)

estimate = c(div.OM.na_pH$estimate, div.OM.na_OM$estimate, div.OM.na_NO3$estimate,
             div.OM.na_oC$estimate, div.OM.na_tP$estimate, div.OM.na_NH4$estimate,
             div.OM.na_H2O$estimate, div.OM.na_FungDiv$estimate, div.OM.na_FungiRich$estimate, div.OM.na_BactDiv$estimate, div.OM.na_BactRich$estimate,
             div.OM.na_quibit$estimate)

corr_OM.na_div <- data.frame(variable, signif, estimate)
corr_OM.na_div$group <- "OM Diversity (Native Apron)"

# Run correlations 
# Native Around
div.OM.no_pH = cor.test(x = OM_nativearound$OM_div, y = OM_nativearound$pH, method = "spearman", exact=FALSE)
div.OM.no_OM = cor.test(x = OM_nativearound$OM_div, y = OM_nativearound$OM, method = "spearman", exact=FALSE)
div.OM.no_NO3 = cor.test(x = OM_nativearound$OM_div, y = OM_nativearound$NO3, method = "spearman", exact=FALSE)
div.OM.no_oC = cor.test(x = OM_nativearound$OM_div, y = OM_nativearound$oC, method = "spearman", exact=FALSE)
div.OM.no_tP = cor.test(x = OM_nativearound$OM_div, y = OM_nativearound$tP, method = "spearman", exact=FALSE)
div.OM.no_NH4 = cor.test(x = OM_nativearound$OM_div, y = OM_nativearound$NH4, method = "spearman", exact=FALSE)
div.OM.no_H2O = cor.test(x = OM_nativearound$OM_div, y = OM_nativearound$H2O, method = "spearman", exact=FALSE)
div.OM.no_FungDiv = cor.test(x = OM_nativearound$OM_div, y = OM_nativearound$FungDiv, method = "spearman", exact=FALSE)
div.OM.no_FungiRich = cor.test(x = OM_nativearound$OM_div, y = OM_nativearound$FungiRich, method = "spearman", exact=FALSE)
div.OM.no_BactDiv = cor.test(x = OM_nativearound$OM_div, y = OM_nativearound$BactDiv, method = "spearman", exact=FALSE)
div.OM.no_BactRich = cor.test(x = OM_nativearound$OM_div, y = OM_nativearound$BactRich, method = "spearman", exact=FALSE)
div.OM.no_quibit = cor.test(x = OM_nativearound$OM_div, y = OM_nativearound$Qubit, method = "spearman", exact=FALSE)

# make dataframe for OM model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Fungal Diversity", "FungalRichness", "Bacterial Diversity", "Bacterial Richness",  
             "Quibit")
signif = c(div.OM.no_pH$p.value, div.OM.no_OM$p.value, div.OM.no_NO3$p.value,
           div.OM.no_oC$p.value, div.OM.no_tP$p.value, div.OM.no_NH4$p.value,
           div.OM.no_H2O$p.value, div.OM.no_FungDiv$p.value, div.OM.no_FungiRich$p.value, div.OM.no_BactDiv$p.value, div.OM.no_BactRich$p.value,
           div.OM.no_quibit$p.value)

estimate = c(div.OM.no_pH$estimate, div.OM.no_OM$estimate, div.OM.no_NO3$estimate,
             div.OM.no_oC$estimate, div.OM.no_tP$estimate, div.OM.no_NH4$estimate,
             div.OM.no_H2O$estimate, div.OM.no_FungDiv$estimate, div.OM.no_FungiRich$estimate, div.OM.no_BactDiv$estimate, div.OM.no_BactRich$estimate,
             div.OM.no_quibit$estimate)

corr_OM.no_div <- data.frame(variable, signif, estimate)
corr_OM.no_div$group <- "OM Diversity (Native Around)" 

# calc nem richness in each sample for major groups
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

# make sub tables for each TRT 

PP_degradedapron = PP_map %>%
  filter(TRT == "DegradedApron")
PP_degradedaround  = PP_map %>%
  filter(TRT == "DegradedAround")
PP_nativeapron = PP_map %>%
  filter(TRT == "NativeApron")
PP_nativearound = PP_map %>%
  filter(TRT == "NativeAround")

mean(PP_nativeapron$PP_div)
sd(PP_nativeapron$PP_div)
mean(PP_nativearound$PP_div)
sd(PP_nativearound$PP_div)
mean(PP_degradedapron$PP_div)
sd(PP_degradedapron$PP_div)
mean(PP_degradedaround$PP_div)
sd(PP_degradedaround$PP_div)

mean(PP_nativeapron$PP_rich)
sd(PP_nativeapron$PP_rich)
mean(PP_nativearound$PP_rich)
sd(PP_nativearound$PP_rich)
mean(PP_degradedapron$PP_rich)
sd(PP_degradedapron$PP_rich)
mean(PP_degradedaround$PP_rich)
sd(PP_degradedaround$PP_rich)

# Richness
# Run correlations 
# Degraded Apron
rich.PP.da_pH = cor.test(x = PP_degradedapron$PP_rich, y = PP_degradedapron$pH, method = "spearman", exact=FALSE)
rich.PP.da_OM = cor.test(x = PP_degradedapron$PP_rich, y = PP_degradedapron$OM, method = "spearman", exact=FALSE)
rich.PP.da_NO3 = cor.test(x = PP_degradedapron$PP_rich, y = PP_degradedapron$NO3, method = "spearman", exact=FALSE)
rich.PP.da_oC = cor.test(x = PP_degradedapron$PP_rich, y = PP_degradedapron$oC, method = "spearman", exact=FALSE)
rich.PP.da_tP = cor.test(x = PP_degradedapron$PP_rich, y = PP_degradedapron$tP, method = "spearman", exact=FALSE)
rich.PP.da_NH4 = cor.test(x = PP_degradedapron$PP_rich, y = PP_degradedapron$NH4, method = "spearman", exact=FALSE)
rich.PP.da_H2O = cor.test(x = PP_degradedapron$PP_rich, y = PP_degradedapron$H2O, method = "spearman", exact=FALSE)
rich.PP.da_quibit = cor.test(x = PP_degradedapron$PP_rich, y = PP_degradedapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Quibit")
signif = c(rich.PP.da_pH$p.value, rich.PP.da_OM$p.value, rich.PP.da_NO3$p.value,
           rich.PP.da_oC$p.value, rich.PP.da_tP$p.value, rich.PP.da_NH4$p.value,
           rich.PP.da_H2O$p.value, 
           rich.PP.da_quibit$p.value)

estimate = c(rich.PP.da_pH$estimate, rich.PP.da_OM$estimate, rich.PP.da_NO3$estimate,
             rich.PP.da_oC$estimate, rich.PP.da_tP$estimate, rich.PP.da_NH4$estimate,
             rich.PP.da_H2O$estimate, 
             rich.PP.da_quibit$estimate)

corr_PP.da_rich <- data.frame(variable, signif, estimate)
corr_PP.da_rich$group <- "PP Richness (Degraded Apron)"

# Run correlations 
# Degraded Around
rich.PP.do_pH = cor.test(x = PP_degradedaround$PP_rich, y = PP_degradedaround$pH, method = "spearman", exact=FALSE)
rich.PP.do_OM = cor.test(x = PP_degradedaround$PP_rich, y = PP_degradedaround$OM, method = "spearman", exact=FALSE)
rich.PP.do_NO3 = cor.test(x = PP_degradedaround$PP_rich, y = PP_degradedaround$NO3, method = "spearman", exact=FALSE)
rich.PP.do_oC = cor.test(x = PP_degradedaround$PP_rich, y = PP_degradedaround$oC, method = "spearman", exact=FALSE)
rich.PP.do_tP = cor.test(x = PP_degradedaround$PP_rich, y = PP_degradedaround$tP, method = "spearman", exact=FALSE)
rich.PP.do_NH4 = cor.test(x = PP_degradedaround$PP_rich, y = PP_degradedaround$NH4, method = "spearman", exact=FALSE)
rich.PP.do_H2O = cor.test(x = PP_degradedaround$PP_rich, y = PP_degradedaround$H2O, method = "spearman", exact=FALSE)
rich.PP.do_quibit = cor.test(x = PP_degradedaround$PP_rich, y = PP_degradedaround$Qubit, method = "spearman", exact=FALSE)
rich.PP.do_PlantDiv = cor.test(x = PP_degradedaround$PP_rich, y = PP_degradedaround$PlantDiv, method = "spearman", exact=FALSE)
rich.PP.do_PlantRich = cor.test(x = PP_degradedaround$PP_rich, y = PP_degradedaround$PlantRich, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Quibit", "Plant Diversity", "Plant Richness")
signif = c(rich.PP.do_pH$p.value, rich.PP.do_OM$p.value, rich.PP.do_NO3$p.value,
           rich.PP.do_oC$p.value, rich.PP.do_tP$p.value, rich.PP.do_NH4$p.value,
           rich.PP.do_H2O$p.value, 
           rich.PP.do_quibit$p.value, rich.PP.do_PlantDiv$p.value, rich.PP.do_PlantRich$p.value)

estimate = c(rich.PP.do_pH$estimate, rich.PP.do_OM$estimate, rich.PP.do_NO3$estimate,
             rich.PP.do_oC$estimate, rich.PP.do_tP$estimate, rich.PP.do_NH4$estimate,
             rich.PP.do_H2O$estimate, 
             rich.PP.do_quibit$estimate, rich.PP.do_PlantDiv$estimate, rich.PP.do_PlantRich$estimate)

corr_PP.do_rich <- data.frame(variable, signif, estimate)
corr_PP.do_rich$group <- "PP Richness (Degraded Around)"

# Run correlations 
# Native Apron
rich.PP.na_pH = cor.test(x = PP_nativeapron$PP_rich, y = PP_nativeapron$pH, method = "spearman", exact=FALSE)
rich.PP.na_OM = cor.test(x = PP_nativeapron$PP_rich, y = PP_nativeapron$OM, method = "spearman", exact=FALSE)
rich.PP.na_NO3 = cor.test(x = PP_nativeapron$PP_rich, y = PP_nativeapron$NO3, method = "spearman", exact=FALSE)
rich.PP.na_oC = cor.test(x = PP_nativeapron$PP_rich, y = PP_nativeapron$oC, method = "spearman", exact=FALSE)
rich.PP.na_tP = cor.test(x = PP_nativeapron$PP_rich, y = PP_nativeapron$tP, method = "spearman", exact=FALSE)
rich.PP.na_NH4 = cor.test(x = PP_nativeapron$PP_rich, y = PP_nativeapron$NH4, method = "spearman", exact=FALSE)
rich.PP.na_H2O = cor.test(x = PP_nativeapron$PP_rich, y = PP_nativeapron$H2O, method = "spearman", exact=FALSE)
rich.PP.na_quibit = cor.test(x = PP_nativeapron$PP_rich, y = PP_nativeapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Quibit")
signif = c(rich.PP.na_pH$p.value, rich.PP.na_OM$p.value, rich.PP.na_NO3$p.value,
           rich.PP.na_oC$p.value, rich.PP.na_tP$p.value, rich.PP.na_NH4$p.value,
           rich.PP.na_H2O$p.value, 
           rich.PP.na_quibit$p.value)

estimate = c(rich.PP.na_pH$estimate, rich.PP.na_OM$estimate, rich.PP.na_NO3$estimate,
             rich.PP.na_oC$estimate, rich.PP.na_tP$estimate, rich.PP.na_NH4$estimate,
             rich.PP.na_H2O$estimate, 
             rich.PP.na_quibit$estimate)

corr_PP.na_rich <- data.frame(variable, signif, estimate)
corr_PP.na_rich$group <- "PP Richness (Native Apron)"

# Run correlations 
# Native Around
rich.PP.no_pH = cor.test(x = PP_nativearound$PP_rich, y = PP_nativearound$pH, method = "spearman", exact=FALSE)
rich.PP.no_OM = cor.test(x = PP_nativearound$PP_rich, y = PP_nativearound$OM, method = "spearman", exact=FALSE)
rich.PP.no_NO3 = cor.test(x = PP_nativearound$PP_rich, y = PP_nativearound$NO3, method = "spearman", exact=FALSE)
rich.PP.no_oC = cor.test(x = PP_nativearound$PP_rich, y = PP_nativearound$oC, method = "spearman", exact=FALSE)
rich.PP.no_tP = cor.test(x = PP_nativearound$PP_rich, y = PP_nativearound$tP, method = "spearman", exact=FALSE)
rich.PP.no_NH4 = cor.test(x = PP_nativearound$PP_rich, y = PP_nativearound$NH4, method = "spearman", exact=FALSE)
rich.PP.no_H2O = cor.test(x = PP_nativearound$PP_rich, y = PP_nativearound$H2O, method = "spearman", exact=FALSE)
rich.PP.no_quibit = cor.test(x = PP_nativearound$PP_rich, y = PP_nativearound$Qubit, method = "spearman", exact=FALSE)
rich.PP.no_PlantDiv = cor.test(x = PP_nativearound$PP_rich, y = PP_nativearound$PlantDiv, method = "spearman", exact=FALSE)
rich.PP.no_PlantRich = cor.test(x = PP_nativearound$PP_rich, y = PP_nativearound$PlantRich, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture",  
             "Quibit", "Plant Diversity", "Plant Richness")
signif = c(rich.PP.no_pH$p.value, rich.PP.no_OM$p.value, rich.PP.no_NO3$p.value,
           rich.PP.no_oC$p.value, rich.PP.no_tP$p.value, rich.PP.no_NH4$p.value,
           rich.PP.no_H2O$p.value, 
           rich.PP.no_quibit$p.value, rich.PP.no_PlantDiv$p.value, rich.PP.no_PlantRich$p.value)

estimate = c(rich.PP.no_pH$estimate, rich.PP.no_OM$estimate, rich.PP.no_NO3$estimate,
             rich.PP.no_oC$estimate, rich.PP.no_tP$estimate, rich.PP.no_NH4$estimate,
             rich.PP.no_H2O$estimate, 
             rich.PP.no_quibit$estimate, rich.PP.no_PlantDiv$estimate, rich.PP.no_PlantRich$estimate)

corr_PP.no_rich <- data.frame(variable, signif, estimate)
corr_PP.no_rich$group <- "PP Richness (Native Around)"

# Diversity
# Run correlations 
# Degraded Apron
div.PP.da_pH = cor.test(x = PP_degradedapron$PP_div, y = PP_degradedapron$pH, method = "spearman", exact=FALSE)
div.PP.da_OM = cor.test(x = PP_degradedapron$PP_div, y = PP_degradedapron$OM, method = "spearman", exact=FALSE)
div.PP.da_NO3 = cor.test(x = PP_degradedapron$PP_div, y = PP_degradedapron$NO3, method = "spearman", exact=FALSE)
div.PP.da_oC = cor.test(x = PP_degradedapron$PP_div, y = PP_degradedapron$oC, method = "spearman", exact=FALSE)
div.PP.da_tP = cor.test(x = PP_degradedapron$PP_div, y = PP_degradedapron$tP, method = "spearman", exact=FALSE)
div.PP.da_NH4 = cor.test(x = PP_degradedapron$PP_div, y = PP_degradedapron$NH4, method = "spearman", exact=FALSE)
div.PP.da_H2O = cor.test(x = PP_degradedapron$PP_div, y = PP_degradedapron$H2O, method = "spearman", exact=FALSE)
div.PP.da_quibit = cor.test(x = PP_degradedapron$PP_div, y = PP_degradedapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Quibit")
signif = c(div.PP.da_pH$p.value, div.PP.da_OM$p.value, div.PP.da_NO3$p.value,
           div.PP.da_oC$p.value, div.PP.da_tP$p.value, div.PP.da_NH4$p.value,
           div.PP.da_H2O$p.value, 
           div.PP.da_quibit$p.value)

estimate = c(div.PP.da_pH$estimate, div.PP.da_OM$estimate, div.PP.da_NO3$estimate,
             div.PP.da_oC$estimate, div.PP.da_tP$estimate, div.PP.da_NH4$estimate,
             div.PP.da_H2O$estimate, 
             div.PP.da_quibit$estimate)

corr_PP.da_div <- data.frame(variable, signif, estimate)
corr_PP.da_div$group <- "PP Diversity (Degraded Apron)"

# Run correlations 
# Degraded Around
div.PP.do_pH = cor.test(x = PP_degradedaround$PP_div, y = PP_degradedaround$pH, method = "spearman", exact=FALSE)
div.PP.do_OM = cor.test(x = PP_degradedaround$PP_div, y = PP_degradedaround$OM, method = "spearman", exact=FALSE)
div.PP.do_NO3 = cor.test(x = PP_degradedaround$PP_div, y = PP_degradedaround$NO3, method = "spearman", exact=FALSE)
div.PP.do_oC = cor.test(x = PP_degradedaround$PP_div, y = PP_degradedaround$oC, method = "spearman", exact=FALSE)
div.PP.do_tP = cor.test(x = PP_degradedaround$PP_div, y = PP_degradedaround$tP, method = "spearman", exact=FALSE)
div.PP.do_NH4 = cor.test(x = PP_degradedaround$PP_div, y = PP_degradedaround$NH4, method = "spearman", exact=FALSE)
div.PP.do_H2O = cor.test(x = PP_degradedaround$PP_div, y = PP_degradedaround$H2O, method = "spearman", exact=FALSE)
div.PP.do_quibit = cor.test(x = PP_degradedaround$PP_div, y = PP_degradedaround$Qubit, method = "spearman", exact=FALSE)
div.PP.do_PlantDiv = cor.test(x = PP_degradedaround$PP_div, y = PP_degradedaround$PlantDiv, method = "spearman", exact=FALSE)
div.PP.do_PlantRich = cor.test(x = PP_degradedaround$PP_div, y = PP_degradedaround$PlantRich, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Quibit", "Plant Diversity", "Plant Richness")
signif = c(div.PP.do_pH$p.value, div.PP.do_OM$p.value, div.PP.do_NO3$p.value,
           div.PP.do_oC$p.value, div.PP.do_tP$p.value, div.PP.do_NH4$p.value,
           div.PP.do_H2O$p.value, 
           div.PP.do_quibit$p.value, div.PP.do_PlantDiv$p.value, div.PP.do_PlantRich$p.value)

estimate = c(div.PP.do_pH$estimate, div.PP.do_OM$estimate, div.PP.do_NO3$estimate,
             div.PP.do_oC$estimate, div.PP.do_tP$estimate, div.PP.do_NH4$estimate,
             div.PP.do_H2O$estimate, 
             div.PP.do_quibit$estimate, div.PP.do_PlantDiv$estimate, div.PP.do_PlantRich$estimate)

corr_PP.do_div <- data.frame(variable, signif, estimate)
corr_PP.do_div$group <- "PP Diversity (Degraded Around)"

# Run correlations 
# Native Apron
div.PP.na_pH = cor.test(x = PP_nativeapron$PP_div, y = PP_nativeapron$pH, method = "spearman", exact=FALSE)
div.PP.na_OM = cor.test(x = PP_nativeapron$PP_div, y = PP_nativeapron$OM, method = "spearman", exact=FALSE)
div.PP.na_NO3 = cor.test(x = PP_nativeapron$PP_div, y = PP_nativeapron$NO3, method = "spearman", exact=FALSE)
div.PP.na_oC = cor.test(x = PP_nativeapron$PP_div, y = PP_nativeapron$oC, method = "spearman", exact=FALSE)
div.PP.na_tP = cor.test(x = PP_nativeapron$PP_div, y = PP_nativeapron$tP, method = "spearman", exact=FALSE)
div.PP.na_NH4 = cor.test(x = PP_nativeapron$PP_div, y = PP_nativeapron$NH4, method = "spearman", exact=FALSE)
div.PP.na_H2O = cor.test(x = PP_nativeapron$PP_div, y = PP_nativeapron$H2O, method = "spearman", exact=FALSE)
div.PP.na_quibit = cor.test(x = PP_nativeapron$PP_div, y = PP_nativeapron$Qubit, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Quibit")
signif = c(div.PP.na_pH$p.value, div.PP.na_OM$p.value, div.PP.na_NO3$p.value,
           div.PP.na_oC$p.value, div.PP.na_tP$p.value, div.PP.na_NH4$p.value,
           div.PP.na_H2O$p.value, 
           div.PP.na_quibit$p.value)

estimate = c(div.PP.na_pH$estimate, div.PP.na_OM$estimate, div.PP.na_NO3$estimate,
             div.PP.na_oC$estimate, div.PP.na_tP$estimate, div.PP.na_NH4$estimate,
             div.PP.na_H2O$estimate, 
             div.PP.na_quibit$estimate)

corr_PP.na_div <- data.frame(variable, signif, estimate)
corr_PP.na_div$group <- "PP Diversity (Native Apron)"

# Run correlations 
# Native Around
div.PP.no_pH = cor.test(x = PP_nativearound$PP_div, y = PP_nativearound$pH, method = "spearman", exact=FALSE)
div.PP.no_OM = cor.test(x = PP_nativearound$PP_div, y = PP_nativearound$OM, method = "spearman", exact=FALSE)
div.PP.no_NO3 = cor.test(x = PP_nativearound$PP_div, y = PP_nativearound$NO3, method = "spearman", exact=FALSE)
div.PP.no_oC = cor.test(x = PP_nativearound$PP_div, y = PP_nativearound$oC, method = "spearman", exact=FALSE)
div.PP.no_tP = cor.test(x = PP_nativearound$PP_div, y = PP_nativearound$tP, method = "spearman", exact=FALSE)
div.PP.no_NH4 = cor.test(x = PP_nativearound$PP_div, y = PP_nativearound$NH4, method = "spearman", exact=FALSE)
div.PP.no_H2O = cor.test(x = PP_nativearound$PP_div, y = PP_nativearound$H2O, method = "spearman", exact=FALSE)
div.PP.no_quibit = cor.test(x = PP_nativearound$PP_div, y = PP_nativearound$Qubit, method = "spearman", exact=FALSE)
div.PP.no_PlantDiv = cor.test(x = PP_nativearound$PP_div, y = PP_nativearound$PlantDiv, method = "spearman", exact=FALSE)
div.PP.no_PlantRich = cor.test(x = PP_nativearound$PP_div, y = PP_nativearound$PlantRich, method = "spearman", exact=FALSE)

# make dataframe for RA model results
variable = c("pH", "Organic Matter",  "NO3",
             "Organic C", "Total P", "NH4",
             "Soil Moisture", 
             "Quibit", "Plant Diversity", "Plant Richness")
signif = c(div.PP.no_pH$p.value, div.PP.no_OM$p.value, div.PP.no_NO3$p.value,
           div.PP.no_oC$p.value, div.PP.no_tP$p.value, div.PP.no_NH4$p.value,
           div.PP.no_H2O$p.value, 
           div.PP.no_quibit$p.value, div.PP.no_PlantDiv$p.value, div.PP.no_PlantRich$p.value)

estimate = c(div.PP.no_pH$estimate, div.PP.no_OM$estimate, div.PP.no_NO3$estimate,
             div.PP.no_oC$estimate, div.PP.no_tP$estimate, div.PP.no_NH4$estimate,
             div.PP.no_H2O$estimate, 
             div.PP.no_quibit$estimate, div.PP.no_PlantDiv$estimate, div.PP.no_PlantRich$estimate)

corr_PP.no_div <- data.frame(variable, signif, estimate)
corr_PP.no_div$group <- "PP Diversity (Native Around)"

drivers <- rbind(corr_BF.da_div, corr_BF.do_div, corr_BF.na_div, corr_BF.no_div,
                 corr_BF.da_rich, corr_BF.do_rich, corr_BF.na_rich, corr_BF.no_rich,
                 corr_FF.da_div, corr_FF.do_div, corr_FF.na_div, corr_FF.no_div,
                 corr_FF.da_rich, corr_FF.do_rich, corr_FF.na_rich, corr_FF.no_rich,
                 corr_RA.da_div, corr_RA.do_div, corr_RA.na_div, corr_RA.no_div,
                 corr_RA.da_rich, corr_RA.do_rich, corr_RA.na_rich, corr_RA.no_rich,
                 corr_OM.da_div, corr_OM.do_div, corr_OM.na_div, corr_OM.no_div,
                 corr_OM.da_rich, corr_OM.do_rich, corr_OM.na_rich, corr_OM.no_rich,
                 corr_PP.da_div, corr_PP.do_div, corr_PP.na_div, corr_PP.no_div,
                 corr_PP.da_rich, corr_PP.do_rich, corr_PP.na_rich, corr_PP.no_rich)

write.csv(drivers, "drivers_of_trophic_diversity.csv")

figure <- read.csv("drivers_of_trophic_diversity_f.csv", row.names = 1)

figure_s <- figure %>%
  filter(variable != "Quibit")

fake_plant_1 <- c("Plant Diversity", 0, 0, "RA", "Diversity", "Degraded", "Apron")
fake_plant_2 <- c("Plant Diversity", 0, 0, "RA", "Diversity", "Native", "Apron")
fake_plant_3 <- c("Plant Diversity", 0, 0, "RA", "Richness", "Degraded", "Apron")
fake_plant_4 <- c("Plant Diversity", 0, 0, "RA", "Richness", "Native", "Apron")
fake_plant_5 <- c("Plant Richness", 0, 0, "RA", "Diversity", "Degraded", "Apron")
fake_plant_6 <- c("Plant Richness", 0, 0, "RA", "Diversity", "Native", "Apron")
fake_plant_7 <- c("Plant Richness", 0, 0, "RA", "Richness", "Degraded", "Apron")
fake_plant_8 <- c("Plant Richness", 0, 0, "RA", "Richness", "Native", "Apron")

fake_plant_9 <- c("Plant Diversity", 0, 0, "PP", "Diversity", "Degraded", "Apron")
fake_plant_10 <- c("Plant Diversity", 0, 0, "PP", "Diversity", "Native", "Apron")
fake_plant_11 <- c("Plant Diversity", 0, 0, "PP", "Richness", "Degraded", "Apron")
fake_plant_12 <- c("Plant Diversity", 0, 0, "PP", "Richness", "Native", "Apron")
fake_plant_13 <- c("Plant Richness", 0, 0, "PP", "Diversity", "Degraded", "Apron")
fake_plant_14 <- c("Plant Richness", 0, 0, "PP", "Diversity", "Native", "Apron")
fake_plant_15 <- c("Plant Richness", 0, 0, "PP", "Richness", "Degraded", "Apron")
fake_plant_16 <- c("Plant Richness", 0, 0, "PP", "Richness", "Native", "Apron")


figure_s_a <- rbind(figure_s, fake_plant_1, fake_plant_2, fake_plant_3,
                    fake_plant_4, fake_plant_5, fake_plant_6, fake_plant_7, fake_plant_8,
                    fake_plant_9, fake_plant_10, fake_plant_11, fake_plant_12, fake_plant_13,
                    fake_plant_14, fake_plant_15, fake_plant_16)
figure_s_a$signif <- as.numeric(figure_s_a$signif)
figure_s_a$estimate <- as.numeric(figure_s_a$estimate)
figure_s_a$habitat = factor(figure_s_a$habitat, 
                    levels = c("Native", "Degraded"))
figure_s_a$activity = factor(figure_s_a$activity, 
                            levels = c("Around", "Apron"))
figure_s_a$trophicgroup = factor(figure_s_a$trophicgroup, 
                        levels = c("BF", "FF", "OM", "PP", "RA"))
figure_s_a$variable = factor(figure_s_a$variable, 
                             levels = c("Soil Moisture", "pH", "NO3", "NH4", 
                                        "Organic C", "Total P", "Organic Matter",  
                                        "Fungal Diversity", "FungalRichness", "Bacterial Diversity", "Bacterial Richness",
                                        "Plant Diversity", "Plant Richness"))

richness_plot = figure_s_a %>%
  filter(metric == "Richness")
ggplot(data = richness_plot, aes(x = variable, y = estimate, fill = activity)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(trophicgroup ~ habitat) +
  scale_fill_manual(values = c("green4","khaki")) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) 

div_plot = figure_s_a %>%
  filter(metric == "Diversity")
ggplot(data = div_plot, aes(x = variable, y = estimate, fill = activity)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(trophicgroup ~ habitat) +
  scale_fill_manual(values = c("green4","khaki")) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) 

