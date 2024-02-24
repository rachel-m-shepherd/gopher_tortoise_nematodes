# abundance of Nematode families

library(mctoolsr)
library(tidyverse)

# input data
tax_table_nem = 'nem_asv_gt.txt'
map_fp = 'metadata_fall_gt.txt'

input_nem = load_taxa_table(tax_table_nem, map_fp)

input_nem_ra = convert_to_relative_abundances(input_nem)
input_nem_ra$map_loaded$sampleID <- row.names(input_nem_ra$map_loaded)
input_nem_ra = filter_data(input_nem_ra, filter_cat = 'sampleID',
                           filter_vals = c("F.10.O_F", "F.7.A_F", "S.3.A_F","S.6.O_F"))

map = input_nem_ra$map_loaded
## BY Treatment
table_nem_level4_fp = summarize_taxonomy(input_nem_ra, level = 4,
                                         report_higher_tax = F)
total_nem = table_nem_level4_fp %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sampleID") %>%
  left_join(map) %>%
  pivot_longer(cols=c(2:62), names_to = "Family", values_to = 'abund') %>%
  as.data.frame() 


total_nem$Family = factor(total_nem$Family, 
                                 levels = c("Tylenchidae", "Tobrilidae",        
                                            "Anatonchidae", "Actinolaimidae",     
                                            "Nygolaimidae" , "Discolaimidae",      
                                            "Trischistomatidae", "Ironidae",           
                                            "Mylonchulidae",       "Mononchidae",        
                                            "Anguinidae",          "Tylenchulidae",      
                                            "Hoplolaimidae",       "Merliniidae",        
                                            "Criconematidae",      "Meloidogynidae",     
                                            "Pratylenchidae",     "Belondiridae",       
                                            "Trichodoridae",    "Telotylenchidae",    
                                            "Belonolaimidae",   "Longidoridae",       
                                            "Hemicycliophoridae",  "Qudsianematidae",    
                                            "Aporcelaimidae",  "Nordiidae",          
                                            "Dorylaimidae", "Chrysonematidae",    
                                            "Mydonomidae",         "Thornenematidae",    
                                            "Aphelenchoididae",    "Diphtherophoridae",  
                                            "Leptonchidae",        "Aphelenchidae",      
                                            "Tylencholaimellidae", "Tylencholaimidae",   
                                            "Sphaerulariidae",     "Allantonematidae",   
                                            "Cephalobidae",        "Bastianiidae",       
                                            "Prismatolaimidae",    "Alaimidae",          
                                            "Monhysteridae",       "Chronogastridae",    
                                            "Microlaimidae",       "Comesomatidae",      
                                            "Rhabditidae",         "Axonolaimidae",      
                                            "Plectidae",           "Drilonematidae",     
                                            "Rhabdolaimidae",      "Diplopeltidae",      
                                            "Haliplectidae",       "Bunonematidae",      
                                            "Neodiplogasteridae",  "Steinernematidae",   
                                            "Isolaimiidae",        "Ohridiidae",         
                                            "Odontolaimidae",      "Diplogasteridae",    
                                            "Cryptonchidae",       "Mermithidae",        
                                             "Trichostrongylidae",  "Cyatholaimidae"))


total_nem$sampleID = factor(total_nem$sampleID, 
                                    levels = c("S.1.O_F","S.2.O_F","S.3.O_F","S.4.O_F","S.5.O_F","S.7.O_F","S.8.O_F","S.9.O_F","S.10.O_F","S.11.O_F","S.12.O_F","S.13.O_F","S.14.O_F","S.15.O_F",
                                               "S.1.A_F","S.2.A_F","S.4.A_F","S.5.A_F","S.6.A_F","S.7.A_F","S.8.A_F","S.9.A_F","S.10.A_F","S.11.A_F","S.12.A_F","S.13.A_F","S.14.A_F","S.15.A_F",
                                               "F.1.O_F","F.2.O_F","F.3.O_F","F.4.O_F","F.5.O_F","F.6.O_F","F.7.O_F","F.8.O_F", "F.9.O_F","F.11.O_F","F.12.O_F","F.13.O_F","F.14.O_F","F.15.O_F",
                                               "F.1.A_F","F.2.A_F", "F.3.A_F", "F.4.A_F","F.5.A_F","F.6.A_F", "F.8.A_F", "F.9.A_F","F.10.A_F", "F.11.A_F","F.12.A_F","F.13.A_F", "F.14.A_F", "F.15.A_F"))



ggplot(total_nem, aes(x=sampleID, y=abund)) + 
  geom_bar(aes(fill = Family), position = "stack", stat="identity") +
  scale_fill_manual(values = c("black",	"#7393B3",	"#00FFFF",	"#F0FFFF",	"#0000FF",	"#0096FF",	"#CCCCFF","#96DED1",	"#40E0D0",	"#1434A4",	
                               "#880808",	"#EE4B2B",	"#6E260E",	"#D22B2B",	"#F88379",	"#9A2A2A",	"#C04000",	"#FAA0A0",	"#E30B5C",	"#FA8072",	"#A42A04",	"#E34234",	"#913831",	
                               "#454B1B",	"#AAFF00",	"#AFE1AF",	"#50C878",	"#7CFC00",	"#4F7942",	"#C9CC3F",
                               "#BF40BF",	"#AA336A",	"#5D3FD3", "#E6E6FA",	"#E0B0FF",	"#7F00FF",	"grey",
                               "#FFBF00",	"#FBCEB1",	"#CD7F32",	"#FFC000",	"#FFEA00",	"#F28C28",	"#E97451",	"#FAD5A5",	"#DAA06D",	"#F4BB44",	"#FFAA33",	"#B87333",	"#FF7518",	"#CC5500",	"#FFFF8F",	"#F0E68C",	"#FADA5E",	"#FFFAA0",	"#F8DE7E",	"#FDDA0D",	"#E1C16E",
                               "#FF00FF",	"#FF69B4",	"#FFB6C1"))+
  labs(
    x = NULL,
    y = "Relative Abundance of Nematode Families",
    fill = NULL
  ) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(fill = guide_legend(nrow = 61))

trt_nem = total_nem %>%
  group_by(TRT, Family) %>% 
  summarise(mean_RA = mean(abund))

trt_nem$Family = factor(trt_nem$Family, 
                          levels = c("Tylenchidae", "Tobrilidae",        
                                     "Anatonchidae", "Actinolaimidae",     
                                     "Nygolaimidae" , "Discolaimidae",      
                                     "Trischistomatidae", "Ironidae",           
                                     "Mylonchulidae",       "Mononchidae",        
                                     "Anguinidae",          "Tylenchulidae",      
                                     "Hoplolaimidae",       "Merliniidae",        
                                     "Criconematidae",      "Meloidogynidae",     
                                     "Pratylenchidae",     "Belondiridae",       
                                     "Trichodoridae",    "Telotylenchidae",    
                                     "Belonolaimidae",   "Longidoridae",       
                                     "Hemicycliophoridae",  "Qudsianematidae",    
                                     "Aporcelaimidae",  "Nordiidae",          
                                     "Dorylaimidae", "Chrysonematidae",    
                                     "Mydonomidae",         "Thornenematidae",    
                                     "Aphelenchoididae",    "Diphtherophoridae",  
                                     "Leptonchidae",        "Aphelenchidae",      
                                     "Tylencholaimellidae", "Tylencholaimidae",   
                                     "Sphaerulariidae",     "Allantonematidae",   
                                     "Cephalobidae",        "Bastianiidae",       
                                     "Prismatolaimidae",    "Alaimidae",          
                                     "Monhysteridae",       "Chronogastridae",    
                                     "Microlaimidae",       "Comesomatidae",      
                                     "Rhabditidae",         "Axonolaimidae",      
                                     "Plectidae",           "Drilonematidae",     
                                     "Rhabdolaimidae",      "Diplopeltidae",      
                                     "Haliplectidae",       "Bunonematidae",      
                                     "Neodiplogasteridae",  "Steinernematidae",   
                                     "Isolaimiidae",        "Ohridiidae",         
                                     "Odontolaimidae",      "Diplogasteridae",    
                                     "Cryptonchidae",       "Mermithidae",        
                                     "Trichostrongylidae",  "Cyatholaimidae"))


trt_nem$TRT = factor(trt_nem$TRT, 
                            levels = c("NativeAround", "NativeApron", "DegradedAround", "DegradedApron"))



ggplot(trt_nem, aes(x=TRT, y=mean_RA)) + 
  geom_bar(aes(fill = Family), position = "stack", stat="identity") +
  scale_fill_manual(values = c("black",	"#7393B3",	"#00FFFF",	"#F0FFFF",	"#0000FF",	"#0096FF",	"#CCCCFF","#96DED1",	"#40E0D0",	"#1434A4",	
                               "#880808",	"#EE4B2B",	"#6E260E",	"#D22B2B",	"#F88379",	"#9A2A2A",	"#C04000",	"#FAA0A0",	"#E30B5C",	"#FA8072",	"#A42A04",	"#E34234",	"#913831",	
                               "#454B1B",	"#AAFF00",	"#AFE1AF",	"#50C878",	"#7CFC00",	"#4F7942",	"#C9CC3F",
                               "#BF40BF",	"#AA336A",	"#5D3FD3", "#E6E6FA",	"#E0B0FF",	"#7F00FF",	"grey",
                               "#FFBF00",	"#FBCEB1",	"#CD7F32",	"#FFC000",	"#FFEA00",	"#F28C28",	"#E97451",	"#FAD5A5",	"#DAA06D",	"#F4BB44",	"#FFAA33",	"#B87333",	"#FF7518",	"#CC5500",	"#FFFF8F",	"#F0E68C",	"#FADA5E",	"#FFFAA0",	"#F8DE7E",	"#FDDA0D",	"#E1C16E",
                               "#FF00FF",	"#FF69B4",	"#FFB6C1"))+
  labs(
    x = NULL,
    y = "Relative Abundance of Nematode Families",
    fill = NULL
  ) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(fill = guide_legend(nrow = 61))
