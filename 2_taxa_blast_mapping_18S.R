###### BLAST 18s GBIF differential abundant, excipients removed ######

##### set working directory and paths for data and figures ------

directory="2024-03_merged_with_june2023/asvs_single"


objects=c("table_merged","rep-seqs_merged") 

#from office
first_part="K:/ogden_grp/Conservation_Science/Conservation_Genetics/GENOME_PROJECTS/Carla/sequencing/sequencing_results/data_analyses"

path_work=paste(first_part,directory, sep="/")
path_work

setwd(path_work)

## Load packages, working paths and metadata ========

######function to load the libraries and install them/ load them if not installed ------

ipak <- function(pkg){
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  } #install BiocManager if needed 
  lapply(packages, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      tryCatch({
        install.packages(pkg) #install packages from CRAN
        library(pkg, character.only = TRUE) 
      }, error = function(e) {
        BiocManager::install(pkg, force = TRUE) # if package is not found in CRAN it installs it from Biocmanager. Force=TRUE install the package even if incompatibility of version of R. This needs to be treated with Caution though, as it is not warranty that the package will run smoothly if created in a more recent R version than the one used by the user
        library(pkg, character.only = TRUE)  
      })
    }
  })
  sapply(pkg, require, character.only = TRUE) #tells you if packages have been loaded
}



packages <- c("R.utils","Biostrings","ggplot2","openxlsx","reshape2","qiime2R",
              "stringr","dplyr","Polychrome","Rmisc","rgl","vegan","gridExtra",
              "gplots","Rmisc","rgbif","rworldmap","countrycode","ggmap","hexbin","viridis",
              "grattantheme")

ipak(packages)


# Load blast results -> result is a tsv table in blastn format 6: https://www.metagenomics.wiki/tools/blast/blastn-output-format-6

all_blast=read.table("blast_all/blast_18S_pe0.90_qc0.90_PR2_all_export/blast6.tsv",sep = "\t")
colnames(all_blast)=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

all_blast_ev=all_blast %>%
  dplyr::count(qseqid, sort = TRUE)

lineage_with_metadata=read.xlsx('K:/ogden_grp/Conservation_Science/Conservation_Genetics/GENOME_PROJECTS/Carla/sequencing/databases/pr2/pr2_version_5.0.0_merged.xlsx',sheet=1) 
lineage_with_metadata <- lineage_with_metadata[grep("Eukaryota",lineage_with_metadata$domain),] #keep only the eukaryotic sequences

dir_save=c("blast_deseq_ancomb","blast_manual") #since not ancombc results after removing excipients sequences, only use the deseq results
tables_query=c("/blast_deseq_ancomb/deseq_ancomb_select_18S","/blast_manual/manual_select_asvs_18S")
names_tablets=c("da_all","manual")


asvs_names=read.csv2("asvs_rep_names_18S.csv")

excipients_count_rar=read.csv("excipients_counts_rar_18S.csv", header=TRUE,sep = ";",check.names=FALSE)
colnames(excipients_count_rar)[1] <- "X"

for (h in 1:length(tables_query)) {
  print(tables_query[h])
  print(names_tablets[h])
  
  query_18S=read.xlsx(paste(path_work,tables_query[h],".xlsx",sep=""),sheet=1)
  query_18S=query_18S[!query_18S$Feature.ID %in% excipients_count_rar$X,]
  
  if (names_tablets[h]=="da_all"){
    query_18S =query_18S %>% mutate(log2FoldChange = coalesce(log2FoldChange, lfc_origin_countryEngland))
    query_18S$asv_country=ifelse(query_18S$log2FoldChange > 0, "E","T")}
  
  
  blast_table_18S=all_blast[all_blast$qseqid %in% query_18S$Feature.ID,]
  blast_table_18S=blast_table_18S[blast_table_18S$pident>0,]
  
  asv_ini=length(unique(query_18S$Feature.ID))
  
  query_18S=query_18S[query_18S$Feature.ID %in% blast_table_18S$qseqid,]
  
  asv_nohit=length(unique(query_18S$Feature.ID))
  
  blast_table_18S=blast_table_18S[blast_table_18S$sseqid %in% lineage_with_metadata$pr2_accession,]
  query_18S=query_18S[query_18S$Feature.ID %in% blast_table_18S$qseqid,]
  
  asv_inPR2=length(unique(query_18S$Feature.ID))
  
  feature_ids=query_18S$Feature.ID  
  
  best_hits_selection=list()
  
  for (i in 1:length(feature_ids)) {
    print(feature_ids[i])
    
    max_hits <- blast_table_18S %>%
      filter(qseqid==feature_ids[i]) %>%
      filter(pident==max(pident))
    
    best_hits_selection [[i]] <- max_hits
    
  }
  
  all_best_hits_selection=dplyr::bind_rows(best_hits_selection)
  
  best_hits_metadata=merge(all_best_hits_selection,lineage_with_metadata,by.x="sseqid",by.y="pr2_accession")
  
  #continue filtering
  best_hits_metadata=best_hits_metadata[!best_hits_metadata$qseqid %in% best_hits_plants$qseqid,] #remove plants for now
  best_hits_metadata <- best_hits_metadata[!grepl("_sp.",best_hits_metadata$species),]
  best_hits_metadata$species=gsub("_"," ",best_hits_metadata$species)
  
  length(unique(best_hits_metadata$qseqid))
  
  best_hits_metadata2 <- best_hits_metadata  %>% 
    dplyr::select(qseqid,pident,species) %>%
    group_by(qseqid,species) %>%
    dplyr::summarise(mean_pident=mean(pident)) %>%
    mutate(species=gsub('Gibberella intermedia','Fusarium fujikuroi', species)) %>% #appears as synonym in gbif, if not changing the command gives an error
    filter(species!="Homo sapiens")
  
   
  #select a minimum threshold for minimum percent identity
  threshold=97
  
  best_hits_metadata2=best_hits_metadata2[best_hits_metadata2$mean_pident>=threshold,]
  
  asv_id=unique(best_hits_metadata2$qseqid)
  
  
  write.csv(cbind(asv_ini=asv_ini, asv_blasthits=asv_nohit, asv_PR2=asv_inPR2, asv_pidenttrhe=length(asv_id)), paste(path_work,"/",dir_save[h],"/",names_tablets[h],"_remaining_asvs_after_filtering.csv",sep=""))
  
  
  all_gbif_results2=data.frame() # do it with a dataframe because couldnt manage to store element within element in list
  
  for (i in (1:length(asv_id))) {
    print(asv_id[i])
    
    species_to_query <- best_hits_metadata2$species[best_hits_metadata2$qseqid==asv_id[i]]
    mean_pident=best_hits_metadata2$mean_pident[best_hits_metadata2$qseqid==asv_id[i]]
    
    for (j in (1:length(species_to_query))) {
      print(species_to_query[j])
      
      name=name_backbone(species_to_query[j])$species
      
      gbif_results<-occ_data(
        scientificName=name, #taxon key was giving a mistake for a query so I searched by name instead
        hasGeospatialIssue = FALSE,
        occurrenceStatus = "PRESENT")
      
      gbif_results_df=gbif_results[['data']]
      
      if(is.null(gbif_results_df)) {
        print ("not found") } else { 
          
          gbif_results_df <- gbif_results_df[!(gbif_results_df$basisOfRecord %in% c("FOSSIL_SPECIMEN","LIVING_SPECIMEN","PRESERVED_SPECIMEN")),]
          gbif_results_df <- gbif_results_df %>% mutate(asv_id=asv_id[i],mean_pident=mean_pident[j])
          
          all_gbif_results2<-dplyr::bind_rows(all_gbif_results2, gbif_results_df)
        }
    }
  }
  
  unique(all_gbif_results2[colnames(all_gbif_results2) %in% c("asv_id","species","mean_pident")]) #to check identities in gbif.com
  
  all_gbif_results2$log2foldchange=query_18S$log2FoldChange[match(all_gbif_results2$asv_id,query_18S$Feature.ID)] #ignore ancombc results for now
  all_gbif_results2$asv_country=query_18S$asv_country[match(all_gbif_results2$asv_id,query_18S$Feature.ID)]
  
  assign(paste("all_gbif_results",names_tablets[h],sep="_"),all_gbif_results2) 
}



##### REPRESENTATION OF RESULTS --------------------

tables_gbif=c("all_gbif_results_da_all","all_gbif_results_manual")

countries_separation=c("E","T")
names_plot=c("England", "Thailand")


#create a base map for later ----------

world_map <- map_data("world")

base_world <- ggplot() + coord_fixed() +
  xlab("") + ylab("") +
  geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
               colour="grey80", fill="grey90") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "white"), legend.position="left",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())


globalfung_results_list <- list.files(path = paste(path_work,"/global_fungi_results/",sep = ""), pattern = "*.txt", full.names = TRUE)
globalfung_results <- do.call(rbind, lapply(globalfung_results_list, function(file) {
  data <- read.delim(file, header = TRUE, sep = "\t")  # adjust sep as needed
  data$species <- gsub(".txt","",basename(file))
  return(data)}))

theme_rep=theme(legend.position = "bottom",legend.key.size = unit(0.2, 'cm'),legend.key.width = unit(0.2, 'cm'),legend.text = element_text(size=10)) 

size_numbers=10


#### for deseq data ####

table=all_gbif_results_da_all
  
asv_species_viz <- table %>% 
   dplyr::select(asv_id,species) %>% 
   unique() %>% 
   mutate(asvs=asvs_names$rep_name[match(asv_id,asvs_names$Feature.ID)])
  
write.xlsx(asv_species_viz, paste(path_work,"/",dir_save[1],"/",names_tablets[1],"_species_represented.xlsx",sep="")) 

table2 <- table %>%
  select(asv_id,species,mean_pident, country) %>%
  filter(!is.na(country)) %>%
  distinct() %>%
  group_by(asv_id,species,mean_pident) %>%
  dplyr::summarize(countries = paste(country, collapse = ","))

query_18S=read.xlsx(paste(path_work,"/blast_deseq_ancomb/deseq_ancomb_select_18S.xlsx",sep=""),sheet=1)

table3=merge(table2,query_18S,by.x="asv_id",by.y="Feature.ID")
write.xlsx(table3, paste(path_work,"/",dir_save[1],"/",names_tablets[1],"_species_geoinfo_abundancesample.xlsx",sep="")) 


#represent taxa using GBIF and Globalgfungi for deseq data  --------

gbif_results_geoinfo <- table %>% 
  mutate(occurrence=rep(1,nrow(table))) %>% 
  distinct() %>%
  group_by(species,countryCode,country,decimalLatitude,decimalLongitude,continent,asv_id,asv_country,log2foldchange,mean_pident) %>% 
  dplyr::summarise(sum_occurrence = sum(occurrence)) %>% 
  dplyr::select(species,asv_country,asv_id,decimalLatitude,decimalLongitude,log2foldchange,mean_pident)

for (i in (1:length(countries_separation[1]))) {
    print(countries_separation[i])
    
    gbif_results_geoinfo_subset <- gbif_results_geoinfo %>% 
      filter(asv_country==countries_separation[i])
    
    gbif_results_geoinfo_subset2 <- gbif_results_geoinfo_subset %>%
      group_by(species,decimalLatitude,decimalLongitude) %>%
      na.omit() %>%
      dplyr::summarise(mean_log2foldchange=mean(round(log2foldchange,2)),mean_pident2=mean(round(mean_pident),2), n_asv=n_distinct(asv_id)) %>%
      mutate(species2=paste(species, " (",n_asv," ASVs, ",mean_pident2,"%, ",mean_log2foldchange,")",sep = "")) %>%
      filter(species!="Zea mays")
    
    levels_species <- na.omit(unique(gbif_results_geoinfo_subset2$species2[order(c(-gbif_results_geoinfo_subset2$mean_log2foldchange,-gbif_results_geoinfo_subset2$n_asv))]))

    viridis_colors <- viridis(length(levels_species), option = "F")

### using GBIF ------    
        
    map_data <- 
      base_world +
      geom_point(data=gbif_results_geoinfo_subset2, 
                 aes(x=decimalLongitude, y=decimalLatitude, fill=factor(species2, levels = levels_species)), colour="black", 
                 pch=21, size=5, alpha=I(0.7)) +
      scale_fill_manual(values = setNames(viridis_colors, levels_species), name="Possible species identity of da-ASVs") +
      ggtitle(paste(names_plot[i])) +
      theme_rep +
      guides(fill = guide_legend(ncol = 2, title.position="top"))
      
    
    assign(paste(countries_separation[i],"da_all_gbif_map",sep ="_"),map_data)
    grattan_save_pptx(map_data, filename=paste(path_work,"/pp_maps/da_all_",countries_separation[i],'_map.pptx',sep=""))
    ggsave(paste(path_work,"/blast_deseq_ancomb/da_all_",countries_separation[i],'_map.pdf',sep=""),map_data,device = cairo_pdf,width = 30, height = 15, units = "cm")
 
    
    subtitle_rep=paste("Distribution of",length(unique(gbif_results_geoinfo_subset2$species)),"possible species, represented by",length(unique(gbif_results_geoinfo_subset$asv_id)), "da-ASVs")
    
    map_data2 <- 
      base_world + 
      stat_bin_hex(data = gbif_results_geoinfo_subset2, aes(x = decimalLongitude, y = decimalLatitude), bins = 10, color = "black") +
      stat_bin_hex(data = gbif_results_geoinfo_subset2, aes(x = decimalLongitude, y = decimalLatitude, label=after_stat(count)), bins = 10, color = "black", geom="text",size=size_numbers) +
      scale_fill_gradient(low="white",high="yellow") +
      ggtitle(paste(names_plot[i])) +
      labs(subtitle=subtitle_rep) +
      theme(legend.position="right")
    
    assign(paste(countries_separation[i],"da_all_gbif_map2",sep ="_"),map_data2)
    grattan_save_pptx(map_data2, filename=paste(path_work,"/pp_maps/da_all_",countries_separation[i],'_map2.pptx',sep=""))
    ggsave(paste(path_work,"/blast_deseq_ancomb/da_all_",countries_separation[i],'_map2.pdf',sep=""),map_data2,device = cairo_pdf,
           width = 15, height = 10, units = "cm")
    ##### IMPORTANT -> da-ASV remaining of THAILAND match to Pichia sp. More than 100 species within this genus. Distribution of Pichia can be retrieved from GBIF
    

### using globalFungi ------    
    
    fung_select=unique(globalfung_results$species[(globalfung_results$species %in% gbif_results_geoinfo_subset2$species)])   

    globalfung_results2 <- globalfung_results[globalfung_results$species %in% fung_select,]
    globalfung_results2$species2=gbif_results_geoinfo_subset2$species2[match(globalfung_results2$species,gbif_results_geoinfo_subset2$species)]

    globalfung_map <- 
      base_world +
      geom_point(data=globalfung_results2, aes(x = longitude, y = latitude, fill=factor(species2, levels = levels_species)), colour="black",pch=21,size=3, alpha=I(0.7)) +  
      scale_fill_manual(values = setNames(viridis_colors, levels_species), name="Possible species identity of ASVs")+
      theme_rep +
      guides(fill = guide_legend(ncol = 2, title.position="top"))
    
    grattan_save_pptx(globalfung_map, filename=paste(path_work,"/pp_maps/",countries_separation[i],"_da_all_globalfungi_map.pptx",sep=""))          
    assign(paste(countries_separation[i],"da_all_globalfungi_map",sep ="_"),globalfung_map)
    ggsave(paste(path_work,"/blast_deseq_ancomb/da_all_",countries_separation[i],'_map_globalfung.pdf',sep=""),globalfung_map,device = cairo_pdf,
           width = 30, height = 15, units = "cm")
    
    subtitle_rep=paste("Distribution of",length(unique(globalfung_results2$species)),"possible species, represented by",length(unique(gbif_results_geoinfo$asv_id[gbif_results_geoinfo$species %in% globalfung_results2$species])), "da-ASVs")

    
    globalfung_results2 <- globalfung_results2[rep(1:nrow(globalfung_results2), globalfung_results2$abundances),]
    
    globalfung_map2 <- 
      base_world + 
      stat_bin_hex(data = globalfung_results2, aes(x = longitude, y = latitude), bins = 10, color = "black") +
      stat_bin_hex(data = globalfung_results2, aes(x = longitude, y = latitude, label=after_stat(count)), bins = 10, color = "black", geom="text") +
      scale_fill_gradient(low="white",high="yellow") +
      ggtitle(paste(names_plot[i])) +
      labs(subtitle=subtitle_rep) +
      theme(legend.position="right") 
  
  grattan_save_pptx(globalfung_map2, filename=paste(path_work,"/pp_maps/",countries_separation[i],"_da_all_globalfungi_map2.pptx",sep=""))          
  assign(paste(countries_separation[i],"da_all_globalfungi_map2",sep ="_"),globalfung_map2)  
  ggsave(paste(path_work,"/blast_deseq_ancomb/da_all_",countries_separation[i],'_map2_globalfung.pdf',sep=""),globalfung_map,device = cairo_pdf,
           width = 15, height = 10, units = "cm")

  rm(globalfung_map, globalfung_map2)
  
  }



#represent Pichia for Thailand using GBIF ####
  
pichia=read.csv("blast_deseq_ancomb/pichia.csv", sep="\t") #looking with the gbif function gives a different species than pichia, so directly downloaded from gbif

pichia2 <- pichia  %>%
dplyr::select(species,decimalLatitude,decimalLongitude,issue) 

map_pichia <- 
  base_world +
  geom_point(data=pichia2, 
             aes(x=decimalLongitude, y=decimalLatitude),fill="gold3",colour="black", 
             pch=21, size=5, alpha=I(0.7)) +
  ggtitle(paste(names_plot[i]))

grattan_save_pptx(map_pichia, filename=paste(path_work,"/pp_maps/da_all_da_all_pichia_map.pptx",sep=""))          
ggsave(paste(path_work,"/blast_deseq_ancomb/da_all_pichia_map.pdf",sep=""),map_pichia,device = cairo_pdf,
       width = 30, height = 15, units = "cm")

map_pichia2 <- 
  base_world + 
  stat_bin_hex(data = pichia2, aes(x = decimalLongitude, y = decimalLatitude), bins = 10, color = "black") +
  stat_bin_hex(data = pichia2, aes(x = decimalLongitude, y = decimalLatitude, label=after_stat(count)), bins = 10, color = "black", geom="text",size=size_numbers) +
  scale_fill_gradient(low="white",high="yellow") +
  theme(legend.position="right")  

grattan_save_pptx(map_pichia2, filename=paste(path_work,"/pp_maps/da_all_da_all_pichia_map2.pptx",sep=""))
ggsave(paste(path_work,"/blast_deseq_ancomb/da_all_pichia_map2.pdf",sep=""),map_pichia2,device = cairo_pdf,
       width = 15, height = 10, units = "cm")


#represent Pichia for Thailand using globalFungi ####

pichia_glofun=read.csv("blast_deseq_ancomb/global_fungi_results/pichia.txt", sep="\t")

map_pichia_glofun <- 
  base_world +
  geom_point(data=pichia_glofun, 
             aes(x = longitude, y = latitude),fill="gold3",colour="black", 
             pch=21, size=5, alpha=I(0.7)) +
  ggtitle(paste(names_plot[i]))

ggsave(paste(path_work,"/blast_deseq_ancomb/da_all_pichia_map_glofun.pdf",sep=""),map_pichia_glofun,device = cairo_pdf,
       width = 30, height = 15, units = "cm")


pichia_glofun2 <- pichia_glofun[rep(1:nrow(pichia_glofun), pichia_glofun$abundances),]

map_pichia_glofun2 <- 
  base_world + 
  stat_bin_hex(data = pichia_glofun2, aes(x = longitude, y = latitude), bins = 10, color = "black") +
  stat_bin_hex(data = pichia_glofun2, aes(x = longitude, y = latitude, label=after_stat(count)), bins = 10, color = "black", geom="text",size=size_numbers) +
  scale_fill_gradient(low="white",high="yellow") +
  theme(legend.position="right")  

ggsave(paste(path_work,"/blast_deseq_ancomb/da_all_pichia_map2_glofun2.pdf",sep=""),map_pichia_glofun2,device = cairo_pdf,
       width = 15, height = 10, units = "cm")


grattan_save_pptx(map_pichia_glofun, filename=paste(path_work,"/pp_maps/map_pichia_glofun.pptx",sep=""))
grattan_save_pptx(map_pichia_glofun2, filename=paste(path_work,"/pp_maps/map_pichia_glofun2.pptx",sep=""))

rm(map_pichia_glofun,map_pichia_glofun2)



#### for manual data ####

table = all_gbif_results_manual

asv_species_viz_manual <- table %>% 
  dplyr::select(asv_id,species) %>% 
  unique() %>% 
  mutate(asvs=asvs_names$rep_name[match(asv_id,asvs_names$Feature.ID)])

write.xlsx(asv_species_viz_manual, paste(path_work,"/",dir_save[2],"/",names_tablets[2],"_species_represented.xlsx",sep="")) 

table2 <- table %>%
  select(asv_id,species,mean_pident, country) %>%
  filter(!is.na(country)) %>%
  distinct() %>%
  group_by(asv_id,species,mean_pident) %>%
  dplyr::summarize(countries = paste(country, collapse = ","))

query_18S=read.xlsx(paste(path_work,"/blast_manual/manual_select_asvs_18S.xlsx",sep=""),sheet=1)

table3=merge(table2,query_18S,by.x="asv_id",by.y="Feature.ID")
write.xlsx(table3, paste(path_work,"/",dir_save[2],"/",names_tablets[2],"_species_geoinfo_abundancesample.xlsx",sep="")) 


gbif_results_geoinfo <- table %>% 
  mutate(occurrence=rep(1,nrow(table))) %>% 
  distinct() %>%
  group_by(species,countryCode,country,decimalLatitude,decimalLongitude,continent,asv_id,asv_country,mean_pident) %>% 
  dplyr::summarise(sum_occurrence = sum(occurrence)) %>% 
  dplyr::select(species,asv_country,asv_id,decimalLatitude,decimalLongitude,mean_pident)


# representation manual Gbif and globalFungi#####


for (i in (1:length(countries_separation))) {
  print(countries_separation[i])
  
  gbif_results_geoinfo_subset <- gbif_results_geoinfo %>% 
    filter(asv_country==countries_separation[i])
  
  gbif_results_geoinfo_subset2 <- gbif_results_geoinfo_subset %>%
    group_by(species,decimalLatitude,decimalLongitude) %>%
    na.omit() %>%
    dplyr::summarise(mean_pident2=mean(round(mean_pident),2), n_asv=n_distinct(asv_id)) %>%
    mutate(species2=paste(species, " (",n_asv," ASVs, ",mean_pident2,"%)",sep = "")) %>%
    filter(species!="Zea mays")
  
  levels_species <- na.omit(unique(gbif_results_geoinfo_subset2$species2[order(c(-gbif_results_geoinfo_subset2$n_asv))]))
  viridis_colors <- viridis(length(levels_species), option = "C")
  
  ### using GBIF ------    
  
  map_data <- 
    base_world +
    geom_point(data=gbif_results_geoinfo_subset2, 
               aes(x=decimalLongitude, y=decimalLatitude, fill=factor(species2, levels = levels_species)), colour="black", 
               pch=21, size=5, alpha=I(0.7)) +
    scale_fill_manual(values = setNames(viridis_colors, levels_species), name="Possible species identity of ASVs") +
    theme_rep +
    guides(fill = guide_legend(ncol = 2, title.position="top"))
            
            
    ggsave(paste(path_work,"/blast_manual/manual",countries_separation[i],'_map.pdf',sep=""),map_data,device = cairo_pdf,width = 30, height = 15, units = "cm")
            
            
  subtitle_rep=paste("Distribution of",length(unique(gbif_results_geoinfo_subset2$species)),"possible species, represented by",length(unique(gbif_results_geoinfo_subset$asv_id)), "da-ASVs")
            
  map_data2 <- 
        base_world + 
        stat_bin_hex(data = gbif_results_geoinfo_subset2, aes(x = decimalLongitude, y = decimalLatitude), bins = 10, color = "black") +
        stat_bin_hex(data = gbif_results_geoinfo_subset2, aes(x = decimalLongitude, y = decimalLatitude, label=after_stat(count)), bins = 10, color = "black", geom="text") +
        scale_fill_gradient(low="white",high="yellow") +
        theme(legend.position="right")
            
  ggsave(paste(path_work,"/blast_manual/manual",countries_separation[i],'_map2.pdf',sep=""),map_data2,device = cairo_pdf,
                width = 15, height = 10, units = "cm")
            ##### IMPORTANT -> da-ASV remaining of THAILAND match to Pichia sp. More than 100 species within this genus. Distribution of Pichia can be retrieved from GBIF
            
            
            ### using globalFungi ------    
            
  fung_select=unique(globalfung_results$species[(globalfung_results$species %in% gbif_results_geoinfo_subset2$species)])   
            
  globalfung_results2 <- globalfung_results[globalfung_results$species %in% fung_select,]
  globalfung_results2$species2=gbif_results_geoinfo_subset2$species2[match(globalfung_results2$species,gbif_results_geoinfo_subset2$species)]

  globalfung_map <- 
      base_world +
      geom_point(data=globalfung_results2, aes(x = longitude, y = latitude, fill=factor(species2, levels = levels_species)), colour="black",pch=21,size=5, alpha=I(0.7)) +  
      scale_fill_manual(values = setNames(viridis_colors, levels_species), name="Possible species identity of ASVs") +
      theme_rep +
      guides(fill = guide_legend(ncol = 2, title.position="top"))
            
             
  subtitle_rep=paste("Distribution of",length(unique(globalfung_results2$species)),"possible species, represented by",length(unique(gbif_results_geoinfo$asv_id[gbif_results_geoinfo$species %in% globalfung_results2$species])), "da-ASVs")
  
  globalfung_results3 <- globalfung_results2[rep(1:nrow(globalfung_results2), globalfung_results2$abundances),]
            
  globalfung_map2 <- 
    base_world + 
    stat_bin_hex(data = globalfung_results3, aes(x = longitude, y = latitude), bins = 10, color = "black") +
    stat_bin_hex(data = globalfung_results3, aes(x = longitude, y = latitude, label=after_stat(count)), bins = 10, color = "black", geom="text") +
    scale_fill_gradient(low="white",high="yellow") +
    theme(legend.position="right") 
            
 
  assign(paste(countries_separation[i],"manual_gbif_map",sep ="_"),map_data)
  assign(paste(countries_separation[i],"manual_gbif_map2",sep ="_"),map_data2)
  assign(paste(countries_separation[i],"manual_globalfungi_map",sep ="_"),globalfung_map)
  assign(paste(countries_separation[i],"manual_globalfungi_map2",sep ="_"),globalfung_map2)
  rm(globalfung_results2,globalfung_results3)
  
}


grattan_save_pptx(E_manual_gbif_map, filename=paste(path_work,"/pp_maps/manual_",countries_separation[i],'_map.pptx',sep=""))
grattan_save_pptx(E_manual_gbif_map2, filename=paste(path_work,"/pp_maps/manual_",countries_separation[i],'_map2.pptx',sep=""))
grattan_save_pptx(E_manual_globalfungi_map, filename=paste(path_work,"/pp_maps/manual_",countries_separation[i],'_globalfungi_map.pptx',sep=""))
grattan_save_pptx(E_manual_globalfungi_map2, filename=paste(path_work,"/pp_maps/manual_",countries_separation[i],'_globalfungi_map2.pptx',sep=""))



####### ALTERNATIVELY ANNOTATE TO EXACT TAXONOMY, SPECIFIC GENUSES OR SPECIES AND REPRESENT IN UNIFIED MAP #######

taxonomy_18S=read_qza("taxonomy_merged2_18S.qza")
taxonomy_18S=taxonomy_18S$data

taxonomy_18S[c("domain", "supergroup", "division", "subdivision", "class", "order","family","genus","species")] <-str_split_fixed(taxonomy_18S$Taxon,";",9) 
taxonomy_18S=taxonomy_18S[taxonomy_18S$domain=="Eukaryota",]

for (h in 1:length(tables_query)) {
  print(tables_query[h]) 
  
  query_18S=read.xlsx(paste(path_work,tables_query[h],".xlsx",sep=""),sheet=1)
  query_18S=query_18S[!query_18S$Feature.ID %in% excipients_count_rar$X,]
  
  if (names_tablets[h]=="da_all"){
    query_18S$asv_country=ifelse(query_18S$log2FoldChange > 0, "E","T") }
  
  query_18S_tax=merge(query_18S,taxonomy_18S, by="Feature.ID")
  query_18S_tax=query_18S_tax[query_18S_tax$species!="",]
  
  assign(paste("exact_species_18s",names_tablets[h],sep ="_"),query_18S_tax)
}  



tables_gbif2=c("exact_species_18s_da_all","exact_species_18s_manual")

for (h in 1:length(tables_gbif)){
  get(tables_gbif2[h]) -> table
  
  gbif_results_geoinfo <- table %>% 
    mutate(occurrence=rep(1,nrow(table))) %>% 
    distinct() %>%
    group_by(species,countryCode,country,decimalLatitude,decimalLongitude,continent,asv_id,asv_country,log2foldchange,mean_pident) %>% 
    dplyr::summarise(sum_occurrence = sum(occurrence))
  
  for (i in (1:length(countries_separation))) {
    print(countries_separation[i])
    
    gbif_results_geoinfo_subset <- gbif_results_geoinfo %>% 
      filter(asv_country==countries_separation[i])
    
    gbif_results_geoinfo_subset2 <- gbif_results_geoinfo_subset %>%
      group_by(species,decimalLatitude,decimalLongitude) %>%
      na.omit() %>%
      dplyr::summarise(mean_log2foldchange=mean(round(log2foldchange,2)),mean_pident2=mean(round(mean_pident),2), n_asv=n_distinct(asv_id)) %>%
      mutate(species2=paste(species, " (",n_asv," ASVs, ",mean_pident2,"%, ",mean_log2foldchange,")",sep = "")) %>%
      filter(species!="Zea mays")
    
    levels_species <- unique(gbif_results_geoinfo_subset2$species2[order(-gbif_results_geoinfo_subset2$n_asv)])
    
    #colors for representation
    set.seed(723451)
    pal = createPalette(length(levels_species),  c("#ff0000", "#00ff00", "#0000ff"))
    print(pal)
    
    names(pal) = levels_species
    
    #color_vector=c("brown","blue","gold","chartreuse3","darkorchid","deepskyblue","peachpuff","seagreen")
    
    map_data <- 
      base_world +
      geom_point(data=gbif_results_geoinfo_subset2, 
                 aes(x=decimalLongitude, y=decimalLatitude, fill=factor(species2, levels = levels_species)), colour="black", 
                 pch=21, size=3, alpha=I(0.7)) +  scale_fill_manual(values=pal, name="Possible species identity of ASVs") +
      ggtitle(names_plot[i]) #+
    #theme(legend.position = c(0.1, 0.2))
    
    print(map_data)
    
    ggsave(paste(path_work,"/gbif_figures_qiimeblast/map_occurrences_",names_tablets[h],countries_separation[i],'.pdf',sep=""),map_data,device = cairo_pdf,
           width = 30, height = 15, units = "cm")
    
    }
}




