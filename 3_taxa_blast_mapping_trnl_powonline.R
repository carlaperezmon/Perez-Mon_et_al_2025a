###### BLAST trnl POW differential abundant, excipients removed ######

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


packages <- c("R.utils","Biostrings","phyloseq","Biostrings","ggplot2","openxlsx","ape","reshape2","qiime2R",
              "stringr","dplyr","Polychrome","Rmisc","rgl","EcolUtils","vegan","plot3D","gridExtra","scatterplot3d",
              "indicspecies", "DESeq2","gplots","Rmisc","rworldmap","countrycode",'rWCVP','BIEN',"kewr","viridis",
              "ggpubr","gridExtra","rWCVPdata","maps","rnaturalearth")

ipak(packages)


### load tables ##### 

all_blast=read.table("blast_all/blast_trnl_pe0.90_qc0.90_NCBItrnl_all_export/blast6.tsv",sep = "\t")
colnames(all_blast)=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

taxonomy_identities=read.csv(paste(first_part,"Conservation_Genetics/GENOME_PROJECTS/Carla/sequencing/databases/trnl_ncbi/taxonomy_r/taxa_tableforqiime.txt",sep=""),sep='\t', header=FALSE)
taxonomy_identities[c("domain", "phylum", "class", "order", "family", "genus","species")] <- str_split_fixed(taxonomy_identities$V2,";",7) 

all_blast$species=taxonomy_identities$species[match(all_blast$sseqid,taxonomy_identities$V1)]
all_blast$species=gsub("s__","", all_blast$species)

all_blast_ev=all_blast %>%
  dplyr::count(qseqid, sort = TRUE)


asvs_names=read.csv2("asvs_rep_names_trnl.csv")

excipients_count_rar=read.csv("excipients_counts_rar_trnl.csv", header=TRUE,sep = ";",check.names=FALSE)

dir_save=c("blast_deseq","blast_manual") #since not ancombc results after removing excipients sequences, only use the deseq results
tables_query=c("/tables_deseq/asv_level/deseq_table_noextrnl_all","/blast_manual/manual_select_asvs_trnl")
names_tablets=c("da_all","manual")
asv_subtitle=c("da-ASVs","ASVs")


#after removing ASVs containing Zea (with at least 0.9% identity, so pretty sure that it is not zea mays)

for (h in 1:length(tables_query)) {
  print(tables_query[h])
  
  #load data and filter tables
  query_trnl=read.xlsx(paste(path_work,tables_query[h],".xlsx",sep=""),sheet=1)
  query_trnl=query_trnl[!query_trnl$Feature.ID %in% excipients_count_rar$X,]
  
  if (names_tablets[h]=="da_all"){
  query_trnl$asv_country=ifelse(query_trnl$log2FoldChange > 0, "E","T") }
  
  query_trnl$country_seq=asvs_names$rep_name[match(query_trnl$Feature.ID,asvs_names$Feature.ID)]

  trnl_blast_table=all_blast[all_blast$qseqid %in% query_trnl$Feature.ID,]
  
  asv_ini=length(unique(query_trnl$Feature.ID))
  
  all(unique(trnl_blast_table$qseqid) %in% unique(query_trnl$Feature.ID))
  all(unique(query_trnl$Feature.ID) %in% unique(trnl_blast_table$qseqid))
  
  ####remove zea
  zea_remove=unique(trnl_blast_table$qseqid[grep( "Zea",trnl_blast_table$species)]) 
  trnl_blast_table=trnl_blast_table[!(trnl_blast_table$qseqid %in% zea_remove),]
  
  query_trnl=query_trnl[query_trnl$Feature.ID %in% trnl_blast_table$qseqid,]
  
  asv_zearm=length(unique(query_trnl$Feature.ID))
  
  ### remove not hits
  
  trnl_blast_table=trnl_blast_table[trnl_blast_table$pident>0,]
  query_trnl=query_trnl[query_trnl$Feature.ID %in% trnl_blast_table$qseqid,]
  
  asv_nohit=length(query_trnl$Feature.ID)
  
  write.csv(cbind(asv_ini=asv_ini, asv_zearm=asv_zearm, asv_nohit=asv_nohit), paste(path_work,"/",dir_save[h],"/",names_tablets[h],"_remaining_asvs_after_zea_removed.csv",sep=""))

  feature_ids=query_trnl$Feature.ID
  
  best_hits_selection=list()
  
  for (i in 1:length(feature_ids)) {
    print(feature_ids[i])
    
    max_hits <- trnl_blast_table %>%
      filter(qseqid==feature_ids[i]) %>%
      filter(pident==max(pident))
    
    best_hits_selection [[i]] <- max_hits
    
  }
  
  all_best_hits_selection=dplyr::bind_rows(best_hits_selection)
  
  all_best_hits_selection2 <- all_best_hits_selection  %>% 
    dplyr::select(qseqid,pident,species) %>%
    dplyr::group_by(qseqid,species) %>%
    dplyr::summarise(mean_pident=mean(pident))
  
  ### merge blast information with botanical information -> geographic area and distribution =====
  
  trnl_blast_map=merge(query_trnl,all_best_hits_selection2,by.x="Feature.ID",by.y="qseqid")
  
  ### select a similarity threshold to choose from and scale the values according to the selected threshold ====
  
  threshold=99

  trnl_blast_map2=trnl_blast_map[trnl_blast_map$mean_pident>=threshold,]
  trnl_blast_map2$occurrence=rep(1,nrow(trnl_blast_map2))
  
  assign(paste("trnl_blast_map",names_tablets[h],sep ="_"),trnl_blast_map2)
  write.table(trnl_blast_map2, paste(path_work,"/",dir_save[h],"/",names_tablets[h],"_trnl_blast_map.txt",sep=""), append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
}  
  
  

#### Build distribution tables with POW using tools online and store dist maps of species separately ------------  

names <- rWCVPdata::wcvp_names
distributions <- rWCVPdata::wcvp_distributions

countries_separation=c("E","T")
names_plot=c("England", "Thailand")

world_map <- ne_countries(scale = "medium", returnclass = "sf")

base_world <- ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey80")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "white"), legend.position="left",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())

asv_subtitle=c("da-ASVS","ASV")

for (h in 1:length(names_tablets)) {
  print(names_tablets[h])
  
  for (i in (1:length(countries_separation))) {
    print(countries_separation[i])
    
    get(paste("trnl_blast_map_",names_tablets[h],sep ="")) -> trnl_blast_map2  
    
    trnl_blast_map2 <-trnl_blast_map2 %>% 
      filter(asv_country==countries_separation[i])
    
    asv_id=unique(trnl_blast_map2$Feature.ID)
    
    all_distributions=data.frame()

  for (j in (1:length(asv_id))) {
      print(asv_id[j])
    
      asv_id2=asvs_names$rep_name[match(asv_id[j],asvs_names$Feature.ID)]
      
      species_to_query <- unique(trnl_blast_map2$species[trnl_blast_map2$Feature.ID==asv_id[j]])
      
      plot_analysis <- list()
      
      for (k in (1:length(species_to_query))) {
        print(species_to_query[k])

        powo_results=purrr::map_dfr(search_powo(species_to_query[k])$results, unlist)
        
        if(nrow(powo_results)==0) {
        print(species_to_query[k])
        
        print ("not found")} else {
      
        species_to_query2=powo_results$name[powo_results$accepted==TRUE]
        
        try({
          
        distribution <- wcvp_distribution(species_to_query2[1], taxon_rank = "species", wcvp_names = names, 
                                          wcvp_distributions = distributions, extinct = FALSE,
                                          location_doubtful = FALSE)
          
        plot_analysis[[paste0(species_to_query[k],'_plot')]] <-wcvp_distribution_map(distribution) + ggtitle(paste(species_to_query[k], "=",species_to_query2))
        
        distribution$species=rep(species_to_query[k], nrow(distribution))
        distribution$asv_id=rep(asv_id[j], nrow(distribution))
        distribution$country_seq=asv_id2
        
        all_distributions=rbind(all_distributions,distribution)
        
        },silent = TRUE) }
      
      }

      if (length(plot_analysis) ==0) { print ("nothing to plot")} else {
      if (length(plot_analysis) > 4) {
      range_loop=seq(1,length(plot_analysis),by=4) 
      
      pdf(file = paste(path_work,"/",dir_save[h],"/pow_ind_asvs/",countries_separation[i],"/",asv_id2,"_noZea_species.pdf",sep=""), paper="a4r",width=10, height=7)
       for(o in 1:length(range_loop)) {
        l=range_loop[o]+3
        
        if (l < length(plot_analysis)) {
          plot_analysis2=plot_analysis[range_loop[o]:l]
          
        } else { 
          plot_analysis2=plot_analysis[range_loop[o]:length(plot_analysis)]
        }
        do.call('grid.arrange',c(plot_analysis2, ncol= 2,nrow=2))
      }
      
      dev.off()      
  
      
      
      } else { 
      pdf(file = paste(path_work,"/",dir_save[h],"/pow_ind_asvs/",countries_separation[i],"/",asv_id2,"_noZea_species.pdf",sep=""), paper="a4r",width=10, height=7)
      do.call('grid.arrange',c(plot_analysis, ncol= 2,nrow=2))  
      dev.off()
      }      
      }
  }

    assign(paste("all_distributions_asvs",names_tablets[h],countries_separation[i],sep ="_"),all_distributions)
    
    }
      
}
    
    


### REPRESENTATION POW for asv separately (aggregation of all possible sps per ASV) ####  

for (h in 1:length(names_tablets)) {
  print(names_tablets[h]) 
  
  for (i in (1:length(countries_separation))) {
  print(countries_separation[i]) 
  
  get(paste("all_distributions_asvs",names_tablets[h],countries_separation[i],sep ="_")) -> all_distributions
  
  all_distributions=all_distributions[!grepl("Adansonia", all_distributions$species),]
  
  asv_id=unique(all_distributions$asv_id)
  
  for (j in (1:length(asv_id))) {
    print(asv_id[j])
    
    asv_id2=asvs_names$rep_name[match(asv_id[j],asvs_names$Feature.ID)]
    
    possible_sps=paste(unique(all_distributions$species[all_distributions$asv_id==asv_id[j]]), collapse = ",")
    
    all_distributions_subset=all_distributions[all_distributions$asv_id==asv_id[j],]
    all_distributions_subset$occurrence=rep(1,nrow(all_distributions_subset))
    all_distributions_subset$rel_occurrence=all_distributions_subset$occurrence/nrow(all_distributions_subset)
    
    all_distributions_subset <- all_distributions_subset %>%
      group_by(LEVEL3_COD) %>%
      dplyr::summarise(sum_occurrence = sum(occurrence), sum_rel_occurrence=sum(rel_occurrence)) %>%
      arrange(sum_occurrence)
    
    
    color_palette <- viridis(10)
    
    # Map numeric values to colors
    all_distributions_subset$numeric_colors <- color_palette[cut(all_distributions_subset$sum_occurrence, 10)]
    
    bot_regions=unique(all_distributions_subset$LEVEL3_COD)
    
    p <- base_world
    
    
    for (l in (1:length(bot_regions))) {
      print(unique(bot_regions[l]))
      
      all_distributions_subset2=all_distributions_subset[all_distributions_subset$LEVEL3_COD==bot_regions[l],]
      
      p=p + geom_sf(data = all_distributions_subset2$geometry[1], fill=all_distributions_subset2$numeric_colors)
      
    }
    
    p2=p+ggtitle(paste("asv_id:",asv_id2)) + labs(subtitle=str_wrap(paste('Aggregated distribution of species:', possible_sps), 100)) + theme(plot.subtitle=element_text(face = "italic"))
    
    color_bar=cbind(all_distributions_subset$sum_occurrence,all_distributions_subset$numeric_colors)
    color_bar=data.frame(unique(color_bar))
    color_bar$y=rep(1,nrow(color_bar))
    levels=as.character(color_bar$X1)
    
    color_bar_plot=ggplot(color_bar, aes(x=factor(X1, levels=levels), y=y)) + 
      geom_bar(stat = "identity",  fill=color_bar$X2) + 
      scale_y_continuous(expand = c(0,0)) + 
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.x=element_blank())
    
    
    all_plot=ggpubr::ggarrange(p2, color_bar_plot, heights = c(10, 0.7),widths = c(10, 0.7),
                       ncol = 1, nrow = 2, align = "none")
    
    
    ggsave(paste(path_work,"/",dir_save[h],"/pow_ind_asvs/",countries_separation[i],"/",names_tablets[h],"_",asv_id2,"_noZea_sps_aggregated.pdf",sep=""),all_plot,device = cairo_pdf,
           width = 30, height = 15, units = "cm")
}
  
}
}



### Representation of all species distributions for all ASVs aggregated ####  

for (h in (1:length(names_tablets))) {
  print(names_tablets[h])
  
  for (i in (1:length(countries_separation))) {
  print(countries_separation[i])
  
  get(paste("all_distributions_asvs",names_tablets[h],countries_separation[i],sep ="_")) -> all_distributions_subset

  all_distributions_subset=all_distributions_subset[!grepl("Adansonia", all_distributions_subset$species),]
  
  all_distributions_subset$occurrence=rep(1,nrow(all_distributions_subset))
  all_distributions_subset$rel_occurrence=all_distributions_subset$occurrence/nrow(all_distributions_subset)
  
  possible_sps=unique(gsub("^(\\w)\\w*\\s", "\\1. ", all_distributions_subset$species))

  genera_table <- all_distributions_subset %>%
    mutate(genera=word(species,1))  %>%
    group_by(genera) %>%
    dplyr::summarise(sum_species=n_distinct(species)) %>%
    mutate(gen_sp=paste(genera," (",sum_species, " sps)",sep=""))
  
  genera=paste(genera_table$gen_sp, collapse = ", ")
  subtitle_rep=paste("Distribution of",length(unique(all_distributions_subset$species)),"possible species, represented by",length(unique(all_distributions_subset$asv_id)), asv_subtitle[h])
  
  subtitle_rep2=paste(subtitle_rep,genera,sep="\n")

  write.csv(genera,paste(path_work,"/pp_maps/",countries_separation[i],"_",names_tablets[h],"_forrep.csv",sep=""))
  
  all_distributions_subset$genera <- word(all_distributions_subset$species,1)
  
  all_distributions_subset <- all_distributions_subset %>%
    group_by(LEVEL3_COD) %>%
    dplyr::summarise(sum_occurrence=sum(occurrence), sum_rel_occurrence=sum(rel_occurrence)) %>%
    arrange(sum_occurrence)
    
  
  color_palette <- viridis(10)
  
  # Map numeric values to colors
  all_distributions_subset$numeric_colors <- color_palette[cut(all_distributions_subset$sum_occurrence, 10)]
  
  bot_regions=unique(all_distributions_subset$LEVEL3_COD)
  
  p <- base_world
  
 
  for (l in (1:length(bot_regions))) {
    print(unique(bot_regions[l]))
    
    all_distributions_subset2=all_distributions_subset[all_distributions_subset$LEVEL3_COD==bot_regions[l],]
    
    p=p + geom_sf(data = all_distributions_subset2$geometry[1], fill=all_distributions_subset2$numeric_colors, color="grey20")
    
  }
  
  
  p2=p+ggtitle(names_plot[i]) + 
    labs(caption=cat(subtitle_rep2)) + 
    theme(plot.caption = element_text(hjust = 0))

  color_bar=data.frame(x=seq(min(all_distributions_subset$sum_occurrence),max(all_distributions_subset$sum_occurrence),round(max(all_distributions_subset$sum_occurrence)/10,0)),
                       y=0.5,
                       color=seq(min(all_distributions_subset$sum_occurrence),max(all_distributions_subset$sum_occurrence),round(max(all_distributions_subset$sum_occurrence)/10,0)))
  
  color_bar_plot=ggplot(color_bar, aes(x=x, y = y, fill=color)) +
    geom_tile() +
    coord_flip()+
    scale_x_continuous(expand = c(0,0), position="top")+
    scale_y_continuous(breaks = c(0,0.5))+
    scale_fill_viridis(option = "viridis") +  # Use the viridis color palette
    theme_void() +  # Remove axis, grid, and background
    theme(legend.position = "none",
          axis.text.y=element_text(color="black", size=14),plot.margin = unit(c(5, 9.5, 5, 10), "cm"))
  

  grattan_save_pptx(p2, filename=paste(path_work,"/pp_maps/",countries_separation[i],"_",names_tablets[h],"_map.pptx",sep=""))
  grattan_save_pptx(color_bar_plot, filename=paste(path_work,"/pp_maps/",countries_separation[i],"_",names_tablets[h],"_legend.pptx",sep=""))
  ggsave(paste(path_work,"/",dir_save[h],"/pow_ind_asvs/",countries_separation[i],"/",names_tablets[h],"allasvs_noZea_sps_aggregated.pdf",sep=""),p2,device = cairo_pdf,
         width = 30, height = 15, units = "cm")
  
}  
} 
  

### Representation for excipients #### 

metadata=read.xlsx(paste(first_part,"Conservation_Genetics/GENOME_PROJECTS/Carla/sequencing/sequencing_results/data_analyses/2024-03_merged_with_june2023/metadata_bvsh_allset.xlsx",sep="/"), sheet=1)
metadata_exc=metadata[metadata$material=="excipient",]

excipients_count_rar2=melt(excipients_count_rar)
excipients_count_rar2$excip=metadata_exc$lot_id2[match(excipients_count_rar2$variable, metadata_exc$`sample-id`)]
colnames(excipients_count_rar2)[1]<-"Feature.ID"

excipients_count_rar2=excipients_count_rar2 %>%
  group_by(Feature.ID,excip) %>%
  dplyr::summarise(Total = sum(value), .groups = 'drop')

names_excipients=unique(excipients_count_rar2$excip)


for (h in 1:length(names_excipients)) {
  print(names_excipients[h])
  
  excipients_subset=excipients_count_rar2[as.character(excipients_count_rar2$excip)==names_excipients[h],] 
  excipients_subset=excipients_subset[excipients_subset$Total > 9,]
    
  excipients_blast_table=all_blast[all_blast$qseqid %in% excipients_subset$Feature.ID,]

  asv_ini=length(unique(excipients_subset$Feature.ID))
  
  all(unique(excipients_blast_table$qseqid) %in% unique(excipients_subset$Feature.ID))
  all(unique(excipients_subset$Feature.ID) %in% unique(excipients_blast_table$qseqid))
  
  ####remove zea
  zea_remove=unique(excipients_blast_table$qseqid[grep( "Zea",excipients_blast_table$species)]) 
  excipients_blast_table=excipients_blast_table[!(excipients_blast_table$qseqid %in% zea_remove),]
  
  excipients_subset=excipients_subset[excipients_subset$Feature.ID %in% excipients_blast_table$qseqid,]
  
  asv_zearm=length(unique(excipients_subset$Feature.ID))
  
  ### remove not hits
  
  excipients_blast_table=excipients_blast_table[excipients_blast_table$pident>0,]
  excipients_subset=excipients_subset[excipients_subset$Feature.ID %in% excipients_blast_table$qseqid,]
  
  asv_nohit=length(excipients_subset$Feature.ID)
  
  write.csv(cbind(asv_ini=asv_ini, asv_zearm=asv_zearm, asv_nohit=asv_nohit), paste(path_work,"/blast_excipients/",names_excipients[h],"_remaining_asvs_after_zea_removed.csv",sep=""))
  
  feature_ids=excipients_subset$Feature.ID
  
  best_hits_selection=list()
  
  for (i in 1:length(feature_ids)) {
    print(feature_ids[i])
    
    max_hits <- excipients_blast_table %>%
      filter(qseqid==feature_ids[i]) %>%
      filter(pident==max(pident))
    
    best_hits_selection [[i]] <- max_hits
    
  }
  
  all_best_hits_selection=dplyr::bind_rows(best_hits_selection)
  
  all_best_hits_selection2 <- all_best_hits_selection  %>% 
    dplyr::select(qseqid,pident,species) %>%
    dplyr::group_by(qseqid,species) %>%
    dplyr::summarise(mean_pident=mean(pident))
  
  ### merge blast information with botanical information -> geographic area and distribution =====
  
  trnl_blast_map=merge(excipients_subset,all_best_hits_selection2,by.x="Feature.ID",by.y="qseqid")
  
  
  ### select a similarity threshold to choose from and scale the values according to the selected threshold ====
  
  threshold=99
  

  trnl_blast_map2=trnl_blast_map[trnl_blast_map$mean_pident>=threshold,]
  trnl_blast_map2$occurrence=rep(1,nrow(trnl_blast_map2))
  
  assign(paste("trnl_blast_map",names_excipients[h],sep ="_"),trnl_blast_map2)
  write.table(trnl_blast_map2, paste(path_work,"/blast_excipients/",names_excipients[h],"_trnl_blast_map.txt",sep=""), append = FALSE, sep = " ", dec = ".",
              row.names = TRUE, col.names = TRUE)
}  


### POWO collection of samples for excipients to create aggregation maps

for (h in 1:length(names_excipients)) {
  print(names_excipients[h])
  
  get(paste("trnl_blast_map_",names_excipients[h],sep ="")) -> trnl_blast_map2  
    
  asv_id=unique(trnl_blast_map2$Feature.ID)
    
  all_distributions=data.frame()
    
  for (j in (1:length(asv_id))) {
    print(asv_id[j])
      
    asv_id2=asvs_names$rep_name[match(asv_id[j],asvs_names$Feature.ID)]
      
    species_to_query <- unique(trnl_blast_map2$species[trnl_blast_map2$Feature.ID==asv_id[j]])
      
    for (k in (1:length(species_to_query))) {
        print(species_to_query[k])
        
        powo_results=purrr::map_dfr(search_powo(species_to_query[k])$results, unlist)
        
        if(nrow(powo_results)==0) {
          print(species_to_query[k])
          
          print ("not found")} else {
            
            species_to_query2=powo_results$name[powo_results$accepted==TRUE]
            
            try({
              
              distribution <- wcvp_distribution(species_to_query2[1], taxon_rank = "species", wcvp_names = names, 
                                                wcvp_distributions = distributions, extinct = FALSE,
                                                location_doubtful = FALSE)
              
              distribution$species=rep(species_to_query[k], nrow(distribution))
              distribution$asv_id=rep(asv_id[j], nrow(distribution))
              distribution$country_seq=asv_id2
              
              all_distributions=rbind(all_distributions,distribution)
              
            },silent = TRUE) }
        
      }
      
  }
  
  assign(paste("all_distributions_asvs",names_excipients[h],sep ="_"),all_distributions)
    
}


### Representation of all species distributions for all ASVs aggregated ####  

for (h in 1:length(names_excipients)) {
  print(names_excipients[h])
  
  get(paste("all_distributions_asvs",names_excipients[h],sep ="_")) -> all_distributions_subset
    
    all_distributions_subset=all_distributions_subset[!grepl("Adansonia", all_distributions_subset$species),]
    
    subtitle_rep=paste(length(unique(all_distributions_subset$asv_id)),"ASVs, representing a total number of",length(unique(all_distributions_subset$species)), "possible species")
    
    all_distributions_subset$occurrence=rep(1,nrow(all_distributions_subset))
    all_distributions_subset$rel_occurrence=all_distributions_subset$occurrence/nrow(all_distributions_subset)
    
    possible_sps=paste(unique(all_distributions_subset$species), collapse = ",")
    
    all_distributions_subset <- all_distributions_subset %>%
      group_by(LEVEL3_COD) %>%
      dplyr::summarise(sum_occurrence = sum(occurrence), sum_rel_occurrence=sum(rel_occurrence)) %>%
      arrange(sum_occurrence)
    
    
    color_palette <- viridis(10)
    
    # Map numeric values to colors
    all_distributions_subset$numeric_colors <- color_palette[cut(all_distributions_subset$sum_occurrence, 10)]
    
    bot_regions=unique(all_distributions_subset$LEVEL3_COD)
    
    
    p <- base_world
    
    
    for (l in (1:length(bot_regions))) {
      print(unique(bot_regions[l]))
      
      all_distributions_subset2=all_distributions_subset[all_distributions_subset$LEVEL3_COD==bot_regions[l],]
      
      p=p + geom_sf(data = all_distributions_subset2$geometry[1], fill=all_distributions_subset2$numeric_colors)
      
    }
    
    
    p2=p+ggtitle(paste(names_excipients[h],"all asvs aggregated")) + labs(subtitle=subtitle_rep) + theme(plot.subtitle=element_text(face = "italic"))
    
    color_bar=cbind(all_distributions_subset$sum_occurrence,all_distributions_subset$numeric_colors)
    color_bar=data.frame(color_bar[!duplicated(color_bar),])
    color_bar$y=rep(1,nrow(color_bar))
    levels=unique(color_bar$X1)
    
    color_bar_plot=ggplot(color_bar, aes(x=factor(X1, levels=levels), y=y)) + 
      geom_bar(stat = "identity",  fill=color_bar$X2) + 
      scale_y_continuous(expand = c(0,0)) + 
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.x=element_blank())
    
    
    all_plot=ggarrange(p2, color_bar_plot, heights = c(10, 0.7),widths = c(10, 0.7),
                       ncol = 1, nrow = 2, align = "none")
    
    all_plot
    
    grattan_save_pptx(p2, filename=paste(path_work,"/pp_maps/",names_excipients[h],"_map.pptx",sep=""))
    grattan_save_pptx(color_bar_plot, filename=paste(path_work,"/pp_maps/",names_excipients[h],"_legend.pptx",sep=""))
    ggsave(paste(path_work,"/blast_excipients/",names_excipients[h],"allasvs_noZea_sps_aggregated.pdf",sep=""),all_plot,device = cairo_pdf,
           width = 30, height = 15, units = "cm")
    
  }  
 



###### tables all ASVs and distribution for tablets and excipients #####

for (h in (1:length(names_tablets))) {
  print(names_tablets[h])
  get(paste("trnl_blast_map",names_tablets[h], sep="_")) -> blast_map
  
  for (i in (1:length(countries_separation))) {
    print(countries_separation[i])
    
    get(paste("all_distributions_asvs",names_tablets[h],countries_separation[i],sep ="_")) -> all_distributions
    
    all_distributions=all_distributions[!grepl("Adansonia", all_distributions$species),]
    
    species_allregions <- all_distributions %>%
      as.data.frame() %>%
      select(asv_id,species,LEVEL3_NAM)  %>%
      unique() %>%
      group_by(asv_id,species) %>%
      dplyr::summarise(botanic_regions = paste(LEVEL3_NAM, collapse = ", "), .groups = "drop")
    
    species_allregions$p_ident=blast_map$mean_pident[match(species_allregions$asv_id,blast_map$Feature.ID)]
    
  write.xlsx(species_allregions,paste(path_work,"/",dir_save[h],"/",countries_separation[i],"_",names_tablets[h],"_species_allregions.xlsx",sep=""))  
}      
}

for (h in (1:length(names_excipients))) {
  print(names_excipients[h])

  get(paste("trnl_blast_map",names_excipients[h], sep="_")) -> blast_map
  get(paste("all_distributions_asvs",names_excipients[h],sep ="_")) -> all_distributions_subset
  
  all_distributions_subset=all_distributions_subset[!grepl("Adansonia", all_distributions_subset$species),]
  
  species_allregions <- all_distributions_subset %>%
      as.data.frame() %>%
      select(asv_id,species,LEVEL3_NAM)  %>%
      unique() %>%
      group_by(asv_id,species) %>%
      dplyr::summarise(botanic_regions = paste(LEVEL3_NAM, collapse = ", "), .groups = "drop")
  
  species_allregions$p_ident=blast_map$mean_pident[match(species_allregions$asv_id,blast_map$Feature.ID)]
  
  write.xlsx(species_allregions,paste(path_work,"/blast_excipients/",names_excipients[h],"_species_allregions.xlsx",sep=""))  
}      




####### ALTERNATIVELY ANNOTATE TO EXACT TAXONOMY #######

taxonomy_trnl=read_qza("taxonomy_merged_trnl.qza")
taxonomy_trnl=taxonomy_trnl$data

taxonomy_trnl[c("domain", "phylum", "class", "order", "family", "genus","species")] <-str_split_fixed(taxonomy_trnl$Taxon,";",7) 

for (h in 1:length(tables_query)) {
  print(tables_query[h]) 

  query_trnl=read.xlsx(paste(path_work,tables_query[h],".xlsx",sep=""),sheet=1)
  query_trnl=query_trnl[!query_trnl$Feature.ID %in% excipients_count_rar$X,]
  
  if (names_tablets[h]=="da_all"){
    query_trnl$asv_country=ifelse(query_trnl$log2FoldChange > 0, "E","T") }
  
  query_trnl$country_seq=asvs_names$rep_name[match(query_trnl$Feature.ID,asvs_names$Feature.ID)]  
  query_trnl_tax=merge(query_trnl,taxonomy_trnl, by="Feature.ID")

  assign(paste("trnl_exact_species",names_tablets[h],sep ="_"),query_trnl_tax)
}  



for (h in 1:length(names_excipients)) {
  print(names_excipients[h])
  
  excipients_subset=excipients_count_rar2[as.character(excipients_count_rar2$excip)==names_excipients[h],] 
  excipients_subset=excipients_subset[excipients_subset$Total > 9,]
  
  excipients_subset_tax=merge(excipients_subset,taxonomy_trnl, by="Feature.ID")
  
  assign(paste("trnl_exact_species",names_excipients[h],sep ="_"),excipients_subset_tax)
}





