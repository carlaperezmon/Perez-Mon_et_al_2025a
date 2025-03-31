#######################################################
##### DIVERSITY AND DA ANALYSES OF merged march 2024 and june 2023 SAMPLES ========
#######################################################

#######################
####### LOADING DATA ####### 
###########################

##### set working directory and paths for data and figures ------

directory="2024-03_merged_with_june2023/asvs_single"


objects=c("table_merged","rep-seqs_merged") 

#from office
first_part="K:/ogden_grp/Conservation_Science/Conservation_Genetics/GENOME_PROJECTS/Carla/sequencing/sequencing_results/data_analyses"

path_work=paste(first_part,directory, sep="/")
path_work

setwd(path_work)

load("../1_diversity_da_analyses_16S_18S_trnl_march2024june2023_hist.RData")

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
              "stringr","dplyr","Polychrome","Rmisc","rgl","EcolUtils","vegan","gridExtra","scatterplot3d",
              "DESeq2","gplots","Rmisc","ANCOMBC","WRS2","ggforce","ggpubr",'rstatix','klaR','caret','plot3D',"viridis","adegenet","DA")

ipak(packages)


#load metadata ----------
marker=c("16S","18S","trnl")
directories=marker

metadata=read.xlsx("../metadata_bvsh_allset.xlsx", sheet= "Sheet1")
metadata=metadata[metadata$material!="multiplex",]
metadata=metadata[metadata$lot_id2!="PCR_neg_control_barcode",] # 0 reads in this control, so really good
colnames(metadata)[colnames(metadata)=="sample-id"] <-"sample.id"
metadata=metadata[!(metadata$sample.id %in% c("27360OR0031L01","27360OR0032L01","27360OR0063L01","27360OR0064L01",'30610OR0085L01')),] #to eliminate low reads samples and outliers (e.g. S2A, S5A, A-bk23 trnl and A-hd23 trnl from jun 2023 run)


rownames(metadata)=metadata$sample.id

for (i in 1:length(marker)) {
  print(marker[i]) 
  metadata_subset=metadata[metadata$marker==marker[i],]
  metadata_subset=metadata_subset[rowSums(is.na(metadata_subset)) != ncol(metadata_subset), ]
  assign(paste("metadata_",marker[i],sep=""),metadata_subset)
  rm(metadata_subset)
}
  

#load asvs tables and repseqs (needed for later) -------


for (i in 1:length(marker)) {
  print(marker[i]) 
  
  for (j in 1:length(objects)) {
    
    data_qza=read_qza(paste(path_work,"/",directories[i],"/",objects[j],"_",marker[i],".qza", sep=""))
    data=data_qza$data
    
    assign(paste(objects[j],"_",marker[i],sep=""),data)
    rm(data_qza,data)
  }
}

## set levels for reordering tables

levels_mat=c("tablet","excipient","water","air","control")
levels_country=c("Thailand","England","Water Thailand, air England", "DC excipient", "WG excipient")
levels_condition=c("no water, closed windows","no water, open windows","water, open windows","wet granulation formulation, open windows","cellulose (WG excipient, USA)","cellulose (DC excipient, USA)", "starch (WG excipient, Thailand)","starch (DC excipient, USA)","env. air (dust)","tap water")



##################################
# REMOVING READS FROM CONTROLS AND CREATE RAREFIED COUNTS =======
###############################


#process negative controls ----

number_filtering=9

subDir <- "tables_general_stats"

for (i in 1:length(marker)){
  get(paste(objects[1],marker[i],sep = "_")) -> count
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  
  count_unfilt=count
  
  mainDir <- paste(path_work,marker[i],sep="/")
  dir.create(file.path(mainDir, subDir), showWarnings = TRUE)
  
  #remove swabs and water filter blanks from environmental samples
  
  swab_control_asv=count[,colnames(count) %in% metadata_subset$sample.id[metadata_subset$material=="air" & metadata_subset$excipient_features=="control"]]
  swab_control_asv=swab_control_asv[rowSums(swab_control_asv) > number_filtering,]
    
  count[rownames(count) %in% rownames(swab_control_asv), colnames(count) %in% metadata_subset$sample.id[metadata_subset$material=="air" & metadata_subset$excipient_features=="dust"]] <- 0
    
  filter_control_asv=count[,colnames(count) %in% metadata_subset$sample.id[metadata_subset$material=="water" & metadata_subset$excipient_compound=="control"]]
  filter_control_asv=filter_control_asv[filter_control_asv > number_filtering]
    
  count[rownames(count) %in% names(filter_control_asv), colnames(count) %in% metadata_subset$sample.id[metadata_subset$material=="water" & metadata_subset$excipient_compound=="water"]] <- 0
    
  count=count[,!colnames(count) %in% metadata_subset$sample.id[metadata_subset$material=="air" & metadata_subset$excipient_features=="control"]]
  count=count[,!colnames(count) %in% metadata_subset$sample.id[metadata_subset$material=="water" & metadata_subset$excipient_compound=="control"]]
  
  #check=count[,"27360OR0031L01"]
  #check=check[names(check) %in% rownames(swab_control_asv)]
  
  #explore the negative controls
  names_neg_control=na.omit(metadata_subset$sample.id[metadata_subset$tablet_id=="PCR_neg_control"])
  
  neg_control=count[,names_neg_control]
  neg_control=neg_control[rowSums(neg_control) > number_filtering,]
  
  counts_filt1=count[!(rownames(count) %in% rownames(neg_control)),]
  
  names_ext_control=na.omit(metadata_subset$sample.id[metadata_subset$tablet_id=="ext_control"])
  
  ext_control=count[,names_ext_control]
  ext_control=ext_control[rowSums(ext_control) > number_filtering,]
  
  counts_filt2=counts_filt1[!(rownames(counts_filt1) %in% rownames(ext_control)),]
  
  counts_filt2=counts_filt2[rowSums(counts_filt2) > 9,]
  
  write.csv(counts_filt2, paste(mainDir,"/",subDir,"/count_filt_",marker[i],".csv",sep="")) 
  
  ##### save the negative controls for later
  
  neg_ext_controls=count[,c(names_ext_control,names_neg_control)]
  neg_ext_controls=neg_ext_controls[rowSums(neg_ext_controls) > number_filtering,]
  
  colnames(neg_ext_controls)=metadata_subset$lot_id2[match(colnames(neg_ext_controls),metadata_subset$sample.id)]
  
  write.csv(neg_ext_controls, paste(mainDir,"/",subDir,"/neg_ext_control_",marker[i],".csv",sep="")) 
  
  
  #### create filt stats table
  
  controls_filt_stats=rbind(initial_counts=colSums(count), 
                            filt_pcr_control_counts=colSums(counts_filt1), 
                            percent_remaining_counts1=round(colSums(counts_filt1)/colSums(count)*100),
                            filt_ext_control_counts=colSums(counts_filt2),
                            percent_remaining_counts2=round(colSums(counts_filt2)/colSums(count)*100))
  
  
  counts_ap=count
  counts_ap[counts_ap > 0] <- 1 
  colSums(counts_ap)
  
  asv_filt1=counts_filt1
  asv_filt1[asv_filt1 > 0] <- 1 
  colSums(asv_filt1)
  
  
  asv_filt2=counts_filt2
  asv_filt2[asv_filt2 > 0] <- 1 
  colSums(asv_filt2)
  
  
  controls_filt_stats=rbind(controls_filt_stats, 
                            inital_ASV=colSums(counts_ap),
                            filt_pcr_control_ASV=colSums(asv_filt1), 
                            percent_remaining_ASV1=round(colSums(asv_filt1)/colSums(counts_ap)*100),
                            filt_ext_control_ASV=colSums(asv_filt2),
                            percent_remaining_ASV2=round(colSums(asv_filt2)/colSums(counts_ap)*100))
  
  colnames(controls_filt_stats)=metadata_subset$lot_id2[match(colnames(controls_filt_stats),metadata_subset$sample.id)]
  
  order_colnames=metadata_subset[order(factor(metadata_subset$material, levels=levels_mat),factor(metadata_subset$origin_country, levels=levels_country),factor(metadata_subset$condition,levels = levels_condition)),]$lot_id2
  
  controls_filt_stats=controls_filt_stats[,order(match(colnames(controls_filt_stats), order_colnames))]
  
  write.csv(controls_filt_stats, paste(mainDir,"/",subDir,"/controls_filt_stats_",marker[i],".csv",sep="")) #no much more differences between discarding all ASV reads from controls or the ones with more than 10 reads so lets be strict and remove all ASV reads from controls
  
  assign(paste("count_filt",marker[i],sep="_"),counts_filt2)
  assign(paste("neg_ext_controls",marker[i],sep="_"),neg_ext_controls)
  
  #check waters
  
  check=as.data.frame(count_unfilt[,colnames(count_unfilt) %in% c('30610OR0038L01','30610OR0039L01',"30610OR0040L01","30610OR0034L01",'30610OR0035L01')])
  check=check[rowSums(check) > 0,]
  colnames(check)=metadata_subset$lot_id2[match(colnames(check),metadata_subset$sample.id)]
  check_ag=sapply(split.default(check, names(check)), rowSums, na.rm = TRUE)
  
  
} 




# create rarefied counts for diversity analyses----  

abundance_threshold=5000
rarDir=paste(abundance_threshold,"_rarefaction", sep="")

for (i in 1:length(marker)){
  get(paste("count_filt",marker[i],sep="_")) -> count
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  
  mainDir <- paste(path_work,marker[i],sep="/")
  dir.create(file.path(mainDir, rarDir), showWarnings = TRUE)
  
  count=count[,colSums(count) > abundance_threshold]
  
  t_count=t(count)
  print(min(rowSums(t_count))) #5264 for 16S, 12016 for 18S and 9064 for trnL
  
 
  rar_count= rrarefy.perm(t_count, sample=min(rowSums(t_count)), n=100, round.out=T)
  print(rowSums(rar_count))

  assign(paste("rar_count",marker[i],sep="_"),rar_count)
  
}  
  
  
#####################################
### ALPHA DIVERSITY =====
####################################

alpha_indeces=c("richness","shannon")
countries=c("England","Thailand")

subDir4 <- paste(abundance_threshold,"_rarefaction/figures_alphadiv",sep="")
subDir5 <- paste(abundance_threshold,"_rarefaction/tables_alphadiv",sep="")


for (i in 1:length(marker)){
  get(paste("rar_count",marker[i],sep="_")) -> count
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  
  mainDir <- paste(path_work,marker[i],sep="/")
  dir.create(file.path(mainDir, subDir4), showWarnings = FALSE)
  dir.create(file.path(mainDir, subDir5), showWarnings = FALSE)
  
  richness=decostand(count, method = "pa")
  richness=cbind.data.frame(samples=rownames(count), richness=(rowSums(richness)))
  
  shannon=diversity(count, index="shannon")
  shannon=cbind.data.frame(samples=names(shannon), shannon=shannon)
  
  for (m in 1:length(alpha_indeces)) {
    
    div_table=get(alpha_indeces[m])  
    
    colnames(div_table)[colnames(div_table) == alpha_indeces[m]] <- "index"   
    
    div_table$lot_id=as.factor(metadata_subset$lot_id2[match(div_table$samples,metadata_subset$sample.id)])
    div_table$origin_country=factor(metadata_subset$origin_country[match(div_table$samples,metadata_subset$sample.id)], levels=levels_country)
    div_table$condition=factor(metadata_subset$condition[match(div_table$samples,metadata_subset$sample.id)], levels=levels_condition)
    div_table$material=factor(metadata_subset$material[match(div_table$samples,metadata_subset$sample.id)], levels=levels_mat)
    div_table$env=factor(ifelse(div_table$material %in% c("air"),"Air","Tablets, excipients & water"), levels=c("Tablets, excipients & water","Air"))
    
    div_table_summ=summarySE(div_table, measurevar="index", groupvars=c("material","origin_country","condition","lot_id", "env"))
    
    div_table_summ=div_table_summ[order(div_table_summ$material,div_table_summ$origin_country, div_table_summ$condition,div_table_summ$env),]
    
    order_lotid=unique(div_table_summ$lot_id) 
    
    x_labels=c("A-BK1"="T-1","A-BK2"="T-2","B-BK"="T-3","A-HD1"="E-1","A-HD2"="E-2","B-HD"="E-3","B-HD3"="E-4","B-HD2"="E-5",
               "cellulose (DC excipient)"="Cel-DC", "starch (DC excipient)"="Str-DC","cellulose (WG excipient)"="Cel-WG", "starch (WG excipient)"="Str-WG",
               "water_bangkok"="T-WTR","water_huddersfield"="E-WTR", "swab_bangkok"="T-dust","swab_huddersfield"="E-dust")
    
    
    
    div_table_rep=ggplot(na.omit(div_table_summ), aes(factor(lot_id, levels=order_lotid),index,fill=origin_country))+ 
      geom_errorbar(aes(ymin=index-sd, ymax=index+sd), color="black", position=position_dodge(0.9), width=0.5) +
      geom_point(size=3, aes(shape=factor(condition, levels=levels_shape), 
                             color=factor(origin_country, levels=levels_country))) +
      scale_fill_manual(values=color_vector2) +
      ylab(alpha_indeces[m]) + 
      xlab("") +
      theme_classic() +
      scale_color_manual(values=color_vector2, name="country of origin") +
      scale_shape_manual(values=shapes_vector2, name="manufacturing conditions") +
      scale_x_discrete(labels=x_labels) +
      facet_row(vars(env), scales = 'free', space = 'free') + theme(legend.position="none")
    
    ggsave (paste(mainDir,"/",subDir4,"/",marker[i],"_",alpha_indeces[m],".pdf", sep=""),div_table_rep, device = cairo_pdf,
            width = 30, height = 15, units = "cm")
    
    assign(paste("rep",alpha_indeces[m],marker[i],sep="_"), div_table_rep)
    
    rm(div_table_rep)
    
    #Statistics analyses -> split in countries and use the t1way which is robust to non normal, non homocedastic data: https://cran.r-project.org/web/packages/WRS2/vignettes/WRS2.pdf
    
    for (n in 1:length(countries)) {
      div_table2=div_table[div_table$material=="tablet",]
      
      #general, without splitting in country
      
      #shapiro.test(residuals(aov(index ~ condition, data=div_table2)))
      #bartlett.test(index ~ condition, data=div_table2)
      
      kruskal=kruskal.test(index ~ condition, data=div_table2) 
      dunn=dunn_test(index ~ condition, p.adjust.method = "holm",data=div_table2)
      
      write.csv(dunn, paste(mainDir,"/",subDir5,"/div_dunn_general",alpha_indeces[m],marker[i],".csv",sep="")) 
      
      div_table3=div_table2[div_table2$origin_country==countries[n] | div_table2$origin_country=="Water Thailand, air England",]
      
      kruskal2=kruskal.test(index ~ condition, data=div_table3) #non 100% sure if it is the best for non-normal data
      dunn2=dunn_test(index ~ condition, p.adjust.method = "holm", data=div_table3) #holm is less restringent than bonferroni
      colnames(dunn2)[colnames(div_table3) == 'index'] <- alpha_indeces[m]
      
      write.csv(dunn2, paste(mainDir,"/",subDir5,"/div_dunn_split",alpha_indeces[m],marker[i],countries[n],".csv",sep="")) 
      
      sink(paste(mainDir,"/",subDir5,"/div_kruskal_dunn",alpha_indeces[m],marker[i],".csv",sep=""), append = TRUE, type = c("output"),split = TRUE)  
      print(Sys.time())
      print(paste("general",marker[i],alpha_indeces[m]))
      print(kruskal)
      print(dunn)
      print(paste("country_split",marker[i],alpha_indeces[m],countries[n]))
      print(kruskal2)
      print(dunn2)
      closeAllConnections()
      
 
    }
  }
  
  assign(paste("shannon",marker[i],sep="_"), shannon)
  assign(paste("richness",marker[i],sep="_"), richness)
  
}  



##################################
# BETA DIVERSITY ANALYSES =======
###############################

samples=c("alltabletsexcipients","alltables","drytablets","wettablets","alltabletsexcipients_noair","wettablets_water","wettablets_water_excipientswg","drytablets_excipientsdc","alltables_noex")
title_plots=c("tablets,excipients & env.controls","all tablets","dry compresion tablets","wet granulation tablets", "tablets,excipients & water","wet granulation tablets & water", "wet granulation tablets, wg excipients & water", "dry compresion tablets & dc excipients","all tablets, excipients sequences removed")

color_vector2=c("Thailand"="red","England"="blue","Water Thailand, air England"="violet", "DC excipient"="gold4", "WG excipient"="gold")
shapes_vector2=c("no water, closed windows"=15,"no water, open windows"=16,"water, open windows"=17,"wet granulation formulation, open windows"=2,"cellulose (DC excipient, USA)"=10,"starch (DC excipient, USA)"=4,
                 "cellulose (WG excipient, USA)"=13,"starch (WG excipient, Thailand)"=8,"tap water"=9,"env. air (dust)"=6)

levels_shape=levels_condition
  
subDir2 <- paste(rarDir,"/figures_ordination",sep="")

##### PCoA ------              

for (i in 1:length(marker)){
  get(paste("rar_count",marker[i],sep="_")) -> count
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  
  mainDir <- paste(path_work,marker[i],sep="/")
  dir.create(file.path(mainDir, subDir2), showWarnings = TRUE)
  
  alltabletsexcipients=count
  alltables=count[rownames(count) %in% metadata_subset$sample.id[metadata_subset$material=="tablet"],]
  drytablets=count[rownames(count) %in% metadata_subset$sample.id[metadata_subset$material=="tablet" & metadata_subset$excipient_formulation=="dc"],]
  wettablets=count[rownames(count) %in% metadata_subset$sample.id[metadata_subset$material=="tablet" & metadata_subset$excipient_formulation=="wg"],] 
  alltabletsexcipients_noair=count[rownames(count) %in% metadata_subset$sample.id[metadata_subset$material!="air"],]
  wettablets_water=count[rownames(count) %in% c(metadata_subset$sample.id[metadata_subset$material=="tablet" & metadata_subset$excipient_formulation=="wg"], metadata_subset$sample.id[metadata_subset$material=="water"]),]
  
  wettablets_water_excipientswg=count[rownames(count) %in% c(metadata_subset$sample.id[metadata_subset$material=="tablet" & metadata_subset$excipient_formulation=="wg"], metadata_subset$sample.id[metadata_subset$material=="water"], metadata_subset$sample.id[metadata_subset$condition=="cellulose (WG excipient, USA)" | metadata_subset$condition=="starch (WG excipient, Thailand)"]),]
  drytablets_excipientsdc=count[rownames(count) %in% c(metadata_subset$sample.id[metadata_subset$material=="tablet" & metadata_subset$excipient_formulation=="dc"], metadata_subset$sample.id[metadata_subset$condition=="cellulose (DC excipient, USA)" | metadata_subset$condition=="starch (DC excipient, USA)"]),]
  
  excipients=count[rownames(count) %in% metadata_subset$sample.id[metadata_subset$material=="excipient"],]
  excipients=excipients[,colSums(excipients) > 0]
  
  alltables_noex=alltables[,!colnames(alltables) %in% colnames(excipients)]
  
   
  for (j in 1:length(samples)) {
    
    pcoa=cmdscale(vegdist(get(samples[j]),method="bray"),eig=TRUE, k=3)
    
    pcoa_df=as.data.frame(pcoa$points[, 1:3])
    
    Eigenvalues <- eigenvals(pcoa) 
    Variance <- Eigenvalues / sum(Eigenvalues) 
    Variance1 <- 100 * signif(Variance[1], 2)
    Variance2 <- 100 * signif(Variance[2], 2)
    Variance3 <- 100 * signif(Variance[3], 2)
    
    
    pcoa_df$lot_id2=metadata_subset$lot_id2[match(rownames(pcoa_df),metadata_subset$sample.id)]
    pcoa_df$origin_country=metadata_subset$origin_country[match(rownames(pcoa_df),metadata_subset$sample.id)]
    pcoa_df$condition=metadata_subset$condition[match(rownames(pcoa_df),metadata_subset$sample.id)]
    pcoa_df$lot_id2_rep=paste(metadata_subset$lot_id2,metadata_subset$replicate)[match(rownames(pcoa_df),metadata_subset$sample.id)]
    
    pcoa_rep1=ggplot(pcoa_df, aes(V1,V2))+ 
      geom_point(size=5, aes(shape=factor(condition, levels=levels_shape), 
                             color=factor(origin_country, levels=levels_country))) +
      scale_color_manual(values=color_vector2, name="Country of origin") +
      scale_shape_manual(values=shapes_vector2, name="Manufacturing conditions")+
      xlab(paste("Axis 1 (", Variance1, "% )")) + 
      ylab(paste("Axis 2 (", Variance2, "% )")) +  
      theme_classic() +
      ggtitle(paste(marker[i],title_plots[j])) 
    
    pcoa_rep2=ggplot(pcoa_df, aes(V1,V3))+ 
      geom_point(size=5, aes(shape=factor(condition, levels=levels_shape), 
                             color=factor(origin_country, levels=levels_country))) +
      scale_color_manual(values=color_vector2, name="Country of origin") +
      scale_shape_manual(values=shapes_vector2, name="Manufacturing conditions")+
      xlab(paste("Axis 1 (", Variance1, "% )")) + 
      ylab(paste("Axis 3 (", Variance3, "% )")) +  
      theme_classic() +
      ggtitle(paste(marker[i],title_plots[j]))
    
    
    pcoa_rep3=ggplot(pcoa_df, aes(V2,V3))+ 
      geom_point(size=5, aes(shape=factor(condition, levels=levels_shape), 
                             color=factor(origin_country, levels=levels_country))) +
      scale_color_manual(values=color_vector2, name="Country of origin") +
      scale_shape_manual(values=shapes_vector2, name="Manufacturing conditions")+
      xlab(paste("Axis 2 (", Variance2, "% )")) + 
      ylab(paste("Axis 3 (", Variance3, "% )")) +  
      theme_classic()+
      ggtitle(paste(marker[i],title_plots[j]))
      
    assign(paste("pcoa_rep1",marker[i],samples[j],sep="_"),pcoa_rep1)
    assign(paste("pcoa_rep2",marker[i],samples[j],sep="_"),pcoa_rep2)
    assign(paste("pcoa_rep3",marker[i],samples[j],sep="_"),pcoa_rep3)
   
  }
}





#### PCoAs 3D of selected samples ####

#samples_comp=c("alltabletsexcipients","alltables","alltables_noex")
samples_comp=c("alltabletsexcipients")

for (i in 1:length(marker[1:2])){
  get(paste("rar_count",marker[i],sep="_")) -> count
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  
  alltabletsexcipients=count
  alltables=count[rownames(count) %in% metadata_subset$sample.id[metadata_subset$material=="tablet"],]
  
  excipients=count[rownames(count) %in% metadata_subset$sample.id[metadata_subset$material=="excipient"],]
  excipients=excipients[,colSums(excipients) > 0]
  
  alltables_noex=alltables[,!colnames(alltables) %in% colnames(excipients)]
  
  mainDir <- paste(path_work,marker[i],sep="/")
  dir.create(file.path(mainDir, subDir2), showWarnings = TRUE)
  
  
  for (j in 1:length(samples_comp)) {
    
    pcoa=cmdscale(vegdist(get(samples_comp[j]),method="bray"),eig=TRUE, k=3)
    
    pcoa_df=as.data.frame(pcoa$points[, 1:3])
    
    Eigenvalues <- eigenvals(pcoa) 
    Variance <- Eigenvalues / sum(Eigenvalues) 
    Variance1 <- 100 * signif(Variance[1], 2)
    Variance2 <- 100 * signif(Variance[2], 2)
    Variance3 <- 100 * signif(Variance[3], 2)
    
    
    pcoa_df$lot_id2=metadata_subset$lot_id2[match(rownames(pcoa_df),metadata_subset$sample.id)]
    pcoa_df$origin_country=metadata_subset$origin_country[match(rownames(pcoa_df),metadata_subset$sample.id)]
    pcoa_df$condition=metadata_subset$condition[match(rownames(pcoa_df),metadata_subset$sample.id)]
    pcoa_df$lot_id2_rep=paste(metadata_subset$lot_id2,metadata_subset$replicate)[match(rownames(pcoa_df),metadata_subset$sample.id)]
    
   
   color_vector=sapply(pcoa_df$origin_country,color_function)
   shapes_vector=sapply(pcoa_df$condition,shape_function)
   

   pdf(paste(mainDir,"/pcoa3D_2_",marker[i],"_",samples_comp[j],".pdf", sep=""))
   plot3D_2=scatter3D(pcoa_df$V1,pcoa_df$V3,pcoa_df$V2,bty = "g", pch = shapes_vector, colvar=NULL,
                      phi = 0, theta = 40,
                      col= color_vector,
                      ticktype = "detailed",
                      xlab=paste("Axis 1 (", Variance1, "% )"), 
                      zlab=paste("Axis 2 (", Variance2, "% )"), 
                      ylab=paste("Axis 3 (", Variance3, "% )"),
                      cex=2,
                      main=ggtitle(paste(marker[i],title_plots[j]))) 
   
   dev.off()
   
  
  assign(paste("pcoa_plot3D_2",marker[i],samples_comp[j],sep="_"),plot3D_2)

   }
}



##### PERMANOVAs ====  

subDir3 <- paste(abundance_threshold,"_rarefaction/tables_permanova",sep="")

for (i in 1:length(marker)) {
  get(paste("rar_count",marker[i],sep="_")) -> count
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  
  mainDir <- paste(path_work,marker[i],sep="/")
  dir.create(file.path(mainDir, subDir3), showWarnings = FALSE)
  
  metadata_subset=metadata_subset[metadata_subset$sample.id %in% rownames(count),]
  
  #permanova only on the tablets
  
  metadata_subset1=metadata_subset[metadata_subset$material=="tablet" & metadata_subset$origin_city!="Water Thailand, air England",]
  count1=count[rownames(count) %in% metadata_subset1$sample.id,]
  rownames(count1)=metadata_subset1$lot_id2[match(rownames(count1),metadata_subset1$sample.id)]
  
  otu_bray=vegdist(count1,method="bray")
  
  #check the dispersion of data
  hist(otu_bray)
  print(anova(betadisper(otu_bray,rownames(count1)))) #Anderson and Walsh (2013) work
  
  origin=metadata_subset1$origin_country[match(rownames(count1),metadata_subset1$lot_id2)]
  water_addition=metadata_subset1$water[match(rownames(count1),metadata_subset1$lot_id2)]
  window_opening=metadata_subset1$window[match(rownames(count1),metadata_subset1$lot_id2)]
  excipient_formulation=metadata_subset1$excipient_formulation[match(rownames(count1),metadata_subset1$lot_id2)] 
  
  if(marker[i]!="trnl") {
    perm=adonis2(otu_bray~origin*window_opening*water_addition*excipient_formulation, permutations=9999)
    write.xlsx(perm, paste(mainDir,"/",subDir3,"/general_permanova_",marker[i],".xlsx", sep=""))
  }
  
  perm2=adonis2(otu_bray~water_addition*origin, permutations=9999)
  write.xlsx(perm2, paste(mainDir,"/",subDir3,"/general_permanova_",marker[i],"2.xlsx", sep=""))
  
  
}


### General Permanova no excipient 

for (i in 1:length(marker)) {
  get(paste("rar_count",marker[i],sep="_")) -> count_ini
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  
  mainDir <- paste(path_work,marker[i],sep="/")
  dir.create(file.path(mainDir, subDir3), showWarnings = FALSE)
  
  alltables=count_ini[rownames(count_ini) %in% metadata_subset$sample.id[metadata_subset$material=="tablet"],]
  
  excipients=count_ini[rownames(count_ini) %in% metadata_subset$sample.id[metadata_subset$material=="excipient"],]
  excipients=excipients[,colSums(excipients) > 0]
  
  alltables_noex=alltables[,!colnames(alltables) %in% colnames(excipients)]
  
  count=alltables_noex
  
  metadata_subset=metadata_subset[metadata_subset$sample.id %in% rownames(count),]
  
  #permanova only on the tablets
  
  metadata_subset1=metadata_subset[metadata_subset$material=="tablet" & metadata_subset$origin_city!="Water Thailand, air England",]
  count1=count[rownames(count) %in% metadata_subset1$sample.id,]
  rownames(count1)=metadata_subset1$lot_id2[match(rownames(count1),metadata_subset1$sample.id)]
  
  otu_bray=vegdist(count1,method="bray")
  
  #check the dispersion of data
  hist(otu_bray)
  print(anova(betadisper(otu_bray,rownames(count1)))) #Anderson and Walsh (2013) work
  
  origin=metadata_subset1$origin_country[match(rownames(count1),metadata_subset1$lot_id2)]
  water_addition=metadata_subset1$water[match(rownames(count1),metadata_subset1$lot_id2)]
  window_opening=metadata_subset1$window[match(rownames(count1),metadata_subset1$lot_id2)]
  
  if(marker[i]!="trnl") {
    perm=adonis2(otu_bray~origin*window_opening*water_addition, permutations=9999)
    write.xlsx(perm, paste(mainDir,"/",subDir3,"/general_permanova_noex",marker[i],".xlsx", sep=""))
  }
  
  perm2=adonis2(otu_bray~water_addition*origin, permutations=9999)
  write.xlsx(perm2, paste(mainDir,"/",subDir3,"/general_permanova_noex",marker[i],"2.xlsx", sep=""))
  
}



#### CROSSPOD TABLES #####

for (i in 1:length(marker[1:2])){
  get(paste("rar_count",marker[i],sep="_")) -> rar_count
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  
  rar_count2=as.data.frame(t(rar_count))
  colnames(rar_count2)=metadata_subset$lot_id2[match(colnames(rar_count2),metadata_subset$sample.id)]
  
  rar_count3=sapply(split.default(rar_count2, names(rar_count2)), rowSums, na.rm = TRUE)
  rar_count3=ifelse(rar_count3 == 0, 0, 1)
  
  order_colnames=metadata_subset[order(factor(metadata_subset$material, levels=levels_mat),factor(metadata_subset$origin_country, levels=levels_country),factor(metadata_subset$condition,levels = levels_condition)),]$lot_id2
  
  rar_count3=rar_count3[,order(match(colnames(rar_count3), order_colnames))]
  
  rar_crossprod_counts=crossprod(rar_count3)
  rar_crossprod_counts_tablets=rar_crossprod_counts[,colnames(rar_crossprod_counts) %in% metadata_subset$lot_id2[metadata_subset$material=="tablet"]]
  
  melt_crossprod=melt(rar_crossprod_counts_tablets)
  melt_crossprod=melt_crossprod[order(match(melt_crossprod$Var1, order_colnames)),]
  melt_crossprod$material=metadata_subset$material[match(melt_crossprod$Var1, metadata_subset$lot_id2)]
  
  cols = c("Var2", "Var1")
  melt_crossprod$value2=melt_crossprod$value
  melt_crossprod$value[as.character(melt_crossprod$Var1)==as.character(melt_crossprod$Var2)] <- NA
  melt_crossprod$value[melt_crossprod$value==0] <-NA
  melt_crossprod$value2[melt_crossprod$value2==0] <-NA
  
  x_labels=c("A-BK1"="T-1","A-BK2"="T-2","B-BK"="T-3","A-HD1"="E-1","A-HD2"="E-2","B-HD"="E-3","B-HD3"="E-4","B-HD2"="E-5",
             "cellulose (DC excipient)"="Cel-DC", "starch (DC excipient)"="Str-DC","cellulose (WG excipient)"="Cel-WG", "starch (WG excipient)"="Str-WG",
             "water_bangkok"="T-WTR","water_huddersfield"="E-WTR", "swab_bangkok"="T-dust","swab_huddersfield"="E-dust")
  
  
  rep_cross=ggplot(melt_crossprod, aes(Var2, factor(Var1, levels=rev(as.character(unique(melt_crossprod$Var1)))))) +
    geom_point(aes(size = value, color=value)) + theme_bw() + xlab("") + ylab("") + 
    scale_size_continuous(range=c(5,15), guide='none') +
    scale_colour_viridis_c(name="Share ASVs", alpha=0.5) +
    geom_text(aes(label = value2)) + 
    scale_x_discrete(position = "top",labels=x_labels) +
    scale_y_discrete(labels=x_labels)+
    facet_grid(factor(material, levels=c("tablet","excipient","water","air"))~.,drop = TRUE,scales = "free", space = "free") +
    ggtitle(paste(marker[i])) + theme(plot.title = element_text(hjust = 0.5))
  
  assign(paste("rep_cross",marker[i],sep="_"),rep_cross)
  
}  




### SHANNON, PCOAS AND CROSSPOD JOIN REPRESENTATIONS ######

size_text=15

theme_rep=theme(panel.background= element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.background =element_rect(fill="white"),
                axis.text.y = element_text(size = size_text, color="black"),
                axis.text.x = element_text(size = size_text, color="black", angle=0),
                axis.title = element_text(size = size_text),
                legend.text = element_text(size = size_text),
                legend.title = element_text(size =size_text),
                legend.position = "none",
                axis.line = element_line(colour="black", linewidth= 0.5, linetype = 1), 
                axis.ticks = element_line(colour="black", linewidth= 0.5),
                axis.ticks.length = unit(0.1, "cm"),
                axis.text.y.right = element_text(color = "red"),
                axis.title.y.right = element_text(color= "red")
)


div_figures_together=ggarrange(rep_richness_16S + theme(axis.text.x=element_blank(), strip.text.x = element_blank()) + ylab("16S richness"),
                               rep_shannon_16S+theme(strip.background = element_blank(), strip.text.x = element_blank(),axis.text.x=element_blank())+ ylab("16S Shannon index"), 
                               rep_richness_18S+theme(strip.background = element_blank(), strip.text.x = element_blank(),axis.text.x=element_blank()) + ylab("18S richness"), 
                               rep_shannon_18S+theme(strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x=element_blank()) + ylab("18S Shannon index"),
                               nrow=4, common.legend = TRUE, legend="none",align = "v",widths = c(1.5,2))


ggsave (paste(path_work,"/richness_shannon_allmarkers.pdf", sep=""),div_figures_together, device = cairo_pdf,
        width = 30, height= 25, units = "cm")


pcoa_figures_together1=ggarrange(pcoa_rep1_16S_alltabletsexcipients + theme_rep, pcoa_rep1_18S_alltabletsexcipients + theme_rep, 
                               ncol=2, common.legend = TRUE, legend="right") 

pcoa_figures_together2=ggarrange(pcoa_rep2_16S_alltabletsexcipients, pcoa_rep3_16S_alltabletsexcipients, 
                                 pcoa_rep2_18S_alltabletsexcipients, pcoa_rep3_18S_alltabletsexcipients,
                                 ncol=4, common.legend = TRUE, legend="none") 

pcoa_figures_together3=ggarrange(pcoa_rep1_16S_drytablets+theme_rep, pcoa_rep1_16S_wettablets+theme_rep, 
                                 pcoa_rep1_18S_drytablets+theme_rep, pcoa_rep1_18S_wettablets+theme_rep,
                                 ncol=4, common.legend = TRUE, legend="none") 

ggsave (paste(path_work,"/pcoa_allmarkers1.pdf", sep=""),pcoa_figures_together1, device = cairo_pdf,
        width = 40, height = 15, units = "cm")

ggsave (paste(path_work,"/pcoa_allmarkers2.pdf", sep=""),pcoa_figures_together2, device = cairo_pdf,
        width = 40, height = 10, units = "cm")

ggsave (paste(path_work,"/pcoa_allmarkers3.pdf", sep=""),pcoa_figures_together3, device = cairo_pdf,
        width = 40, height = 10, units = "cm")


crosspods_together=ggarrange(rep_cross_16S, rep_cross_18S, 
                             ncol=2, common.legend = FALSE, legend="right") 



##### DESEQ (no excipients) ######## 

subDir10 <-"tables_deseq/asv_level"
samples2=c("drytablets","wettablets")


for (i in 1:length(marker)) {
  print(marker[i])
  get(paste("count_filt",marker[i],sep="_")) -> count_ini
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  get(paste("rep-seqs_merged",marker[i],sep = "_")) -> rep_seqs
  
  alltables=count_ini[,colnames(count_ini) %in% metadata_subset$sample.id[metadata_subset$material=="tablet"]]
  
  excipients=count_ini[,colnames(count_ini) %in% metadata_subset$sample.id[metadata_subset$material=="excipient"]]
  excipients=excipients[rowSums(excipients) > 0,]
  
  alltables_noex=alltables[!rownames(alltables) %in% rownames(excipients),]
  
  count=alltables_noex
  
  metadata_subset=metadata_subset[metadata_subset$material=="tablet" & metadata_subset$origin_country!="Water Thailand, air England",]
  metadata_subset=metadata_subset[metadata_subset$condition!="wet granulation formulation, open windows",]
  
  metadata_subset$origin_country=factor(metadata_subset$origin_country, levels=c("Thailand","England"))
  metadata_subset$condition=factor(metadata_subset$condition, levels=c("no water, closed windows","no water, open windows","water, open windows"))
  
  count=count[,colnames(count) %in% metadata_subset$sample.id]
  count=count[,colSums(count) > 1000]
  
  metadata_subset=metadata_subset[metadata_subset$sample.id %in% colnames(count),]
  
  count <- count[, order(factor(colnames(count), levels = metadata_subset$sample.id))]
  
  
  #### deseq on whole set
  
  print(all(colnames(count) == metadata_subset$sample.id))
  print(ncol(count) == nrow(metadata_subset))
  
  
  dds_all <- DESeqDataSetFromMatrix(countData=count, 
                                    colData=metadata_subset, 
                                    design= ~ origin_country + origin_country:water)
  
  
  dds_test_all <- DESeq(dds_all, sfType ="poscounts")
  dds_results_all= results(dds_test_all, contrast=c("origin_country", "England", "Thailand"),  pAdjustMethod="BH", alpha = 0.05, cooksCutoff=FALSE) # in the contrast argument we are telling to use England as nominator and Thailand as denominator -> log2 (England/Thailand)
  
  dds_counts=counts(dds_test_all, normalized=TRUE)
  assign(paste("dds_counts_noex",marker[i],"_all",sep=""),dds_counts)
  
  plotDispEsts(dds_test_all)
  
  sink(paste(mainDir,"/",subDir10,"/summary_deseq_all_noex_",marker[i],".csv",sep=""), append = FALSE, type = c("output"),split = TRUE)  
  print(Sys.time())
  print(summary(dds_results_all))
  closeAllConnections()
  
  dds_results_all2=as.data.frame(dds_results_all)
  dds_results_all2 <- na.omit(dds_results_all2[dds_results_all2$padj < 0.05, ])
  dds_results_all2=dds_results_all2[order(-dds_results_all2$log2FoldChange),]
  dds_results_all2$Feature.ID=row.names(dds_results_all2)
  
  writeXStringSet(get("rep_seqs")[rownames(dds_results_all2)],paste(mainDir,"/",subDir10,"/fasta_noex",marker[i],"_all.fa", sep=""))
  
  
  ### deseq splitting
  
  drytablets=count[,colnames(count) %in% metadata_subset$sample.id[metadata_subset$material=="tablet" & metadata_subset$water=="dry"]]
  wettablets=count[,colnames(count) %in% metadata_subset$sample.id[metadata_subset$material=="tablet" & metadata_subset$water=="wet"]]
  
  
  for (m in 1:length(samples2)) {
    print(samples2[m])
    
    count1=get(samples2[m])
    
    metadata_subset1=metadata_subset[metadata_subset$sample.id %in% colnames(count1),]
    metadata_subset1=metadata_subset1[order(metadata_subset1$origin_country,metadata_subset1$condition),]
    
    count1 <- count1[, order(factor(colnames(count1), levels = metadata_subset1$sample.id))]

    print(all(colnames(count1) == metadata_subset1$sample.id))
    print(ncol(count1) == nrow(metadata_subset1))
    print(is.data.frame(metadata_subset1))
    print(is.data.frame(count1))
    
    
    dds <- DESeqDataSetFromMatrix(countData=count1, 
                                  colData=metadata_subset1, 
                                  design= ~ origin_country)
    
    dds_test <- DESeq(dds, sfType ="poscounts")
    dds_results= results(dds_test, contrast=c("origin_country", "England", "Thailand"),  pAdjustMethod="BH", alpha = 0.05, cooksCutoff=FALSE) # in the contrast argument we are telling to use England as nominator and Thailand as denominator -> log2 (England/Thailand)
    
    
    sink(paste(mainDir,"/",subDir10,"/summary_deseq_noex",marker[i],"_",samples2[m],".csv",sep=""), append = FALSE, type = c("output"),split = TRUE)  
    print(Sys.time())
    print(summary(dds_results))
    closeAllConnections()
    
    
    dds_results2=as.data.frame(dds_results)
    dds_results2 <- na.omit(dds_results2[dds_results2$padj < 0.05, ])
    dds_results2=dds_results2[order(-dds_results2$log2FoldChange),]
    dds_results2$Feature.ID=row.names(dds_results2)
    
    writeXStringSet(get("rep_seqs")[rownames(dds_results2)],paste(mainDir,"/",subDir10,"/fasta_noex",marker[i],"_",samples2[m],".fa", sep=""))
    write.xlsx(dds_results2, paste(mainDir,"/",subDir10,"/deseq_table_noex",marker[i],"_",samples2[m],".xlsx",sep=""))
    
    assign(paste("dds_results_noex",marker[i],"_",samples2[m],sep=""),dds_results2)
    
  }
  
  write.xlsx(dds_results_all2, paste(mainDir,"/",subDir10,"/deseq_table_noex",marker[i],"_all.xlsx",sep=""))
  assign(paste("dds_results_noex",marker[i],"_all",sep=""),dds_results_all2)
}  



###### ANACOMB (no excipients) #######

for (i in 1:length(marker)) {
  print(marker[i])
  get(paste("count_filt",marker[i],sep="_")) -> count_ini
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  
  alltables=count_ini[,colnames(count_ini) %in% metadata_subset$sample.id[metadata_subset$material=="tablet"]]
  
  excipients=count_ini[,colnames(count_ini) %in% metadata_subset$sample.id[metadata_subset$material=="excipient"]]
  excipients=excipients[rowSums(excipients) > 0,]
  
  alltables_noex=alltables[!rownames(alltables) %in% rownames(excipients),]
  
  count=alltables_noex

  #ancombc only on the tablets
  
  metadata_subset1=metadata_subset[metadata_subset$material=="tablet" & metadata_subset$origin_country!="Water Thailand, air England",]
  metadata_subset1=metadata_subset1[metadata_subset1$condition!="wet granulation formulation, open windows",]
  
  count1=count[,colnames(count) %in% metadata_subset1$sample.id]
  
  metadata_subset1$origin_country=factor(metadata_subset1$origin_country, levels=c("Thailand","England"))
  metadata_subset1$condition=factor(metadata_subset1$condition, levels=c("no water, closed windows","no water, open windows","water, open windows","wet granulation formulation, open windows"))
  
  otu=otu_table(count1, taxa_are_rows = TRUE)
  sampledata=sample_data(metadata_subset1)
  
  physeq1 = merge_phyloseq(otu, sampledata)
  
  out = ancombc2(data = physeq1, assay_name = "counts", 
                 fix_formula = "origin_country + condition",
                 p_adj_method = "holm", verbose = TRUE)
  
  res=out$res
  res2=as.data.frame(do.call(cbind, res))
  res2=res2[res2$p_origin_countryEngland < 0.1,]
  
  #writeXStringSet(get(paste("rep-seqs-dada2",marker[i],sep = "_"))[res2$lfc.taxon],paste(path_work,"/",directories[i],"/ancombc/asv_level/fasta_",marker[i],"_",samples2[m],".fa", sep=""))
  
  assign(paste("ancombc_results_noex",marker[i],sep=""),res2)
  
}


#no ancombc results for trnl



####### Select/mark ASVs that appeared differentially abundant in Deseq and Ancombc, and record exclusive abundance in England or Thai samples ######

deseq_ancomb_comb_16S=merge(dds_results_noex16S_all, ancombc_results_noex16S, by.x="Feature.ID", by.y="taxon", all=TRUE) ## choose no excipients one
deseq_ancomb_comb_18S=merge(dds_results_noex18S_all, ancombc_results_noex18S, by.x="Feature.ID", by.y="taxon", all=TRUE)
deseq_ancomb_comb_trnl=dds_results_noextrnl_all #no ancombc results available

write.xlsx(deseq_ancomb_comb_18S,paste(path_work,"/18S/tables_deseq_ancomb/deseq_ancomb_comb_18S_all.xlsx",sep = ""))

  
for (i in 2:length(marker)) {
  print(marker[i])
  get(paste("rar_count",marker[i],sep="_")) -> rar_count
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  get(paste("deseq_ancomb_comb",marker[i],sep = "_")) -> deseq_ancomb
  get(paste("rep-seqs_merged",marker[i],sep = "_")) -> rep_seqs 
  
  
  mainDir <- paste(path_work,marker[i],sep="/")
  
  metadata_subset=metadata_subset[metadata_subset$material=="tablet" & metadata_subset$origin_country!="Water Thailand, air England",]
  
  order_sample_id=metadata_subset$sample.id[order(factor(metadata_subset$lot_id2, levels=c("A-BK1","A-BK2","B-BK","A-HD1","A-HD2","B-HD","B-HD2","B-HD3")), metadata_subset$replicate)]
  
  rar_count2=as.data.frame(t(rar_count))
  rar_count2=rar_count2[,colnames(rar_count2) %in% metadata_subset$sample.id]
  rar_count2=rar_count2[,order(factor(colnames(rar_count2),levels=order_sample_id))]
  
  rar_count3=ifelse(rar_count2 == 0, 0, 1)
  
  asvs_select_england=rar_count3[rowSums(rar_count3[,colnames(rar_count3) %in% metadata_subset$sample.id[metadata_subset$origin_country=="England"]]) > 1,] #select asvs that appear in at least two samples
  asvs_select_england=rownames(asvs_select_england)                    
  
  asvs_select_thailand=rar_count3[rowSums(rar_count3[,colnames(rar_count3) %in% metadata_subset$sample.id[metadata_subset$origin_country=="Thailand"]]) > 1,] 
  asvs_select_thailand=rownames(asvs_select_thailand)                    
  
  rar_count4=rar_count2[rownames(rar_count2) %in% c(asvs_select_england,asvs_select_thailand),]  
  
  if(marker[i]=="18S"){
  rar_count4=rar_count4[rownames(rar_count4) %in% c(tax_18S_asvs),]}  
  
  colnames(rar_count4)=metadata_subset$lot_id2[match(colnames(rar_count4),metadata_subset$sample.id)]
  
  deseq_ancomb=deseq_ancomb[,colnames(deseq_ancomb) %in% c("Feature.ID","log2FoldChange","padj","lfc_origin_countryEngland","p_origin_countryEngland","q_origin_countryEngland")]
  deseq_ancomb$da_asvs=asvs_rep_names$rep_name[match(deseq_ancomb$Feature.ID,asvs_rep_names$Feature.ID)]
  
  deseq_ancomb_rar_count=merge(deseq_ancomb, rar_count4, by.x="Feature.ID", by.y="row.names", all=FALSE)
  
  if(marker[i]=="18S"){
  deseq_ancomb_rar_count$lfc_origin_countryEngland=as.numeric(deseq_ancomb_rar_count$lfc_origin_countryEngland)}
  
  deseq_ancomb_rar_count=deseq_ancomb_rar_count[order(deseq_ancomb_rar_count$log2FoldChange),]
  order_asvs=deseq_ancomb_rar_count$da_asvs
  
  write.csv2(deseq_ancomb_rar_count, paste(marker[i],"/blast_deseq_ancomb/deseq_ancomb_select_",marker[i],".csv", sep=""))
  writeXStringSet(get("rep_seqs")[deseq_ancomb_rar_count$Feature.ID],paste(marker[i],"/blast_deseq_ancomb/deseq_ancomb_asvs_",marker[i],".fa", sep=""))

}


#### SELECT ASVS BASED ON RAREFIED ABUNDANCE rather than in da ('manual' selection) ######

for (i in 2:length(marker)) {
  print(marker[i])
  get(paste("rar_count",marker[i],sep="_")) -> count
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  get(paste("deseq_ancomb_comb",marker[i],sep = "_")) -> deseq_ancomb
  get(paste("rep-seqs_merged",marker[i],sep = "_")) -> rep_seqs 
  
  metadata_subset=metadata_subset[metadata_subset$origin_country!="Water Thailand, air England",]
  
  alltables=count[rownames(count) %in% metadata_subset$sample.id[metadata_subset$material=="tablet"],]
  
  excipients=count[rownames(count) %in% metadata_subset$sample.id[metadata_subset$material=="excipient"],]
  excipients=excipients[,colSums(excipients) > 0]
  
  alltables_noex=alltables[,!colnames(alltables) %in% colnames(excipients)]
  
  order_sample_id=metadata_subset$sample.id[order(factor(metadata_subset$lot_id2, levels=c("A-BK1","A-BK2","B-BK","A-HD1","A-HD2","B-HD","B-HD2","B-HD3")), metadata_subset$replicate)]
  
  rar_count2=as.data.frame(t(alltables_noex))
  rar_count2=rar_count2[,colnames(rar_count2) %in% metadata_subset$sample.id]
  rar_count2=rar_count2[,order(factor(colnames(rar_count2),levels=order_sample_id))]
  

  if(marker[i]=="18S"){
    rar_count2=rar_count2[rownames(rar_count2) %in% c(tax_18S_asvs),]}  
  
  
  rar_count3=ifelse(rar_count2 == 0, 0, 1)
  

  asvs_select_england_exclusive=rar_count3[rowSums(rar_count3[,colnames(rar_count3) %in% metadata_subset$sample.id[metadata_subset$origin_country=="England"]]) > 1 & rowSums(rar_count3[,colnames(rar_count3) %in% metadata_subset$sample.id[metadata_subset$origin_country=="Thailand"]]) == 0,] #select asvs that exclusively appears in a country and not other
  asvs_select_england_exclusive=rownames(asvs_select_england_exclusive)                    
  
  if (marker[i]!="trnl") {
  asvs_select_thailand_exclusive=rar_count3[rowSums(rar_count3[,colnames(rar_count3) %in% metadata_subset$sample.id[metadata_subset$origin_country=="Thailand"]]) > 1 & rowSums(rar_count3[,colnames(rar_count3) %in% metadata_subset$sample.id[metadata_subset$origin_country=="England"]]) == 0,] 
  } else {
  asvs_select_thailand_exclusive=rar_count3[rowSums(rar_count3[,colnames(rar_count3) %in% metadata_subset$sample.id[metadata_subset$origin_country=="Thailand"]]) > 0 & rowSums(rar_count3[,colnames(rar_count3) %in% metadata_subset$sample.id[metadata_subset$origin_country=="England"]]) == 0,]} 
  
  asvs_select_thailand_exclusive=rownames(asvs_select_thailand_exclusive)                    
  
  write.csv(asvs_select_england_exclusive,paste(marker[i],"/blast_manual/asvs_select_england_exclusive_",marker[i],".csv", sep = ""))
  write.csv(asvs_select_thailand_exclusive,paste(marker[i],"/blast_manual/asvs_select_thailand_exclusive_",marker[i],".csv", sep = ""))
  
  asvs_select_england=rar_count3[rowSums(rar_count3[,colnames(rar_count3) %in% metadata_subset$sample.id[metadata_subset$origin_country=="England"]]) > 1,] #select asvs that appear in at least two samples of a same country
  asvs_select_england=rownames(asvs_select_england)                    
  
  if (marker[i]!="trnl") {asvs_select_thailand=rar_count3[rowSums(rar_count3[,colnames(rar_count3) %in% metadata_subset$sample.id[metadata_subset$origin_country=="Thailand"]]) > 1, ,drop=FALSE] 
  } else {asvs_select_thailand=rar_count3[rowSums(rar_count3[,colnames(rar_count3) %in% metadata_subset$sample.id[metadata_subset$origin_country=="Thailand"]]) > 0, ,drop=FALSE]} 
  
  asvs_select_thailand=rownames(asvs_select_thailand)                    
  
  
  colnames(rar_count2)=metadata_subset$lot_id2[match(colnames(rar_count2),metadata_subset$sample.id)]
  
  count_england=rar_count2[rownames(rar_count2) %in% c(asvs_select_england),]  
  count_england=count_england[order(-rowSums(count_england)),]
  count_england$sum=rowSums(count_england)
  count_england$asv_country=rep("E",nrow(count_england))
  count_england$asvs=asvs_rep_names$rep_name[match(rownames(count_england),asvs_rep_names$Feature.ID)]
  count_england$exclusive=ifelse(rownames(count_england) %in% asvs_select_england_exclusive,"exclusive_E","non-exclusive")
  count_england=count_england[count_england$exclusive=="exclusive_E",]
  count_england=count_england[count_england$sum > 100,]
  
  count_thailand=rar_count2[rownames(rar_count2) %in% c(asvs_select_thailand),]  
  count_thailand=count_thailand[order(-rowSums(count_thailand)),]
  count_thailand$sum=rowSums(count_thailand)
  count_thailand$asv_country=rep("T",nrow(count_thailand))
  count_thailand$asvs=asvs_rep_names$rep_name[match(rownames(count_thailand),asvs_rep_names$Feature.ID)]
  count_thailand$exclusive=ifelse(rownames(count_thailand) %in% asvs_select_thailand_exclusive,"exclusive_T","non-exclusive")
  count_thailand=count_thailand[count_thailand$exclusive=="exclusive_T",]
  count_thailand=count_thailand[count_thailand$sum > 100,]
  
  count_select_all=rbind(count_england,count_thailand)
  count_select_all=count_select_all[rowSums(is.na( count_select_all)) != ncol(count_select_all), ]
  count_select_all$Feature.ID=rownames(count_select_all)
  
  count_select_all2=count_select_all[!colnames(count_select_all) %in% c("sum")]
   
  write.xlsx(count_select_all2, paste(marker[i],"/blast_manual/manual_select_asvs_",marker[i],".xlsx", sep=""))
  writeXStringSet(get("rep_seqs")[rownames(count_select_all2)],paste(marker[i],"/blast_manual/manual_select_asvs_",marker[i],".fa", sep=""))
  
}  
  

#### save seqs for blast #####

for (i in 2:length(marker)){
  get(paste("rep-seqs_merged",marker[i],sep = "_")) -> rep_seqs 
  get(paste("count_filt",marker[i],sep = "_")) -> count_filt
  get(paste("metadata",marker[i],sep = "_")) -> metadata_subset
  
  alltables=count[rownames(count) %in% metadata_subset$sample.id[metadata_subset$material=="tablet"],]
  
  excipients=count[rownames(count) %in% metadata_subset$sample.id[metadata_subset$material=="excipient"],]
  excipients=excipients[,colSums(excipients) > 0]
  
  alltables_noex=alltables[,!colnames(alltables) %in% colnames(excipients)]
  
  select_asvs=colnames(alltables_noex)
  writeXStringSet(get("rep_seqs")[rownames(count_filt)],paste(marker[i],"/blast_all/fasta_",marker[i],"_all.fa", sep=""))
}




#### PREDICTIVE MODELS ######
#DPAC deseq features from 16S and 18S ########

deseq_ancomb=rbind(deseq_ancomb_comb_16S,deseq_ancomb_comb_18S)

metadata_subset=metadata[metadata$material=="tablet",]
metadata_subset$DNA_aliquot_id[metadata_subset$sample.id=="30610OR0029L01"] <- "PLCB-RVQ5-T_4"
metadata_subset$DNA_aliquot_id[metadata_subset$sample.id=="30610OR0037L01"] <- "PLCB-G0ZS-T_4"
metadata_subset=metadata_subset[order(metadata_subset$origin_country),]

metadata_subset=metadata_subset %>% mutate(lot_id3=x_labels[lot_id2], lot_id2_rep=paste(lot_id3, replicate)) 

count_16s=rar_count_16S[rownames(rar_count_16S) %in% metadata_subset$sample.id,]
count_18s=rar_count_18S[rownames(rar_count_18S) %in% metadata_subset$sample.id,]

rownames(count_16s)=metadata_subset$lot_id2_rep[match(rownames(count_16s),metadata_subset$sample.id)]
rownames(count_18s)=metadata_subset$lot_id2_rep[match(rownames(count_18s),metadata_subset$sample.id)]

rownames(count_16s) %in% rownames(count_18s)

table=merge(count_16s,count_18s, by="row.names", all=TRUE)
rownames(table)=table$Row.names
table=table[,-1]
table[is.na(table)] <- 0    

metadata_subset2=metadata_subset[metadata_subset$origin_country!="Water Thailand, air England",]
metadata_sup=metadata_subset[metadata_subset$origin_country=="Water Thailand, air England",]

supp_ind=table[rownames(table) %in% metadata_sup$lot_id2_rep,,drop = FALSE]

table=table[rownames(table) %in% metadata_subset2$lot_id2_rep,]
table=table[,colnames(table) %in% deseq_ancomb$Feature.ID]
table=table[,colSums(table)>100]
table=table[rowSums(table)>0,]
table=table[order(match(rownames(table), metadata_subset2$lot_id2_rep)),]


### dapc on combined dataset

dapc_model=dapc(table, country,n.pca = 13,n.da=1) #chose 13 to retain >90% total variance
pred.sup <- predict.dapc(dapc_model, newdata=supp_ind)

write.csv(pred.sup, paste(path_work,"/dpca_prediction_16S&18S.csv", sep="") )

pdf(file=paste(path_work,"/dpca_metrics_16S&18S.pdf", sep=""),width=10, height=7)
layout(matrix(c(1,1,2,3), nrow=2, byrow = TRUE))
scatter(dapc_model, legend=FALSE)
assignplot(dapc_model)
compoplot(dapc_model, show.lab = TRUE, posi = "topright", col.pal=c("lightblue","tomato"), legend=FALSE)
dev.off()


dapc_den=data.frame(rbind(dapc_model$ind.coord,pred.sup$ind.scores))

dapc_den$lot_id2_rep=row.names(dapc_den)
dapc_den$origin_country=metadata_subset$origin_country[match(rownames(dapc_den),metadata_subset$lot_id2_rep)]
dapc_den$condition=metadata_subset$condition[match(rownames(dapc_den),metadata_subset$lot_id2_rep)]
dapc_den_levels=dapc_den$lot_id2_rep[order(dapc_den$origin_country)]


rep=ggplot(dapc_den,aes(factor(lot_id2_rep, levels=dapc_den_levels), LD1)) +
  coord_flip()+
  geom_point(size=3, aes(shape=factor(condition, levels=levels_shape), 
                         color=factor(origin_country, levels=levels_country))) +
  scale_color_manual(values=color_vector2, name="Country of origin") +
  scale_shape_manual(values=shapes_vector2, name="Manufacturing conditions")+
  xlab("")+
  ylab("Discriminant function scores")+
  theme_classic()+ theme(legend.position="none") +
  geom_hline(yintercept = dapc_model$grp.coord[1], color="blue")+
  geom_hline(yintercept = dapc_model$grp.coord[2], color="red") +
  ggtitle("DPCA-discriminant function scores for 16S and 18S da-ASVs")

ggsave (paste(path_work,"/dpca_scores_16S&18S.pdf", sep=""),rep, device = cairo_pdf,
        width = 30, height = 15, units = "cm")


loocv_predictions <- data.frame(sample=character(),prediction=character(),England_post=numeric(), Thailand_post=numeric(), real=character())

for(k in 1:nrow(table)) {
  # Leave one out: Training data (all rows except the i-th row)
  train <- table[-k,]
  predict <- table[k, ,drop=FALSE]
  
  country2=metadata_subset2$origin_country[match(rownames(train),metadata_subset2$lot_id2_rep)]
  
  dapc_model2=dapc(train, country2, n.pca=npca_dapc[i],n.da=1)
  loocv_predictions[k,"sample"] <- rownames(predict)
  loocv_predictions[k,"prediction"] <- predict.dapc(dapc_model2, newdata=predict)[[1]]
  loocv_predictions[k,"prediction"] <-ifelse(loocv_predictions[k,"prediction"]==1,"England","Thailand")
  loocv_predictions[k,c("England_post","Thailand_post")] <- round(predict.dapc(dapc_model2, newdata=predict)$posterior,4)
  loocv_predictions[k,"real"] <- metadata_subset2$origin_country[match(rownames(predict),metadata_subset2$lot_id2_rep)]
}

write.csv(loocv_predictions,paste(path_work,"/dpca_loocv_predictions_table_16S&18S.csv", sep=""))

sink(paste(path_work,"/dpca_loocv_accurary_16S&18S.csv", sep=""), append = FALSE, type = c("output"),split = TRUE)  
paste("DAPC correct predictions for total of",nrow(loocv_predictions), "samples for 16S and 18S da-ASVs")
sum(loocv_predictions$prediction == loocv_predictions$real)
paste("model accuracy")
paste((round(sum(loocv_predictions$prediction == loocv_predictions$real) / nrow(loocv_predictions),2))*100, "%")
closeAllConnections()

loocv_predictions2 = loocv_predictions[,colnames(loocv_predictions) %in% c("sample","England_post","Thailand_post")]
loocv_predictions2[,c(2,3)]=loocv_predictions2[,c(2,3)]*100

loocv_predictions2=melt(loocv_predictions2)
loocv_predictions2=loocv_predictions2[order(loocv_predictions2$variable,loocv_predictions2$sample),]

loocv_predictions2[,3]=100-loocv_predictions2[,3]

assign_rep=ggplot(loocv_predictions2, aes(x=sample,y=value,fill=variable))+
  geom_col(alpha=0.7)+
  scale_fill_manual(values=c("red","blue"))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(name="Membership prob. (%)")+
  theme_classic2() +
  theme(
    axis.line.x = element_blank(),   # Remove x-axis line
    axis.ticks.x = element_blank(),   # Remove x-axis ticks
    axis.line.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0,'lines'))+
  geom_segment(aes(x = -1.5, xend = -1.5, y = 0, yend = 100), color = "black")

ggsave("assign_prob.png", assign_rep, units = "cm",width = 10, height = 5,
       dpi = 300)







