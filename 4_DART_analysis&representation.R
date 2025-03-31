##### set working directory and paths for data and figures ------

#from office
path_work="K:/ogden_grp/Conservation_Science/Conservation_Genetics/GENOME_PROJECTS/Carla/dart/simulated/huds_vs_bk"

setwd(path_work)

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

packages <- c('openxlsx','ggbiplot','mt','dplyr','gplots','ggpubr','tidyr','plot3D','reshape2','stringr','vegan')

ipak(packages)


#load variables and to represent the ordinations ######

names_correspondances_hudd_vs_bk=read.xlsx("names_correspondances_hudd_vs_bk.xlsx",sheet = "Sheet1")

ordination=c("pca","kda")

color_vector=c("Thailand"="red","England"="blue","Water Thailand, air England"="violet")
shapes_vector=c("no water, closed windows"=15,"no water, open windows"=16,"water, open windows"=17,"wet granulation formulation, open windows"=2)

pca_cov_15mmu_scores=read.csv("pca_scores_allgroups_cov_15mmu.csv",sep=",",header=FALSE)
kda_cov_15mmu_scores=read.csv("kda_scores_allgroups_cov_15mmu.csv",sep=",",header=FALSE)

selected_peaks=read.xlsx("masses_anovaselect_allgroups.xlsx",sheet = "Sheet1")
selected_peaks$count=rep(0.5,nrow(selected_peaks))
colnames(selected_peaks)[1] <-'mass_selected'
rownames(selected_peaks)=selected_peaks$mass_selected

pca_variances=c(49.11,17.26,10.09) #read from the excel files saved as pca/kda variances allgroups
kda_variances=c(77.54,20.36,1.11)

color_function <- function (x) {if(x=="Thailand") "red"
  else if(x=="England")"blue" else if (x=="Water Thailand, air England") "violet"}

shape_function <- function (x) {if(x=="no water, closed windows") 15
  else if(x=="no water, open windows") 16 
  else if (x=="water, open windows") 17 
  else if (x=="wet granulation formulation, open windows") 2
}

#PCA and KDA on covariances ####

for(i in 1:length(ordination)) {

  get(paste(ordination[i],"_cov_15mmu_scores",sep = "")) -> table_ordination
  get(paste(ordination[i],"_variances",sep = "")) -> variances
  
  table_ordination=table_ordination %>% separate_wider_delim(V1, delim = ",", names = c("Index", "Component1","Component2","Component3")) %>% mutate_if(is.character, as.numeric)
  table_ordination=merge(names_correspondances_hudd_vs_bk,table_ordination, by="Index") #

  Variance1=variances[1]
  Variance2=variances[2]
  Variance3=variances[3]
  
  ord_rep_1_2=ggplot(table_ordination, aes(Component1,Component2))+ 
  geom_point(size=5, aes(shape=condition, 
                         color=country)) +
  scale_color_manual(values=color_vector, name="Country of origin") +
  scale_shape_manual(values=shapes_vector, name="Manufacturing conditions")+
  xlab(paste("PC1 (", round(Variance1,0), "% )")) + 
  ylab(paste("PC2 (", round(Variance2,0), "% )")) +  
  theme_classic()

  ord_rep_1_3=ggplot(table_ordination, aes(Component1,Component3))+ 
    geom_point(size=3, aes(shape=condition, 
                           color=country)) +
    scale_color_manual(values=color_vector, name="Country of origin") +
    scale_shape_manual(values=shapes_vector, name="Manufacturing conditions")+
    xlab(paste("PC1 (", round(Variance1,0), "% )")) + 
    ylab(paste("PC3 (", round(Variance3,0), "% )")) +  
    theme_classic()
  
  ord_rep_2_3=ggplot(table_ordination, aes(Component2,Component3))+ 
    geom_point(size=3, aes(shape=condition, 
                           color=country)) +
    scale_color_manual(values=color_vector, name="Country of origin") +
    scale_shape_manual(values=shapes_vector, name="Manufacturing conditions")+
    xlab(paste("PC2 (", round(Variance2,0), "% )")) + 
    ylab(paste("PC3 (", round(Variance3,0), "% )")) +  
    theme_classic()
  
  allplots=ggarrange(ord_rep_1_2, ord_rep_1_3, ord_rep_2_3, nrow = 1, common.legend = T,legend="right")
  
  pdf(paste(path_work,"/",ordination[i],"_covariance_15mmu_tolerance.pdf", sep=""))
  
  print(annotate_figure(allplots, top = text_grob(paste(ordination[i],"covariance, 15mmu tolerance"),face = "bold", size = 14)))
  
  dev.off()
 
  color_vector2=sapply(table_ordination$country,color_function)
  shapes_vector2=sapply(table_ordination$condition,shape_function)
  
  pdf(paste(path_work,"/",ordination[i],"_covariance_15mmu_3d.pdf", sep=""))
  scatter3D(table_ordination$Component1,table_ordination$Component3,table_ordination$Component2,bty = "g", pch = shapes_vector2, colvar=NULL,
            phi = 0, theta = 40,
            col= color_vector2,
            ticktype = "detailed",
            xlab=paste("Axis 1 (", round(Variance1, 0), "% )"), 
            zlab=paste("Axis 2 (", round(Variance2, 0), "% )"), 
            ylab=paste("Axis 3 (", round(Variance3, 0), "% )"),
            cex=2, type="p")
  dev.off()
  
  
  scatterplot3d(table_ordination$Component1,table_ordination$Component3,table_ordination$Component2, 
                color=color_vector2,
                pch=shapes_vector2,
                angle=180,
                xlab="Axis 1", 
                zlab="Axis 2", 
                ylab="Axis 3",
                grid=TRUE, 
                box=TRUE,
                scale.y = 0.5,
                cex.symbols=1,
                type="p",tick.marks = FALSE)
  
  
}



## load heatmap and represent as colmaps ##### 

heatmap_15mmu=read.xlsx("heatmap_15mmu.xlsx",sheet = "Sheet1")
colnames(heatmap_15mmu)[1] <- "mass"
colnames(heatmap_15mmu) <- gsub("_cleaned","",colnames(heatmap_15mmu))

table=melt(heatmap_15mmu, id='mass')
colnames(table) <-c("mass","samples","intensity")
table$groups=names_correspondances_hudd_vs_bk$group[match(table$samples,names_correspondances_hudd_vs_bk$name)] #some differences between replicates
table$country=names_correspondances_hudd_vs_bk$country[match(table$samples,names_correspondances_hudd_vs_bk$name)]
table$samples=substr(table$samples,start=1,stop=2)

x_labels=c("A-BK1"="T-1","A-BK2"="T-2","B-BK"="T-3","A-HD1"="E-1","A-HD2"="E-2","B-HD"="E-3","B-HD3"="E-4","B-HD2"="E-5")

table2 <- table %>%
  group_by(samples,groups,country,mass) %>%
  dplyr::summarize(avg_intensity = mean(intensity, na.rm = TRUE), .groups = "drop") %>%
  mutate(avg_intensity = na_if(avg_intensity, 0),group2 = x_labels[groups])  %>%
  arrange(factor(country, levels=c("Thailand","England","Water Thailand, air England")),group2,samples) %>%
  mutate(replicate = substr(samples,start=2,stop=2)) %>%
  mutate(samples2 = paste(group2,replicate))

order_samples=unique(table2$samples2)

number_peaks <- table2 %>%
  group_by(samples2) %>%
  dplyr::summarise(count = sum(avg_intensity > 0, na.rm = TRUE))
  
add_table2=cbind.data.frame(samples=rep("Sel. features",nrow(selected_peaks)),groups=rep('none',nrow(selected_peaks)),country=rep("none",nrow(selected_peaks)),
                 mass=as.numeric(selected_peaks$mass_selected), avg_intensity=rep(50,nrow(selected_peaks)), group2=rep("none", nrow(selected_peaks)),
                 replicate=rep("",nrow(selected_peaks)),samples2=rep("Sel. features",nrow(selected_peaks)))

table2=rbind.data.frame(table2,add_table2)
order_samples=unique(table2$samples2)

table2$count=number_peaks$count[match(table2$samples2,number_peaks$samples2)]
table2$count[is.na(table2$count)] <- 188
table2$samples3=paste(table2$samples2, " (", table2$count, ")",sep="")

order_samples=unique(table2$samples3)

rep_colmap=ggplot(table2, aes(x = mass, y=avg_intensity)) +
  geom_col(aes(fill = country),width = 5) +
  scale_fill_manual(values=color_vector, name="Country of origin") +
  scale_y_continuous(breaks=seq(0,100,50)) +
  scale_x_continuous(limits=c(min(table2$mass[table2$avg_intensity > 0]-10),max(table2$mass[table2$avg_intensity > 0])+10), expand = c(0, 0))+
  facet_wrap(~factor(samples3, levels=order_samples), nrow=length(unique(table2$samples3)), strip.position = "right", axis.labels = "margins")+
  xlab("Mass ( m/z )")+
  ylab("Mean Intensity (%)")+
  theme_classic()+
  theme(legend.position = "none",strip.text.y = element_text(angle = 0),strip.background = element_blank())

rep_colmap

ggsave(paste(path_work,"/colmap.pdf", sep=""),rep_colmap, device = cairo_pdf,
        width = 30, height = 30, units = "cm")



#### load KDA predictive model and represent ####

kda_predictive=read.csv("kda_scores_2groups.csv",sep=",",header=FALSE)
kda_predictive=kda_predictive %>% separate_wider_delim(V1, delim = ",", names = c("Index", "Component1","Component2","Component3")) %>% mutate_if(is.character, as.numeric)
kda_predictive$Index[kda_predictive$Index %in% seq(53,61,1)] <- seq(62,70,1)

add7=read.xlsx("kda_group_7ordination.xlsx", sheet = 1)
add7=add7[,-2]

kda_predictive=rbind.data.frame(kda_predictive,add7)

kda_predictive=merge(names_correspondances_hudd_vs_bk,kda_predictive, by="Index") #
kda_predictive=kda_predictive %>% mutate(groups2=x_labels[group],replicate = substr(name,start=2,stop=2),samples3=paste(groups2,replicate))

rep=ggplot(kda_predictive, aes(Component1,Component2))+ 
  geom_point(size=3, aes(shape=condition, 
                         color=country)) +
  scale_color_manual(values=color_vector, name="Country of origin") +
  scale_shape_manual(values=shapes_vector, name="Manufacturing conditions")+
  xlab("Axis 1 (99.99%)") + 
  ylab("Axis 2 (< 0.01%)") +  
  geom_vline(xintercept = mean(kda_predictive[kda_predictive$country=="England",]$Component1), color="blue")+
  geom_vline(xintercept = mean(kda_predictive[kda_predictive$country=="Thailand",]$Component1), color="red") +
  theme_classic()+theme(legend.position = "none")


kda_loocv_scores=read.xlsx("kda_loocv_2groups_edited.xlsx",sheet=1)

kda_loocv_scores=kda_loocv_scores %>%
  rename(Index=File) %>%
  mutate(Index=gsub("_cleaned.txt","",Index)) %>%
  mutate(groups=names_correspondances_hudd_vs_bk$group[match(Index,names_correspondances_hudd_vs_bk$name)])  %>%
  mutate(groups2=x_labels[groups],replicate = substr(Index,start=2,stop=3),samples3=paste(groups2,replicate)) %>%
  select(samples3,Thailand_post,England_post) 

kda_loocv_scores2=melt(kda_loocv_scores)
kda_loocv_scores2$country=ifelse(grepl("T",kda_loocv_scores2$samples3), "Thailand","England")

assign_rep=ggplot(kda_loocv_scores2, aes(x=samples3,y=value,fill=variable))+
  geom_col(alpha=0.7)+
  scale_fill_manual(values=c("red","blue"))+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_continuous(limits=c(0,100), name="Membership prob. (%)")+
  theme_classic2() +
  theme(
    axis.line.x = element_blank(),   # Remove x-axis line
    axis.ticks.x = element_blank(),   # Remove x-axis ticks
    axis.line.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0,'lines'))+
  geom_segment(aes(x = -1.5, xend = -1.5, y = 0, yend = 100), color = "black")

ggsave("assign_prob.png", assign_rep, units = "cm",width = 10, height = 5,
  dpi = 300)
