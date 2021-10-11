

packages <- c("reshape", "ggplot2", "cowplot", "dplyr", "gggenes", "groupdata2", "svglite")
install.packages(setdiff(packages, rownames(installed.packages()))) 

library(reshape)
library(ggplot2)
library(cowplot)
library(dplyr)
library(gggenes)
library(groupdata2)
library(svglite)

args <-  commandArgs(trailingOnly = TRUE)
results_dir <- args[1]
visualizations_dir <- args[2]
image_format <- args[3]

setwd(results_dir)

## Troubleshooting ##

#setwd('/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset2/test2/pdua/results') #debugging
#visualizations_dir <- '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset2/test2/pdua/visualizations' #debugging

#setwd('/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset4/test2/mexb/internal_data') #debugging
#visualizations_dir <- '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset4/test2/mexb/visualizations' #debugging
#image_format <- 'svg'

## ## ## ## ## ## ## ##

# dataframe with scaling parameters for each
df_scaling <- data.frame('unique_region_count' = c(10, 20, 30, 40 ),
                         'v_abh' = c(.8, .4, .27, .2 ), #arrow body height
                         'v_ah' = c(1, .5, .4, .25), # arrowhead height
                         'v_aw' = c(.5, .25, .125, .125), # arrowhead width
                         'v_gtext_vjust' = c(-3.2, -3.1, -2.2, -3.2), # orthocluster ID vjust
                         'v_gtext_size' = c(2, 1, 1, .5), # orthocluster ID size
                         'v_gtext2_vjust' = c(.3, .3, .3, .6), # pseudo gene highlight  vjust
                         'v_gtext2_size' = c(4, 4, 4, 4) ## pseudo gene highlight  size
                         )


#### Main visualization prepare ####
df_regions1 <- read.csv('rtable_region_representatives.csv', na.strings=c('','NA'))

# add text overlays
df_regions1 <- df_regions1%>%
  mutate(refseq_gene=coalesce(refseq_gene,refseq_product))%>%
  mutate(refseq_gene=coalesce(refseq_gene, as.character(ortho_cluster_id)))%>%
  arrange(order)%>%
  mutate(mod_ortho_cluster_id=ortho_cluster_id)%>%
  mutate(text_highlight=case_when(pseudo_check=='p' ~ 'X',
                                  TRUE~NA_character_))

# standardize region order by factor labeling
clustered_list <- df_regions1%>%select(dbscan_label)%>%unique()%>%filter(dbscan_label>-1)

# making factor levels of the dbscan labels. making unclusteed label '-1' be at the end.
df_regions1$dbscan_label_dummy <- factor(df_regions1$dbscan_label,levels=rev(unique(as.character(df_regions1$dbscan_label))))
if ("-1" %in% levels(df_regions1$dbscan_label_dummy) ){
  relevel <- setdiff(levels(df_regions1$dbscan_label_dummy),'-1')
  relevel <- c(relevel,'-1')
  df_regions1$dbscan_label_dummy <- factor(df_regions1$dbscan_label,levels=relevel)
}


df_regions1$middle_label <- (df_regions1$norm_start+df_regions1$norm_end-5)/2


## Chunk up output
# make chunks
region_chunk_size <- 30 #default chunk size if the if statement below is not met
# check if the number of regions is smaller than the defined chunk size
if (length(levels(df_regions1$dbscan_label_dummy )) <= region_chunk_size){
  region_chunk_size = length(levels(df_regions1$dbscan_label_dummy ))
}

df_chunk_order <- data.frame('dbscan_label'= (levels(df_regions1$dbscan_label_dummy )), 
                             'region_order' = seq(1,length(levels(df_regions1$dbscan_label_dummy ))))
df_chunk_order <- groupdata2::group(df_chunk_order, n = region_chunk_size, method = 'greedy')
df_chunk_order <- df_chunk_order%>%rename(chunk_group=3)


## Scaling plotting inputs
if (region_chunk_size <= 10){
  df_scaling <- df_scaling%>%filter(unique_region_count == 10)
} else if (region_chunk_size <= 20) {
  df_scaling <- df_scaling%>%filter(unique_region_count == 20)
} else if (region_chunk_size <= 30) {
  df_scaling <- df_scaling%>%filter(unique_region_count == 30)
} else if (region_chunk_size >= 31) {
  df_scaling <- df_scaling%>%filter(unique_region_count == 40)
} 

v_abh <- df_scaling%>%select(v_abh)%>%pull()
v_ah <- df_scaling%>%select(v_ah)%>%pull()
v_aw <- df_scaling%>%select(v_aw)%>%pull()
v_gtext_vjust <- df_scaling%>%select(v_gtext_vjust)%>%pull()
v_gtext_size <- df_scaling%>%select(v_gtext_size)%>%pull()
v_gtext2_vjust <- df_scaling%>%select(v_gtext2_vjust)%>%pull()
v_gtext2_size <- df_scaling%>%select(v_gtext2_size)%>%pull()


#### Plot Main visualization ####

for (ch in unique(df_chunk_order$chunk_group)){
  #ch = 1
  label_keep <- as.vector(as.integer(df_chunk_order%>%filter(chunk_group == ch)%>%select(dbscan_label)%>%pull()))
  
  p1_ggregions <- ggplot(df_regions1%>%filter(dbscan_label%in%label_keep),
                         aes( xmin=norm_start, xmax=norm_end, y = dbscan_label_dummy, fill=as.factor(ortho_cluster_id), label = as.character(refseq_gene), forward=strand))+
    facet_wrap(~ dbscan_label_dummy,scales='free_y',ncol=1)+
    #geom_vline(xintercept = 0,size=.5,color='grey')+
    geom_gene_arrow(arrow_body_height = unit(v_abh,'cm'),arrowhead_height=unit(v_ah,'cm'), arrowhead_width=unit(v_aw,'cm'))+
    geom_gene_label(align='centre',grow=T, reflow=T, min.size=.5,size=6)+
    geom_text(aes(x=middle_label,y=dbscan_label_dummy,label=as.character(mod_ortho_cluster_id)), vjust=v_gtext_vjust, size=v_gtext_size)+
    geom_text(aes(x=middle_label,y=dbscan_label_dummy,label=as.character(text_highlight)), vjust=v_gtext2_vjust, size=v_gtext2_size,alpha=0.5)+
    theme_genes()+
    theme(legend.position='none',
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=15,face='bold'),
          #axis.text.y=element_blank(),
          strip.background = element_rect(fill = "transparent", colour = NA),
          panel.background =  element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "white"),#, colour = NA),#element_blank()
          )+
    scale_y_discrete(#labels = df_chunk_order%>%filter(chunk_group==ch)%>%mutate(region_order=paste(region_order,'. c',dbscan_label,sep=''))%>%select(region_order)%>%pull(),
      labels = df_chunk_order%>%filter(chunk_group==ch)%>%mutate(region_order=paste('g',dbscan_label,sep=''))%>%select(region_order)%>%pull(),
      breaks = df_chunk_order%>%filter(chunk_group==ch)%>%select(dbscan_label)%>%pull(),
      limits = rev)+
    xlab(label='Normalized gene position (bp)')
  
  
  # dissimilarity plot
  df_crm <- read.csv('rtable_label_dissimilarity.csv')%>%arrange(order)
  df_crm$dbscan_label_dummy <- factor(df_crm$dbscan_label,levels = levels(df_regions1$dbscan_label_dummy))
  
  
  p2_ggdis <- ggplot(df_crm%>%filter(dbscan_label%in%label_keep),
                     aes(x=dissimilarity,y=dbscan_label_dummy))+
    geom_boxplot(outlier.shape=NA,alpha=0.5,position='dodge',aes(fill=measurement))+
    geom_point(position=position_jitterdodge(),aes(color=measurement))+
    
    geom_text(data=df_crm%>%filter(dbscan_label%in%label_keep)%>%
                select(dbscan_label,dbscan_label_dummy,full_index_pos)%>%unique()%>%group_by(dbscan_label)%>%
                mutate(count=n())%>%ungroup()%>%select(dbscan_label_dummy,count)%>%unique()%>%mutate(x_pos=1.05),
              aes(x=x_pos, y=dbscan_label_dummy, label=as.character(count)), colour='grey34')+
    
    facet_wrap(~dbscan_label_dummy,ncol=1,scales='free_y')+
    
    theme_classic()+
    theme(legend.position='none',
          axis.title.y=element_blank(),
          #axis.text.y=element_text(size=15, color = 'transparent'),
          axis.text.y=element_blank(),
          axis.line.y=element_blank(),
          panel.grid.major.y=element_line(size=1.1),
          panel.grid.major.x=element_line(size=0.75,linetype='dashed'),
          axis.ticks.y=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank())+
    scale_x_continuous(limits=c(0,1.05))+
    ylab(label='Cluster label')+
    xlab(label='Dissimilarity')
  
  
  # identity and coverage plot
  df_id <- read.csv('rtable_label_identitycoverage.csv')%>%arrange(order)
  df_id$dbscan_label_dummy <- factor(df_id$dbscan_label, levels = levels(df_regions1$dbscan_label_dummy))
  
  p3_ggid <- ggplot(df_id%>%filter(dbscan_label%in%label_keep),
                    aes(x=value,y=dbscan_label_dummy))+
    geom_boxplot(outlier.shape=NA,alpha=0.5,position='dodge',aes(fill=variable))+
    geom_point(position=position_jitterdodge(),aes(color=variable, shape=representative_region))+
    facet_wrap(~dbscan_label_dummy, ncol=1, scales='free_y')+
    theme_classic()+
    theme(legend.position='none',
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=15, hjust = 0, color = 'grey50'),
          #axis.text.y=element_blank(),
          axis.line.y=element_blank(),
          panel.grid.major.y=element_line(size=1.1),
          panel.grid.major.x=element_line(size=0.75,linetype='dashed'),
          axis.ticks.y=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank())+
    scale_x_continuous(limits=c(0,105))+
    scale_y_discrete(
                     labels = df_chunk_order%>%filter(chunk_group==ch)%>%mutate(region_order=paste(region_order,'.',sep=''))%>%select(region_order)%>%pull(),
                     breaks = df_chunk_order%>%filter(chunk_group==ch)%>%select(dbscan_label)%>%pull(),
                     limits = rev)+
    ylab(label='Cluster label')+
    xlab(label='% value relative to seed gene')
  
  
  p7 <- plot_grid(p3_ggid,p1_ggregions,p2_ggdis,align='hv',nrow=1,rel_widths = c(0.25,0.5,0.25), axis='bt')
  save_plot(paste(visualizations_dir,'/group_summary_', ch,'.',image_format,sep='') ,p7,base_height=8,base_aspect_ratio = 2)
}


#### Plot taxa with regions found vs searched ####
# show breakdown of genomes with at least one region found by taxa
df_meta_orig <- read.csv('rtable_genome_input_metadata.csv')
df_s <- read.csv('rtable_summarized_region_clusters.csv')
df_meta <- merge(df_meta_orig, df_s%>%select(assembly_id)%>%unique()%>%mutate(found=1), all=TRUE, by='assembly_id')
df_meta <- df_meta%>%mutate(found=case_when(is.na(found)==TRUE~0,TRUE~found))%>%
  group_by(genus)%>%mutate(count=sum(found),total=n())%>%select(genus,count,total)%>%unique()
df_meta <- melt(as.data.frame(df_meta),id.vars='genus')
df_meta <- df_meta%>%arrange(desc(value))
df_meta$genus_dummy <- factor(df_meta$genus, levels=(unique(as.character(df_meta$genus))))

p1_taxafound <- ggplot(df_meta, aes(x=genus_dummy,y=value,fill=variable))+
  geom_bar(stat='identity',position='dodge')+
  geom_text(data=df_meta,aes(x=genus_dummy,y=value+5,label=value,group=variable),size=5,position=position_dodge(width=0.75))+
  theme_classic()+
  theme(legend.position='none',
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        axis.text.y=element_text(size=15),
        axis.line.y=element_blank(),
        panel.grid.major.y=element_line(size=1.1),
        panel.grid.major.x=element_line(size=0.75,linetype='dashed'),
        axis.ticks.y=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  ylab(label='Count')
save_plot(paste(visualizations_dir,'/taxa_searched','.',image_format,sep=''),p1_taxafound,base_height=8,base_aspect_ratio = 2)


#### Plot percentage of hits per cluster per taxa ####
# show breakdown by taxa and cluster label, using whole counts
df_s <- read.csv('rtable_summarized_region_clusters.csv')
df_s <- merge(df_s,df_meta_orig,by='assembly_id',all.x=TRUE)

for (ch in unique(df_chunk_order$chunk_group)){
  #ch = 1
  label_keep <- as.vector(as.integer(df_chunk_order%>%filter(chunk_group == ch)%>%select(dbscan_label)%>%pull()))
  
  # show breakdown by taxa and cluster label, shown as a percentage of regions per total number of genomes
  df_s2 <- merge(read.csv('rtable_summarized_region_clusters.csv'),read.csv('rtable_genome_input_metadata.csv'),by='assembly_id',all.x=TRUE)
  df_s2 <- df_s2%>%select(assembly_id,dbscan_label,genus)%>%unique()%>%group_by(genus,dbscan_label)%>%mutate(count=n())%>%ungroup()%>%select(genus,dbscan_label,count)%>%unique()
  df_gtotal <- read.csv('rtable_genome_input_metadata.csv')%>%select(assembly_id,genus)%>%group_by(genus)%>%summarize(total=n())
  df_s2 <- merge(df_s2,df_gtotal,by='genus',all.x=TRUE)
  df_s2 <- df_s2%>%mutate(percentage=round(count/total,4)*100)%>%select(genus,dbscan_label,percentage)%>%unique()
  df_s2 <- cast(df_s2,genus~dbscan_label)
  df_s2 <- df_s2%>%replace(is.na(.),0)
  df_s2 <- melt(df_s2,id.vars='genus')
  df_mr <- merge(read.csv('rtable_summarized_region_clusters.csv'),read.csv('rtable_genome_input_metadata.csv'),by='assembly_id',all.x=TRUE)
  df_mr <- df_mr%>%select(assembly_id,region_id,dbscan_label,genus)%>%unique()%>%
  group_by(genus,dbscan_label)%>%mutate(multiregion=case_when(length(unique(assembly_id))!=length(region_id)~'mr',TRUE~'sr'))%>%ungroup()%>%
  select(genus,dbscan_label,multiregion)%>%unique()
  df_s2 <- merge(df_s2,df_mr,by=c('genus','dbscan_label'),all.x=TRUE)
  
  df_s2$dbscan_label_dummy <- factor(df_s2$dbscan_label,levels = rev(levels(df_regions1$dbscan_label_dummy)))
  df_s2$genus_dummy <- factor(df_s2$genus,levels = levels(df_meta$genus_dummy))
  
  # start
  p3_region_taxa_breakdown <- ggplot(df_s2%>%filter(dbscan_label%in%label_keep),
                                     aes(x=genus_dummy,y=dbscan_label_dummy))+
  geom_tile(aes(height=0.95,width=0.95,fill=value),color='black',alpha=0.85)+
  geom_label(data=df_s2%>%filter(dbscan_label%in%label_keep)%>%
               mutate(value2=case_when(value==0~NA_character_,
                                                  value>0&multiregion=='mr'~paste(as.character(value),'*',sep=' '),
                                                  value>0&multiregion=='sr'~as.character(value))),
             aes(x=genus_dummy,y=dbscan_label_dummy,label=as.character(value2)),size=4,color='black',fill='white')+
  theme_classic()+
  theme(legend.position='none',
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=15, color = 'grey50'),
        axis.line.y=element_blank(),
        panel.grid.major.y=element_line(size=1.1),
        panel.grid.major.x=element_line(size=0.75,linetype='dashed'),
        axis.ticks.y=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())+
  scale_fill_gradient(low='white',high='darkblue',limits=c(0,100))+
  scale_y_discrete(
    labels = df_chunk_order%>%filter(chunk_group==ch)%>%mutate(region_order=paste(region_order,'.',sep=''))%>%select(region_order)%>%pull(),
    breaks = df_chunk_order%>%filter(chunk_group==ch)%>%select(dbscan_label)%>%pull())
  
  
  p4_empty <-  ggplot(df_s2%>%filter(dbscan_label%in%label_keep)%>%select(dbscan_label,dbscan_label_dummy)%>%unique(),
                     aes(y=dbscan_label_dummy))+
    geom_blank()+
    theme_classic()+
    theme(legend.position='none',
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=15, face='bold'),
          axis.line=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_blank(),
          axis.ticks.y=element_blank())+
    scale_y_discrete(
      labels = df_chunk_order%>%filter(chunk_group==ch)%>%mutate(dbscan_label=paste('g',dbscan_label,sep=''))%>%select(dbscan_label)%>%pull(),
      breaks = df_chunk_order%>%filter(chunk_group==ch)%>%select(dbscan_label)%>%pull())
    
  
  p8 <- plot_grid(p3_region_taxa_breakdown, p4_empty, align='hv',nrow=1,rel_widths = c(0.98,0.05), axis='bt')
  save_plot(paste(visualizations_dir,'/groups_by_taxa_',ch,'.',image_format,sep=''),p8,base_height=8,base_aspect_ratio = 2)
  
}



  

