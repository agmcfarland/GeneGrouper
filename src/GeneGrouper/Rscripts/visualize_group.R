
packages <- c("reshape", "ggplot2", "cowplot", "dplyr", "gggenes", "groupdata2", "svglite")
install.packages(setdiff(packages, rownames(installed.packages()))) 

library(reshape)
library(ggplot2)
library(cowplot)
library(dplyr)
library(gggenes)
library(groupdata2)

packageVersion("reshape")
packageVersion("ggplot2")
packageVersion("cowplot")
packageVersion("dplyr")
packageVersion("gggenes")
packageVersion("groupdata2")

args <-  commandArgs(trailingOnly = TRUE)
results_dir <- args[1]
visualizations_dir <- args[2]
image_format <- args[3]

setwd(results_dir)

## Troubleshooting ##
#setwd('/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset1/test1/mexb/results') #debugging mac
#visualizations_dir <- '/Users/owlex/Dropbox/Documents/Northwestern/Hartmann_Lab/syntenease_project/gtr/testbed/dataset1/test1/mexb/visualizations' #debugging mac
#setwd('/projects/b1042/HartmannLab/alex/GeneGrouper_test/testbed/dataset1/test1/pdua/results') #debugging linux
#visualizations_dir <- '/projects/b1042/HartmannLab/alex/GeneGrouper_test/testbed/dataset1/test1/pdua/visualizations' #debugging linux
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
df_regions1 <- read.csv('rtable_cwc_regions.csv', na.strings=c('','NA'))

# add text overlays
df_regions1 <- df_regions1%>%
  mutate(refseq_gene=coalesce(refseq_gene,refseq_product))%>%
  mutate(refseq_gene=coalesce(refseq_gene, as.character(ortho_cluster_id)))%>%
  arrange(order)%>%
  mutate(mod_ortho_cluster_id=ortho_cluster_id)%>%
  mutate(text_highlight=case_when(pseudo_check=='p' ~ 'X',
                                  TRUE~NA_character_))

# standardize region order by factor labeling
df_regions1$cwc_id_dummy <- factor(df_regions1$cwc_id,levels=unique(as.character(df_regions1$cwc_id)))
df_regions1$middle_label <- (df_regions1$norm_start+df_regions1$norm_end-5)/2


## Chunk up output
# make chunks
region_chunk_size <- 30
# check if the number of regions is smaller than the defined chunk size
if (length(levels(df_regions1$cwc_id_dummy )) <= region_chunk_size){
  region_chunk_size = length(levels(df_regions1$cwc_id_dummy ))
}

df_chunk_order <- data.frame('cwc_id'=levels(df_regions1$cwc_id_dummy ), 
                             'region_order' = seq(1,length(levels(df_regions1$cwc_id_dummy ))))
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

#### Begin section for visualization ####

max_count <- read.csv('rtable_cwc_dis_counts.csv')%>%select(count)%>%filter(count==max(count))%>%pull()
dbscan_label_inspected <- unique(df_regions1$dbscan_label)

for (ch in unique(df_chunk_order$chunk_group)){
  #ch <- 1
  label_keep <- as.vector(as.integer(df_chunk_order%>%filter(chunk_group == ch)%>%select(cwc_id)%>%pull()))
  
  
  p1_cwc_ggregions <- ggplot(df_regions1%>%filter(cwc_id%in%label_keep),
                         aes( xmin=norm_start, xmax=norm_end, y = cwc_id_dummy, fill=as.factor(ortho_cluster_id), label = as.character(refseq_gene), forward=strand))+
    facet_wrap(~ cwc_id_dummy,scales='free_y',ncol=1)+
    #geom_vline(xintercept = 0,size=.5,color='grey')+
    geom_gene_arrow(arrow_body_height = unit(v_abh,'cm'),arrowhead_height=unit(v_ah,'cm'), arrowhead_width=unit(v_aw,'cm'))+
    geom_gene_label(align='centre',grow=T, reflow=T, min.size=.5,size=6)+
    geom_text(aes(x=middle_label,y=cwc_id_dummy,label=as.character(mod_ortho_cluster_id)), vjust=v_gtext_vjust, size=v_gtext_size)+
    geom_text(aes(x=middle_label,y=cwc_id_dummy,label=as.character(text_highlight)), vjust=v_gtext2_vjust, size=v_gtext2_size,alpha=0.5)+
    theme_genes()+
    theme(legend.position='none',
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=15, face='bold'),
          strip.background = element_rect(fill = "transparent", colour = NA),
          panel.background =  element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "white"),#, colour = NA),#element_blank()
    )+
    xlab(label='Normalized gene position (bp)')+
    scale_y_discrete(
      labels = df_chunk_order%>%filter(chunk_group==ch)%>%mutate(cwc_id=paste('s',cwc_id,sep=''))%>%select(cwc_id)%>%pull(),
      breaks = df_chunk_order%>%filter(chunk_group==ch)%>%select(cwc_id)%>%pull())
  
  df_cwc_ggdc <- read.csv('rtable_cwc_dis_counts.csv')%>%
    arrange(representative_relative_dissimilarity)
  df_cwc_ggdc$cwc_id_dummy <- factor(df_cwc_ggdc$cwc_id, levels = levels(df_regions1$cwc_id_dummy))  
  
  
  p2_cwc_ggdc <- ggplot(df_cwc_ggdc%>%filter(cwc_id%in%label_keep),
                        aes(x=count,y=cwc_id_dummy))+
    geom_bar(stat='identity',position='dodge')+
    theme_classic()+
    facet_wrap(~cwc_id_dummy,ncol=1,scales='free_y')+
    theme(legend.position='none',
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=15, color = 'grey50'),
          axis.line.y=element_blank(),
          panel.grid.major.y=element_line(size=1.1),
          panel.grid.major.x=element_line(size=0.75,linetype='dashed'),
          axis.ticks.y=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank())+
    ylab(label='CwC label')+
    xlab(label='Count')+
    scale_x_continuous(limits=c(0,max_count))+
    scale_y_discrete(
      labels = df_chunk_order%>%filter(chunk_group==ch)%>%mutate(region_order=paste(region_order,'.',sep=''))%>%select(region_order)%>%pull(),
      breaks = df_chunk_order%>%filter(chunk_group==ch)%>%select(cwc_id)%>%pull(),
      limits = rev)
  
  
  # relative dissimilarity of region sub clusters
  p3_cwc_ggdc <- ggplot(df_cwc_ggdc%>%filter(cwc_id%in%label_keep)
                        ,aes(x=representative_relative_dissimilarity,y=cwc_id_dummy))+
    geom_bar(stat='identity',position='dodge')+
    facet_wrap(~cwc_id_dummy,ncol=1,scales='free_y')+
    theme_classic()+
    theme(legend.position='none',
          axis.title.y=element_blank(),
          #axis.text.y=element_text(size=15),
          axis.text.y=element_blank(),
          axis.line.y=element_blank(),
          panel.grid.major.y=element_line(size=1.1),
          panel.grid.major.x=element_line(size=0.75,linetype='dashed'),
          axis.ticks.y=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank())+
    scale_x_continuous(limits=c(0,1.05))+
    ylab(label='CwC label')+
    xlab(label='Relative dissimilarity to sub cluster 0')
  
  p1mgg <- plot_grid(p2_cwc_ggdc,p1_cwc_ggregions,p3_cwc_ggdc,align='hv',nrow=1,rel_widths = c(1/4,1/2,1/4))
  

  save_plot(paste(visualizations_dir,'/inspect_group_',dbscan_label_inspected,'_',ch,'.',image_format,sep=''),p1mgg,base_height=8,base_aspect_ratio = 2)
  
  }



