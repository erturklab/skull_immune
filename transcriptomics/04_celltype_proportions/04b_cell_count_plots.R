library(rhdf5)
library(Matrix)
library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(yaml)


#set paths and dates
DATA_DIR <- '/lustre/groups/ml01/workspace/louis.kuemmerle/projects/A1/data2/'
WORKING_DIR <- '/home/icb/louis.kuemmerle/projects/a1/code/skull_immune/'
figure_path <- paste0(WORKING_DIR,'figures/cell_counts/')
dir.create(figure_path, recursive=TRUE)
#file_path <- paste0(DATA_DIR, 'table/')
file_path <- paste0(DATA_DIR, 'table_oct22/')
data_path <- DATA_DIR
#full_data_file <- 'cellxgene_april21_wSham'
file_names <- list.files(path = file_path, pattern = '*.csv')
today <- format(Sys.Date(), '%y%m%d')

#print(file_names)
#print(tools::file_path_sans_ext(file_names))

#font sizes
S_SIZE = 17
M_SIZE = 19
L_SIZE = 21

#figure size
WIDTH=6
HEIGTH=10

# Legend:
# When showing the legends the bar plot widths are different from plot to plot 
# We can not arrange them since we don't produce all plots in one round (otherwise it would be possible e.g. with
# cowplot). If you need the legend use
#     guide_legend(title=LEGEND_LABEL,keyheight=0.8,keywidth=0.8,default.unit="cm")
# instead of
#     guides(fill=FALSE)

# Labels
X_LABEL = 'Region'
Y_LABEL = 'Cell proportion'
LEGEND_LABEL = 'Cell types'

#########################################
# Stacked percentage plot full dataset  #
#########################################

for(data_file in file_names){
    
    #read cell composition file
    cell_counts <- read.csv(paste0(file_path, data_file),check.names=FALSE)
    #reorder columns (neutrophils last)
    #cell_counts <- cell_counts[,c(1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,9)]
    
    cell_tmp <- cell_counts[, -1]
    #cell_tmp.a <- cell_tmp[,!colnames(cell_tmp) %in% c('neurons') ]#, "naive.neutrophils",  "progenitor.cells") ]
    cell_tmp.a <- cell_tmp#[,!colnames(cell_tmp)]
    cell_tmp.m <- melt(cell_tmp.a, id.vars=c( 'condition', 'region')) 
    cell_tmp.m <- cell_tmp.m[cell_tmp.m$value>0,]
    
    #compute percentages
    tmp <- aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition), FUN=sum)
    tmp1 <- aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition, region=cell_tmp.m$region), FUN=sum) 
    tmp2 <-  aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition, variable=cell_tmp.m$variable), FUN=sum) 
    tmp2$perc <- 0
    tmp3 <- aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition, region=cell_tmp.m$region,
                                                variable=cell_tmp.m$variable), FUN=sum) 
    tmp3$perc <- 0
    
    for (condi in tmp$condition){
      idx <- tmp2$condition==condi 
      tmp2$perc[idx] = tmp2$x[idx]/ tmp$x[tmp$condition==condi] * 100 
      for (regi in tmp1$region){
        idx2 <- tmp3$region==regi & tmp3$condition ==condi
        tmp3$perc[idx2] = tmp3$x[idx2]/ tmp1$x[tmp1$region==regi & tmp1$condition==condi] * 100 
      }
    }
    
    #set columns to factor instead of char to set the plotted orders
    tmp3$condition <- factor(tmp3$condition,levels = c("Naive", "Sham", "MCAO"))
    # rename to MCAo
    levels(tmp3$condition) <- c("Naive", "Sham", "MCAo")    
    tmp3$region <- factor(tmp3$region,levels = c("Brain", "Meninges", "Bones", "Calvaria", "Vertebra", "Scapula", "Humerus", "Pelvis", "Femur"))
    #tmp3$region <- factor(tmp3$region,levels = c("Brain", "Meninges", "Bones", "Skull", "Vertebra", "Scapula", "Humerus", "Pelvis", "Femur")) 
    
    #set color scheme
    library(RColorBrewer)
    nb.cols <- length(unique(tmp3[["variable"]]))#13
    #mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols)
    #mycolors <- list('neutrophils Ltf+' = '#1f77b4','neutrophils Srgn+' = '#ff7f0e','neutrophils Srgn+/Socs3+' = '#2ca02c','neutrophils Srgn+/Il1r2-' = '#d62728','neutrophils Ltf+/Mki67+' = '#9467bd')
    color_df <- read.csv(file = paste0(file_path,'colors/',data_file))
    mycolors <- split(color_df$color,color_df$celltype)
    
    g1 <- ggplot(tmp3, aes(x=region, y=perc, fill=variable)) +
      scale_fill_manual(values = mycolors) +
      geom_bar(position='stack',stat='identity')  +  
      labs(x=X_LABEL, y=Y_LABEL) +
      facet_wrap(~condition) +
      theme_classic() + 
      guides(fill=FALSE)+#guide_legend(title=LEGEND_LABEL,keyheight=0.8,keywidth=0.8,default.unit="cm")) +
      theme(strip.text = element_text(face="bold", size=L_SIZE),   # size of subtitles ("MCAO", "Naive")
            axis.title.x = element_text(face="bold", size=M_SIZE), # x axis label (Cell type)
            axis.title.y = element_text(face="bold", size=M_SIZE), # y axis label (Cell proportion)
            axis.text.x  = element_text(angle=90, hjust = 1, vjust =0.5,  size=S_SIZE), # x axis ticks (regions)
            axis.text.y = element_text(size=S_SIZE), # y axis tick labels (cell proportion percentage)
            legend.title = element_text(size=M_SIZE),
            legend.text = element_text(size=S_SIZE))
    g1
    ggsave(filename = paste0(figure_path,today, '_stacked_cells_per_condition_', data_file, '.pdf'), 
          plot = g1, width = WIDTH, height=HEIGTH) 
    
    # Stacked percentage plot full dataset bones only #
    tmp4 <- tmp3[!tmp3$region %in% c('Brain','Meninges'),]
    
    g2 <- ggplot(tmp4, aes(x=region, y=perc, fill=variable)) +
      scale_fill_manual(values = mycolors) +
      geom_bar(position='stack',stat='identity')  +  
      labs(x=X_LABEL, y=Y_LABEL) +
      facet_wrap(~condition) +
      theme_classic() + 
      guides(fill=FALSE)+#guide_legend(title=LEGEND_LABEL,keyheight=0.8,keywidth=0.8,default.unit="cm")) +
      theme(strip.text = element_text(face="bold", size=L_SIZE),   # size of subtitles ("MCAO", "Naive")
            axis.title.x = element_text(face="bold", size=M_SIZE), # x axis label (Cell type)
            axis.title.y = element_text(face="bold", size=M_SIZE), # y axis label (Cell proportion)
            axis.text.x  = element_text(angle=90, hjust = 1, vjust =0.5,  size=S_SIZE), # x axis ticks (regions)
            axis.text.y = element_text(size=S_SIZE), # y axis tick labels (cell proportion percentage)
            legend.title = element_text(size=M_SIZE),
            legend.text = element_text(size=S_SIZE))
    g2
    ggsave(filename = paste0(figure_path,today, '_stacked_cells_per_condition_', data_file, '_bones.pdf'), 
           plot = g2, width = WIDTH, height=HEIGTH) 
}


##########################################
# Stacked percentage plots bones pooled  #
##########################################

for(data_file in file_names){
    
    #read cell composition file
    cell_counts <- read.csv(paste0(file_path, 'bones_pooled/', data_file),check.names=FALSE)
    
    cell_tmp <- cell_counts[, -1]
    cell_tmp.a <- cell_tmp
    cell_tmp.m <- melt(cell_tmp.a, id.vars=c( 'condition', 'region')) 
    cell_tmp.m <- cell_tmp.m[cell_tmp.m$value>0,]
    
    #compute percentages
    tmp <- aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition), FUN=sum)
    tmp1 <- aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition, region=cell_tmp.m$region), FUN=sum) 
    tmp2 <-  aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition, variable=cell_tmp.m$variable), FUN=sum) 
    tmp2$perc <- 0
    tmp3 <- aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition, region=cell_tmp.m$region,
                                                variable=cell_tmp.m$variable), FUN=sum) 
    tmp3$perc <- 0
    
    for (condi in tmp$condition){
      idx <- tmp2$condition==condi 
      tmp2$perc[idx] = tmp2$x[idx]/ tmp$x[tmp$condition==condi] * 100 
      for (regi in tmp1$region){
        idx2 <- tmp3$region==regi & tmp3$condition ==condi
        tmp3$perc[idx2] = tmp3$x[idx2]/ tmp1$x[tmp1$region==regi & tmp1$condition==condi] * 100 
      }
    }
    
    #set columns to factor instead of char to set the plotted orders
    tmp3$condition <- factor(tmp3$condition,levels = c("Naive", "Sham", "MCAO"))
    # rename to MCAo
    levels(tmp3$condition) <- c("Naive", "Sham", "MCAo")
    tmp3$region <- factor(tmp3$region,levels = c("Brain", "Meninges", "Bones", "Calvaria", "Vertebra", "Scapula", "Humerus", "Pelvis", "Femur"))
    #tmp3$region <- factor(tmp3$region,levels = c("Brain", "Meninges", "Bones", "Skull", "Vertebra", "Scapula", "Humerus", "Pelvis", "Femur"))
    
    #set color scheme
    library(RColorBrewer)
    nb.cols <- length(unique(tmp3[["variable"]]))
    color_df <- read.csv(file = paste0(file_path,'colors/',data_file))
    mycolors <- split(color_df$color,color_df$celltype)
    
    g1 <- ggplot(tmp3, aes(x=region, y=perc, fill=variable)) +
      scale_fill_manual(values = mycolors) +
      geom_bar(position='stack',stat='identity')  +  
      labs(x=X_LABEL, y=Y_LABEL) +
      facet_wrap(~condition) +
      theme_classic() + 
      guides(fill=FALSE)+#guide_legend(title=LEGEND_LABEL,keyheight=0.8,keywidth=0.8,default.unit="cm")) +
      theme(strip.text = element_text(face="bold", size=L_SIZE),   # size of subtitles ("MCAO", "Naive")
            axis.title.x = element_text(face="bold", size=M_SIZE), # x axis label (Cell type)
            axis.title.y = element_text(face="bold", size=M_SIZE), # y axis label (Cell proportion)
            axis.text.x  = element_text(angle=90, hjust = 1, vjust =0.5,  size=S_SIZE), # x axis ticks (regions)
            axis.text.y = element_text(size=S_SIZE), # y axis tick labels (cell proportion percentage)
            legend.title = element_text(size=M_SIZE),
            legend.text = element_text(size=S_SIZE))
    g1
    ggsave(filename = paste0(figure_path,today, '_stacked_cells_per_condition_', data_file, '_bones_pooled.pdf'), 
          plot = g1, width = WIDTH, height=HEIGTH) 
}


##########################################################################
# Stacked percentage plots with brain and meninges but only CD45+ cells  #
##########################################################################

file_names <- list.files(path = paste0(file_path, 'immune_only/'), pattern = '*.csv')

for(data_file in file_names){
    
    #read cell composition file
    cell_counts <- read.csv(paste0(file_path, 'immune_only/', data_file),check.names=FALSE)    
    #reorder columns (neutrophils last)
    #cell_counts <- cell_counts[,c(1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,9)]
    
    cell_tmp <- cell_counts[, -1]
    #cell_tmp.a <- cell_tmp[,!colnames(cell_tmp) %in% c('neurons') ]#, "naive.neutrophils",  "progenitor.cells") ]
    cell_tmp.a <- cell_tmp#[,!colnames(cell_tmp)]
    cell_tmp.m <- melt(cell_tmp.a, id.vars=c( 'condition', 'region')) 
    cell_tmp.m <- cell_tmp.m[cell_tmp.m$value>0,]
    
    #compute percentages
    tmp <- aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition), FUN=sum)
    tmp1 <- aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition, region=cell_tmp.m$region), FUN=sum) 
    tmp2 <-  aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition, variable=cell_tmp.m$variable), FUN=sum) 
    tmp2$perc <- 0
    tmp3 <- aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition, region=cell_tmp.m$region,
                                                variable=cell_tmp.m$variable), FUN=sum) 
    tmp3$perc <- 0
    
    for (condi in tmp$condition){
      idx <- tmp2$condition==condi 
      tmp2$perc[idx] = tmp2$x[idx]/ tmp$x[tmp$condition==condi] * 100 
      for (regi in tmp1$region){
        idx2 <- tmp3$region==regi & tmp3$condition ==condi
        tmp3$perc[idx2] = tmp3$x[idx2]/ tmp1$x[tmp1$region==regi & tmp1$condition==condi] * 100 
      }
    }
    
    #set columns to factor instead of char to set the plotted orders
    tmp3$condition <- factor(tmp3$condition,levels = c("Naive", "Sham", "MCAO"))
    # rename to MCAo
    levels(tmp3$condition) <- c("Naive", "Sham", "MCAo")    
    tmp3$region <- factor(tmp3$region,levels = c("Brain", "Meninges", "Bones", "Calvaria", "Vertebra", "Scapula", "Humerus", "Pelvis", "Femur"))
    #tmp3$region <- factor(tmp3$region,levels = c("Brain", "Meninges", "Bones", "Skull", "Vertebra", "Scapula", "Humerus", "Pelvis", "Femur")) 
    
    #set color scheme
    library(RColorBrewer)
    nb.cols <- length(unique(tmp3[["variable"]]))#13
    #mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols)
    #mycolors <- list('neutrophils Ltf+' = '#1f77b4','neutrophils Srgn+' = '#ff7f0e','neutrophils Srgn+/Socs3+' = '#2ca02c','neutrophils Srgn+/Il1r2-' = '#d62728','neutrophils Ltf+/Mki67+' = '#9467bd')
    color_df <- read.csv(file = paste0(file_path,'colors/',data_file))
    mycolors <- split(color_df$color,color_df$celltype)
    
    g1 <- ggplot(tmp3, aes(x=region, y=perc, fill=variable)) +
      scale_fill_manual(values = mycolors) +
      geom_bar(position='stack',stat='identity')  +  
      labs(x=X_LABEL, y=Y_LABEL) +
      facet_wrap(~condition) +
      theme_classic() + 
      guides(fill=FALSE)+#guide_legend(title=LEGEND_LABEL,keyheight=0.8,keywidth=0.8,default.unit="cm")) +
      theme(strip.text = element_text(face="bold", size=L_SIZE),   # size of subtitles ("MCAO", "Naive")
            axis.title.x = element_text(face="bold", size=M_SIZE), # x axis label (Cell type)
            axis.title.y = element_text(face="bold", size=M_SIZE), # y axis label (Cell proportion)
            axis.text.x  = element_text(angle=90, hjust = 1, vjust =0.5,  size=S_SIZE), # x axis ticks (regions)
            axis.text.y = element_text(size=S_SIZE), # y axis tick labels (cell proportion percentage)
            legend.title = element_text(size=M_SIZE),
            legend.text = element_text(size=S_SIZE))
    g1
    ggsave(filename = paste0(figure_path,today, '_stacked_cells_per_condition_', data_file, '_immune_only.pdf'), 
          plot = g1, width = WIDTH, height=HEIGTH)
}




#################################################
## Stacked percentage plot neutrophils dataset  #
#################################################
#data_file <- neutro_data_file
#
##read cell composition file
#cell_counts <- read.csv(paste0(file_path, 'cell_counts_', data_file, '.csv'),check.names=FALSE)
##reorder columns
#cell_counts <- cell_counts[,c(1,2,3,6,5,7,4,8)] 
#
#cell_tmp <- cell_counts[, -1]
#cell_tmp.a <- cell_tmp[,!colnames(cell_tmp) %in% c('neurons') ]#, "naive.neutrophils",  "progenitor.cells") ]
#cell_tmp.m <- melt(cell_tmp.a, id.vars=c( 'condition', 'region')) 
#cell_tmp.m <- cell_tmp.m[cell_tmp.m$value>0,]
#
##compute percentages
#tmp <- aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition), FUN=sum)
#tmp1 <- aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition, region=cell_tmp.m$region), FUN=sum) 
#tmp2 <-  aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition, variable=cell_tmp.m$variable), FUN=sum) 
#tmp2$perc <- 0
#tmp3 <- aggregate(cell_tmp.m$value, by=list(condition=cell_tmp.m$condition, region=cell_tmp.m$region,
#                                            variable=cell_tmp.m$variable), FUN=sum) 
#tmp3$perc <- 0
#
#for (condi in tmp$condition){
#  idx <- tmp2$condition==condi 
#  tmp2$perc[idx] = tmp2$x[idx]/ tmp$x[tmp$condition==condi] * 100 
#  for (regi in tmp1$region){
#    idx2 <- tmp3$region==regi & tmp3$condition ==condi
#    tmp3$perc[idx2] = tmp3$x[idx2]/ tmp1$x[tmp1$region==regi & tmp1$condition==condi] * 100 
#  }
#}
#
##set color scheme
#library(RColorBrewer)
#nb.cols <- length(unique(tmp3[["variable"]]))#13
##mycolors <- colorRampPalette(brewer.pal(8, "RdYlBu"))(nb.cols)
##mycolors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd')
#mycolors <- list('neutrophils Ltf+' = '#1f77b4','neutrophils Srgn+' = '#ff7f0e','neutrophils Srgn+/Socs3+' = '#2ca02c','neutrophils Srgn+/Il1r2-' = '#d62728','neutrophils Ltf+/Mki67+' = '#9467bd')
#
#g1 <- ggplot(tmp3, aes(x=region, y=perc, fill=variable)) +
#  scale_fill_manual(values = mycolors) +
#  geom_bar(position='stack',stat='identity')  +  
#  labs(x=X_LABEL, y=Y_LABEL) +
#  facet_wrap(~condition) +
#  theme_classic() + 
#  guides(fill=guide_legend(title=LEGEND_LABEL,keyheight=0.8,keywidth=0.8,default.unit="cm")) +
#  theme(strip.text = element_text(face="bold", size=L_SIZE),   # size of subtitles ("MCAO", "Naive")
#        axis.title.x = element_text(face="bold", size=M_SIZE), # x axis label (Cell type)
#        axis.title.y = element_text(face="bold", size=M_SIZE), # y axis label (Cell proportion)
#        axis.text.x  = element_text(angle=90, hjust = 1, vjust =0.5,  size=S_SIZE), # x axis ticks (regions)
#        axis.text.y = element_text(size=S_SIZE), # y axis tick labels (cell proportion percentage)
#        legend.title = element_text(size=M_SIZE),
#        legend.text = element_text(size=S_SIZE))
#g1
#ggsave(filename = paste0(figure_path,today, '_stacked_cells_per_condition_', data_file, '.pdf'), 
#       plot = g1, width = 10, height=10) 
#
## Stacked percentage plot full dataset bones only #
#tmp4 <- tmp3[!tmp3$region %in% c('Brain','Meninges'),]
#
#g2 <- ggplot(tmp4, aes(x=region, y=perc, fill=variable)) +
#  scale_fill_manual(values = mycolors) +
#  geom_bar(position='stack',stat='identity')  +  
#  labs(x=X_LABEL, y=Y_LABEL) +
#  facet_wrap(~condition) +
#  theme_classic() + 
#  guides(fill=guide_legend(title=LEGEND_LABEL,keyheight=0.8,keywidth=0.8,default.unit="cm")) +
#  theme(strip.text = element_text(face="bold", size=L_SIZE),   # size of subtitles ("MCAO", "Naive")
#        axis.title.x = element_text(face="bold", size=M_SIZE), # x axis label (Cell type)
#        axis.title.y = element_text(face="bold", size=M_SIZE), # y axis label (Cell proportion)
#        axis.text.x  = element_text(angle=90, hjust = 1, vjust =0.5,  size=S_SIZE), # x axis ticks (regions)
#        axis.text.y = element_text(size=S_SIZE), # y axis tick labels (cell proportion percentage)
#        legend.title = element_text(size=M_SIZE),
#        legend.text = element_text(size=S_SIZE))
#g2
#ggsave(filename = paste0(figure_path,today, '_stacked_cells_per_condition_', data_file, '_bones.pdf'), 
#       plot = g2, width = 10, height=10)
