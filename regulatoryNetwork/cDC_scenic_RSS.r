library('dplyr')
library("tidyverse")
library(SCENIC)
library('ggplot2')
library('RColorBrewer')
library(AUCell)
library(ggrepel)


################## read in all scenic results
AUClist <- list()

for (celltype in c('cDC1', 'cDC2', 'cDC3')){
  setwd("/net/bmc-lab5/data/kellis/users/ruiwenfu/SCENIC")
  setwd(celltype)

  scenicOptions <- readRDS("int/scenicOptions.Rds")

  cellInfo <- readRDS('int/cellInfo.Rds')
  regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
  regulonAUC <- getAUC(regulonAUC) %>% as.data.frame %>% rownames_to_column('row')

  AUClist[[celltype]] <- regulonAUC
}

AUCmerged <- AUClist %>% reduce(full_join, by = 'row') %>%
    column_to_rownames('row') %>%
    mutate_all(~replace(., is.na(.), 0))

################### calculate RSS
annotation <- c()
for (celltype in c('cDC1', 'cDC2', 'cDC3')){
  annotation <- c(annotation, rep(c(celltype),each=ncol(AUClist[[celltype]])))
}

RSS <- calcRSS(AUCmerged, annotation)

RSS <- RSS %>% as.data.frame %>% rownames_to_column('row')

write.csv(RSS ,
        'cDC_RSS_score.csv')


################## plotting
for (celltype in c('cDC1', 'cDC2', 'cDC3')){
  setwd("/net/bmc-lab5/data/kellis/users/ruiwenfu/SCENIC")
  setwd(celltype)

  RSS$Regulon_Specificity_Score <-RSS[[celltype]]
  RSS <-RSS  %>%
    mutate(rank = rank(-RSS$Regulon_Specificity_Score),
               label = ifelse(rank <=10, row, ''),
                color = ifelse(rank <=10, 'red', ''))

p <- ggplot(RSS, aes(x = rank, y = Regulon_Specificity_Score, color = color, label = label))+
    geom_point()+
    scale_color_manual(values=c('black',"red"))+
    theme_classic()  +
        theme(legend.position = "none")+
    theme(text = element_text(size = 8)) +
     geom_text_repel(    force = 40,
    force_pull   = 0, # do not pull toward data points
    nudge_y      = -0.5,
    direction    = "x",
    angle        = 90,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1,
                     size = 2
      )

ggsave('RSS_plot.png',plot = p, width = 4, height = 4)

}
