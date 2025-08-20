setwd("~/projects/NPC_SpatialProteomics/")
library(tidyverse)
library(data.table)
source('0.Helpers.R')

MT_CT <- qs::qread('SingleCellDataFrame.qs')

LEVELS <- c('CT0', 'CT1', 'CT2', 'CT3', 'CT4') %>% 
  purrr::set_names() %>% 
  purrr::map(\(X) levels(MT_CT[[X]])) %>% 
  c(list(Channel = c(
    'EpCAM', 'EGFR', 'BetaCatenin', 
    'Vimentin', 'CollagenI', 'AlphaSMA', 'CD31', 
    'CD45', 'CD11b', 'CD15', 'CD14', 'CD68', 'CD163', 'CD11c', 'HLADR', 
    'CD3', 'CD4', 'CD8a', 'CD45RO', 'FoxP3', 'CD56', 'GranzymeB','CD20', 'CD38', 'CD138', 
    'ICOS', 'PD1', 'LAG3', 'Tim3', 'Ki67', 'p53', 'PDL1', 'LMP1', 'Caspase3cleaved')))

ClinicalInfo_ROI <- data.table::fread('ClinicalInfo.csv')

ClinicalInfo_Pat <- ClinicalInfo_ROI %>% 
  dplyr::select(-ROI_ID, -TLS_ROI) %>% 
  dplyr::distinct()

# Figure S5a
pdf('FigS5a.pdf', width = 25, height = 6)
MT_CT %>% 
  dplyr::filter(ROI_ID %in% c('T039_ROI2')) %>% 
  dplyr::mutate(Y_position = -Y_position) %>%  
  {patchwork::wrap_plots(
    nrow = 1,
    guides = 'collect',
    ggplot(., aes(x = X_position, y = Y_position, group = -1L)) +
      geom_point(aes(color = CT0), show.legend = T,
                 bound = data.frame(x = c(0, 500, 500, 0),
                                    y = c(0, 0, -500, -500))) +
      scale_color_manual(values = c('#3B597AFF', 'grey70') %>% 
                           setNames(LEVELS$CT0),
                         guide = 'none'),
    
    ggplot(., aes(x = X_position, y = Y_position, group = -1L)) +
      ggforce::geom_voronoi_tile(aes(fill = SynoraAnnotation), show.legend = T,
                                 bound = data.frame(x = c(0, 500, 500, 0),
                                                    y = c(0, 0, -500, -500))) +
      scale_fill_manual(values = c("#CD3122", "#BEE183", "#2f7ab9") %>% 
                          setNames(c('Outside', 'Boundary', 'Nest')),
                        guide = 'none', na.value = 'grey70', drop = F),
    
    ggplot(., aes(x = X_position, y = Y_position, group = -1L)) +
      ggforce::geom_voronoi_tile(aes(fill = as.factor(Layer)), show.legend = T,
                                 bound = data.frame(x = c(0, 500, 500, 0),
                                                    y = c(0, 0, -500, -500))) +
      scale_fill_manual(values = viridis::turbo(n = 21, direction = -1) %>% 
                          setNames(as.factor(-10:10)), 
                        guide = 'none', na.value = 'grey70', drop = F),
    
    ggplot(., aes(x = X_position, y = Y_position, group = -1L)) +
      ggforce::geom_voronoi_tile(aes(fill = SA), show.legend = T,
                                 bound = data.frame(x = c(0, 500, 500, 0),
                                                    y = c(0, 0, -500, -500))) +
      scale_fill_manual(values = c("#762A83", "#1AE4B6", "#1965B0") %>% 
                          setNames(LEVELS$SA[1:3]), 
                        guide = 'none', na.value = 'grey70', drop = F)
  )}&
  theme_void() &
  theme(plot.margin = unit(c(0,0,0,0), "mm")) & 
  coord_equal()
dev.off()


# Figure S5b
pdf('FigS5b.pdf', height = 4.5, width = 6)
MT_CT %>% 
  dplyr::filter(!is.na(Layer)) %>% 
  dplyr::mutate(Layer = Layer %>% as.factor()) %>% 
  ggplot() +
  geom_bar(aes(x = Layer, fill = SA), stat = 'count') +
  ggprism::theme_prism(border = T, base_size = 12, base_fontface = 'plain') +
  theme(aspect.ratio = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     breaks = c(0, 1e5, 2e5)) +
  scale_x_discrete(expand = expansion(mult = c(0.01, 0.01)),
                   breaks = c('-10', '-5', '0', '5', '10')) +
  scale_fill_manual(values = c("#762A83", "#1AE4B6", "#1965B0")) +
  theme(legend.position = 'none') +
  labs(y = NULL)

dev.off()
