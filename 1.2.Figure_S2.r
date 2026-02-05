setwd("~/projects/NPC_SpatialProteomics/")
library(tidyverse)
library(data.table)
library(ComplexHeatmap)
source('0.Helpers.R')

MT_CT <- qs::qread('SingleCellDataFrame.qs')
LEVELS <- qs::qread('AnnotationLevels.qs')
ClinicalInfo_ROI <- data.table::fread('ClinicalInfo.csv')
ClinicalInfo_Pat <- ClinicalInfo_ROI %>% 
  dplyr::select(-ROI_ID, -TLS_ROI) %>% 
  dplyr::distinct()


PALETTE1 <- c('#005F73', 
              '#04707F', '#0A9396',
              '#9B9B7A', '#baa587', '#d9ae94', '#d08c60', '#b58463',
              '#E9D8A6', '#ECCE90', '#EFC479', '#F0B962', '#F0AF4A', '#EFA52F', '#EE9B00',
              'grey70') %>% 
  setNames(LEVELS$CT3)
PALETTE2 <- c('#6c757d', '#6A3D9A') %>% setNames(c('0', '1'))
PALETTE3 <- c('#046C9A', '#CB2314') %>% setNames(c('0', '1'))
HEIGHT1 <- unit(0.4, "cm")
HEIGHT2 <- unit(10, "cm")
WIDTH1 <- unit(0.15, "cm")

# Figure S2b
pdf('FigS2b.pdf', height = 5, width = 27)
MT_CT %>%
  Helper1('ROI_ID', 'CT3') %>% 
  tidyr::pivot_wider(id_cols = 'ROI_ID', names_from = 'Feature', values_from = 'Value', values_fill = 0) %>% 
  dplyr::left_join(ClinicalInfo_ROI[,.(ROI_ID, EFS, InductChemo)], by = 'ROI_ID') %>% 
  dplyr::arrange(EFS, InductChemo, desc(Epi)) %>% 
  dplyr::nest_by(EFS, InductChemo, .keep = T) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(Plot = data %>% 
                  {Heatmap(matrix = matrix(0, nrow = 0, ncol = nrow(.)), cluster_rows = F, cluster_columns = F,
                           row_title = NULL, height = 0, width = WIDTH1 * nrow(.), show_heatmap_legend = F,
                           top_annotation = HeatmapAnnotation(
                             gap = unit(1, "mm"),
                             '  ' = anno_simple(.$InductChemo, height = HEIGHT1, col = PALETTE2),
                             ' ' = anno_simple(.$EFS, height = HEIGHT1, col = PALETTE3),
                             '   ' = anno_barplot(dplyr::select(., LEVELS$CT3),
                                                  gp = gpar(fill = PALETTE1[LEVELS$CT3],
                                                            col = NA),
                                                  bar_width = 0.8, height = HEIGHT2, axis = F)
                           )
                  )} %>%
                  list()
  ) %>% 
  dplyr::select(InductChemo, Plot) %>%
  dplyr::nest_by(InductChemo) %>%
  dplyr::mutate(Plots = data$Plot %>%
                  purrr::reduce(`+`) %>%
                  draw(annotation_legend_list = list(
                    Legend(title = 'Event', labels =  c('0', '1'),
                           legend_gp = gpar(fill = PALETTE3)),
                    Legend(title = 'Induction\nChemotherapy', labels =  c('0', '1'),
                           legend_gp = gpar(fill = PALETTE2)),
                    Legend(title = 'Cell Type', labels = LEVELS$CT3,
                           legend_gp = gpar(fill = PALETTE1[LEVELS$CT3]))
                  )) %>%
                  list()) %>% 
  dplyr::pull(Plots)
dev.off()
