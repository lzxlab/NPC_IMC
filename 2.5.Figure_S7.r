setwd('~/projects/NPC_SpatialProteomics/')
library(tidyverse)
library(data.table)
library(Seurat)
library(monocle3)

MT_CT <- qs::qread('SingleCellDataFrame.qs')
LEVELS <- qs::qread('AnnotationLevels.qs')
ClinicalInfo_ROI <- data.table::fread('ClinicalInfo.csv')
ClinicalInfo_Pat <- ClinicalInfo_ROI %>% 
  dplyr::select(-ROI_ID, -TLS_ROI) %>% 
  dplyr::distinct()


# Figure S7a
cds_CD8_T <- qs::qread('cds_CD8_T.qs')
pdf('FigS7a.pdf', width = 8, height = 7)
cds_CD8_T %>% 
  monocle3::plot_cells(color_cells_by = "pseudotime",  
             trajectory_graph_color = "black",
             trajectory_graph_segment_size = 1.5,
             label_principal_points = F,
             label_branch_points = F,
             label_roots = F,
             label_leaves = F,
             show_trajectory_graph = T)
dev.off()

# Figure S7b
pdf('FigS7b.pdf', width = 12, height = 3)
c('PDCD1', 'LAG3', 'HAVCR2', 'CCR7') %>% 
  purrr::map(.progress = T, \(X) {
    monocle3::plot_genes_in_pseudotime(cds_CD8_T[X,], label_by_short_name = F)$data
  }) %>% 
  dplyr::bind_rows() %>% 
  ggplot(aes(x = pseudotime, y = expectation)) + 
  facet_wrap(.~ feature_label, nrow = 1, scales = 'free') +
  geom_line() +
  ggprism::theme_prism(base_fontface = 'plain', border = T)
dev.off()



# Figure S7c
IMC_CD8_T <- MT_CT %>%
  dplyr::filter(CT3 == 'CD8_T') %>% 
  dplyr::mutate(APT = PD1 %>%
                  dplyr::dense_rank() %>%
                  scales::rescale(to = c(0, 1))) %>% 
  dplyr::filter(!is.na(Layer))

pdf('FigS7c.pdf', width = 11, height = 3)
IMC_CD8_T %>% 
  tidyr::pivot_longer(cols = c('PD1', 'LAG3', 'Tim3', 'CD45RO'), 
                      names_to = 'Channel', 
                      values_to = 'Expression') %>% 
  dplyr::mutate(Channel = Channel %>% 
                  forcats::fct_relevel(c('PD1', 'LAG3', 'Tim3', 'CD45RO'))) %>% 
  dplyr::nest_by(Channel) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Prediction =
                  speedglm::speedglm(data =  data, Expression ~ 
                                       splines::ns(APT, df = 3),
                                     model = F, y = F, fitted = T) %>%
                  predict(type = "response") %>%
                  list()) %>% 
  tidyr::unnest(cols = c(data, Prediction)) %>% 
  ggplot(aes(x = APT)) +
  facet_wrap(.~Channel, nrow = 1, scales = 'free_y') +
  geom_line(aes(y = Prediction, group = Channel)) +
  labs(x = 'Activation Pseudotime', y = 'Expression') +
  ggprism::theme_prism(base_fontface = 'plain', border = T)
dev.off()





# Figure S7d
IMC_CD8_T_embeddings <- qs::qread('UMAP_DATA.qs') %>% dplyr::as_tibble()

pdf('FigS7d.pdf', width = 9, height = 6)
IMC_CD8_T_embeddings %>% 
  dplyr::mutate(Layer = IMC_CD8_T$Layer) %>% 
  ggplot() +
  geom_point(aes(x = V1, y = V2, color = Layer), size = 0.1) +
  theme_void(base_size = 12) + 
  theme(panel.spacing = unit(10, 'line')) +
  labs(x = NULL, y = NULL, title = NULL) +
  scale_color_gradientn(colors = viridis::turbo(n = 7)) +
  coord_equal()
IMC_CD8_T_embeddings %>% 
  dplyr::mutate(APT = IMC_CD8_T$APT) %>% 
  ggplot() +
  geom_point(aes(x = V1, y = V2, color = APT), size = 0.1) +
  theme_void(base_size = 12) + 
  labs(x = NULL, y = NULL, title = NULL) +
  scale_colour_viridis_c(name = 'Activation Pseudotime', option = "plasma") +
  coord_equal()
dev.off()



# Figure S7e
pdf('FigS7e.pdf', height = 4, width = 4)
IMC_CD8_T %>%
  dplyr::arrange(desc(abs(Layer))) %>% 
  ggplot(aes(x = SA, y = APT, fill = SA)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(method = 'wilcox', 
                             comparisons = list(c('Stroma', 'LE'), c('LE', 'Tumor'), c('Stroma', 'Tumor'))) + 
  scale_fill_manual(values = c("#762A83", "#1AE4B6", "#1965B0") %>% 
                      setNames(LEVELS$SA[1:3]), na.value = 'grey70', drop = F, guide = 'none') +
  scale_y_continuous(expand = expansion(c(0.01, 0.1))) +
  ggprism::theme_prism(base_fontface = 'plain', border = T, base_line_size = 0.5) +
  labs(x = NULL, y = 'Activation Pseudotime')
dev.off()







