setwd('~/projects/NPC_SpatialProteomics/')
library(tidyverse)
library(data.table)


MT_CT <- qs::qread('SingleCellDataFrame.qs')

LEVELS <- c('CT0', 'CT1', 'CT2', 'CT3', 'CT4', 'SA') %>% 
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



# Figure S9a
DDE_ICOS_Tim3 <- qs::qread('DDE_ICOS_Tim3.qs')
pdf('FigS9a.pdf', width = 10, height = 5)
DDE_ICOS_Tim3 %>% 
  dplyr::arrange(CT_Key, CT_Query, Channel) %>% 
  dplyr::group_by(CT_Query, Channel) %>% 
  dplyr::mutate(Label = ifelse(Radius == 8 & p < 0.05 & abs(ES) > 0.25, as.character(CT_Query), NA)) %>% 
  dplyr::mutate(ES = ifelse(p < 0.05, ES, NA)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(Panel = paste0(Channel, ' on ', CT_Key) %>% forcats::fct_inorder()) %>% 
  ggplot(aes(x = Radius, y = ES)) +
  facet_wrap(. ~ Panel, nrow = 1, axes = 'all') + 
  geom_hline(yintercept = c(0), color = 'grey70') +
  geom_hline(yintercept = c(-0.4, -0.2, 0.2, 0.4), linetype = 'dotted', color = 'grey70') +
  geom_line(aes(color = CT_Query, group = CT_Query)) +
  ggrepel::geom_text_repel(aes(label = Label),
                           nudge_x = -2, max.overlaps = Inf, size = 3,
                           min.segment.length = 0,
                           direction = 'y', hjust = 1)+
  scale_x_continuous(breaks = 0:6 * 8, 
                     expand = expansion(c(0.15, 0.05))) +
  scale_y_continuous(breaks = c(-2:2 * 0.2)) +
  scale_color_manual(values = c('#3B597AFF', 
                                paletteer::paletteer_d("Redmonder::sPBIGn", n = 3)[-1],
                                paletteer::paletteer_d("Redmonder::sPBIYl", n = 6)[-1],
                                paletteer::paletteer_d("Redmonder::sPBIRd", n = 8)[-1],
                                'grey70') %>% 
                       setNames(LEVELS$CT3)) +
  labs(x = 'Distance (um)', y = 'Effect Size (rc)') +
  ggprism::theme_prism(base_fontface = 'plain', border = T) + 
  theme(strip.text = element_text(size = 12))
dev.off()



# Figure S9b
ST_Bei <- qs::qread('ST_Bei.qs')

pdf('FigS9b.pdf', width = 12, height = 10)
ST_Bei[c('GSM6248646_NPC_ST5', 'GSM6248657_NPC_ST19')] %>% 
  purrr::imap(\(X, Y) {
    Meta_New <- paste0(
      ifelse(colSums(X@assays$Spatial$count[c('CD3D', 'CD3E'),]) > 0 & 
               colSums(X@assays$Spatial$count[c('CD8A', 'CD8B'),]) > 0, 
             'CD8_T present', 'CD8_T absent'),
      ', ',
      ifelse(colSums(X@assays$Spatial$count[c('EPCAM', 'KRT13', 'KRT8', 'KRT5'),] == 0) == 0, 
             'TCA', 'non-TCA')) %>% 
      forcats::fct_relevel('CD8_T absent, non-TCA', 'CD8_T absent, TCA', 
                           'CD8_T present, non-TCA', 'CD8_T present, TCA')
    X <- X %>% AddMetaData(Meta_New, col.name = 'Annotation')
    
    return(patchwork::wrap_plots(
      ncol = 1,
      SpatialDimPlot(X, 'Annotation', pt.size.factor = ifelse(Y == 'GSM6248646_NPC_ST5', 3, 4), image.alpha = 0) + 
        scale_fill_manual(values = c('#FEF7C7', '#B0CBE7', '#DABD61', '#3A488A')) + 
        theme(plot.title = element_text(hjust = 0.5)) + 
        guides(fill = guide_legend(override.aes = list(size = 6))) + 
        labs(title = Y),
      SpatialFeaturePlot(X, 'HAVCR2', pt.size.factor = ifelse(Y == 'GSM6248646_NPC_ST5', 3, 4), image.alpha = 0, slot = 'counts') + 
        theme(plot.title = element_text(hjust = 0.5))
    ))
  }) %>% 
  patchwork::wrap_plots(nrow = 1, guides = 'collect')
dev.off()


ST_Bei_exp <- ST_Bei %>% 
  purrr::imap(\(X, Y) {
    if (typeof(X) == 'S4') {
      X@assays$Spatial$count[c('EPCAM', 'KRT13', 'KRT8', 'KRT5', 
                               'CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'HAVCR2', 
                               'CD68', 'CD14', 'CD163'),] %>% 
        t() %>%
        as.data.frame.matrix() %>%
        tibble::rownames_to_column('Cell_ID') %>%
        dplyr::mutate(Dataset = Y)
    } else {
      X %>% 
        dplyr::filter(V1 %in% c('EPCAM', 'KRT13', 'KRT8', 'KRT5', 
                                'CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'HAVCR2', 
                                'CD68', 'CD14', 'CD163')) %>% 
        tibble::column_to_rownames('V1') %>%
        t() %>%
        as.data.frame.matrix() %>% 
        tibble::rownames_to_column('Cell_ID') %>% 
        dplyr::mutate(Dataset = Y)
    }
  }) %>% 
  dplyr::bind_rows()%>% 
  dplyr::as_tibble()


# Figure S9c
pdf('FigS9c.pdf', width = 4, height = 4)
ST_Bei_exp %>% 
  dplyr::filter((CD3D + CD3E) > 0 & (CD8A + CD8B) > 0) %>% 
  dplyr::mutate(TCA_group = ifelse((EPCAM > 0) & (KRT13 > 0) & (KRT8 > 0) & (KRT5 > 0), 'TCA', 'non-TCA')) %>% 
  dplyr::group_by(Dataset, TCA_group) %>% 
  dplyr::summarize(`Mean HAVCR2 Count on CD3+CD8+ spots` = mean(HAVCR2)) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = TCA_group, y = `Mean HAVCR2 Count on CD3+CD8+ spots`)) +
  geom_boxplot(aes(fill = TCA_group), outliers = F, width = 0.5) +
  geom_point() +
  geom_line(aes(group = Dataset), alpha = 0.3) + 
  ggpubr::stat_compare_means(comparisons = list(c('TCA', 'non-TCA')), method = 'wilcox', paired = T) + 
  ggprism::theme_prism(border = T, base_fontface = 'plain') + 
  labs(x = NULL, y = 'Mean HAVCR2 Count\non CD3+CD8+ spots') + 
  scale_y_continuous(expand = expansion(0.01, 0.1)) +
  scale_fill_manual(values = c('#FEF7C7', '#5D74A5'), guide = 'none')
dev.off()




# Figure S9d
pdf('FigS9d.pdf', width = 4, height = 4)
ST_Bei_exp %>% 
  dplyr::filter((CD3D + CD3E) > 0 & (CD8A + CD8B) > 0) %>% 
  dplyr::mutate(Mph_profile = dplyr::case_when(
    CD14 > 0 & CD68 > 0 & CD163 == 0 ~ 'CD14+CD68+\nCD163- spots',
    CD14 > 0 & CD68 > 0 & CD163 > 0 ~ 'CD14+CD68+\nCD163+ spots',
    T ~ 'Other spots'
  )) %>% 
  dplyr::group_by(Dataset, Mph_profile) %>% 
  dplyr::summarize(`Mean HAVCR2 Count on CD3+CD8+ spots` = mean(HAVCR2)) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = Mph_profile, y = `Mean HAVCR2 Count on CD3+CD8+ spots`)) +
  geom_boxplot(aes(fill = Mph_profile), outliers = F, width = 0.5) +
  geom_point() +
  geom_line(aes(group = Dataset), alpha = 0.3) + 
  ggpubr::stat_compare_means(comparisons = list(
    c('CD14+CD68+\nCD163- spots', 'CD14+CD68+\nCD163+ spots'), 
    c('CD14+CD68+\nCD163+ spots', 'Other spots'), 
    c('CD14+CD68+\nCD163- spots', 'Other spots')), 
    method = 'wilcox', paired = T) + 
  ggprism::theme_prism(border = T, base_fontface = 'plain') + 
  theme(axis.text.x = element_text(vjust = 1, size = 8)) +
  labs(x = NULL, y = 'Mean HAVCR2 Count\non CD3+CD8+ spots') + 
  scale_y_continuous(expand = expansion(0.01, 0.1)) +
  scale_fill_manual(values = c('#647D4B', '#7D96AF', '#E1E1E1'), guide = 'none')
dev.off()


# Figure S9d
pdf('FigS9e.pdf', width = 8, height = 8)
dplyr::left_join(
  ST_Bei_exp %>% 
    dplyr::mutate(CD8T_profile = dplyr::case_when(
      (CD3D + CD3E) > 0 & (CD8A + CD8B) > 0 & HAVCR2 > 0 ~ 'CD3+CD8+HAVCR2+ spots',
      (CD3D + CD3E) > 0 & (CD8A + CD8B) > 0 & HAVCR2 == 0 ~ 'CD3+CD8+HAVCR2- spots',
      T ~ 'CD8+T absent spots'
    )) %>% 
    dplyr::count(Dataset, CD8T_profile, name = 'n1'),
  ST_Bei_exp %>%  
    dplyr::mutate(TCA_group = ifelse((EPCAM > 0) & (KRT13 > 0) & (KRT8 > 0) & (KRT5 > 0), 'TCA', 'non-TCA')) %>% 
    dplyr::count(., Dataset, TCA_group, name = 'n2')) %>% 
  ggplot(aes(x = n1, y = n2)) +
  facet_wrap(CD8T_profile ~ TCA_group, scales = 'free', ncol = 2, strip.position = "left") +
  geom_point() +
  geom_smooth(color = "#1965B0", fill = "lightgray", method = 'lm') +
  ggprism::theme_prism(base_fontface = 'plain', border = T) +
  theme(strip.placement = 'outside') +
  ggpubr::stat_cor(method = 'spearman') +
  labs(x = NULL, y = NULL)
dev.off()


# Figure S9f
pdf('FigS9f.pdf', width = 12, height = 8)
dplyr::left_join(
  ST_Bei_exp %>% 
    dplyr::mutate(CD8T_profile = dplyr::case_when(
      (CD3D + CD3E) > 0 & (CD8A + CD8B) > 0 & HAVCR2 > 0 ~ 'CD3+CD8+HAVCR2+ spots',
      (CD3D + CD3E) > 0 & (CD8A + CD8B) > 0 & HAVCR2 == 0 ~ 'CD3+CD8+HAVCR2- spots',
      T ~ 'CD8+T absent spots'
    )) %>% 
    dplyr::count(Dataset, CD8T_profile, name = 'n1'),
  ST_Bei_exp %>%  
    dplyr::mutate(Mph_profile = dplyr::case_when(
      CD14 > 0 & CD68 > 0 & CD163 == 0 ~ 'CD14+CD68+CD163- spots',
      CD14 > 0 & CD68 > 0 & CD163 > 0 ~ 'CD14+CD68+CD163+ spots',
      T ~ 'Macrophage absent spots'
    )) %>% 
    dplyr::count(Dataset, Mph_profile, name = 'n2')) %>% 
  ggplot(aes(x = n1, y = n2)) +
  facet_wrap(CD8T_profile ~ Mph_profile, scales = 'free', ncol = 3, strip.position = "left") +
  geom_point() +
  geom_smooth(color = "#1965B0", fill = "lightgray", method = 'lm') +
  ggprism::theme_prism(base_fontface = 'plain', border = T) +
  theme(strip.placement = 'outside') +
  ggpubr::stat_cor(method = 'spearman') +
  labs(x = NULL, y = NULL)
dev.off()
