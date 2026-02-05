setwd("~/projects/NPC_SpatialProteomics/")
library(tidyverse)
library(data.table)
source('0.Helpers.R')

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


# Figure 2a
pdf('Fig2a.pdf', width = 15, height = 6)
MT_CT %>% 
  dplyr::filter(ROI_ID %in% c('T063_ROI2', 'T174_ROI1', 'T057_ROI2')) %>% 
  dplyr::mutate(Y_position = -Y_position) %>%  
  dplyr::nest_by(ROI_ID) %>% 
  dplyr::arrange(ROI_ID %>% match(c('T063_ROI2', 'T174_ROI1', 'T057_ROI2'))) %>% 
  dplyr::mutate(points_sf = data %>% 
                  sf::st_as_sf(coords = c("X_position", "Y_position"), remove = F) %>% 
                  list(),
                bbox = matrix(c(0, 0, 
                                500, 0, 
                                500, -500, 
                                0, -500, 
                                0, 0), 
                              ncol = 2, byrow = TRUE) %>% 
                  list() %>% 
                  sf::st_polygon() %>% 
                  sf::st_sfc() %>% 
                  list(),
                points_voronoi = points_sf %>%
                  sf::st_union() %>%
                  sf::st_voronoi(envelope = bbox) %>%
                  sf::st_collection_extract("POLYGON") %>%
                  sf::st_sf() %>%
                  sf::st_join(points_sf) %>% 
                  list(),
                points_buffered = points_sf %>%
                  sf::st_buffer(dist = 10) %>%
                  sf::st_union() %>% 
                  list(),
                points_final = points_voronoi %>% 
                  sf::st_intersection(points_buffered) %>%
                  sf::st_intersection(bbox) %>%
                  list()
  ) %>% 
  dplyr::arrange(ROI_ID %>% match(c('T109_ROI1', 'T109_ROI2'))) %>% 
  dplyr::mutate(
    PLOT1 = {points_final %>% 
        ggplot() +
        geom_sf(aes(fill = CT3), color = NA) +
        scale_fill_manual(name = ROI_ID,
                          values =
                            c('#3B597AFF',
                              paletteer::paletteer_d("Redmonder::sPBIGn", n = 3)[-1],
                              paletteer::paletteer_d("Redmonder::sPBIYl", n = 6)[-1],
                              paletteer::paletteer_d("Redmonder::sPBIRd", n = 8)[-1],
                              'grey70') %>%
                            setNames(LEVELS$CT3),
                          drop = F)
    } %>% list(),
    PLOT2 = {points_final %>% 
        ggplot() +
        geom_sf(aes(fill = SA), color = NA) +
        scale_fill_manual(values = c("#762A83", "#1AE4B6", "#1965B0") %>%
                            setNames(LEVELS$SA[1:3]),
                          na.value = 'grey70', drop = F)
    } %>% list()
  ) %>% 
  dplyr::mutate(PLOTS = {patchwork::wrap_plots(PLOT1, PLOT2, nrow = 1, guides = "collect") &
      theme_void() &
      guides(fill = guide_legend(ncol = 1)) &
      theme(plot.margin = unit(c(0,0,0,0), "mm"))} %>% 
        list()
  ) %>% 
  dplyr::pull(PLOTS) %>% 
  purrr::walk(print)
dev.off()


# Figure 2b
pdf('Fig2b.pdf', width = 10, height = 5)
patchwork::wrap_plots(
  MT_CT %>% 
    dplyr::filter(CT0 == 'Epi') %>% 
    dplyr::group_by(ROI_ID) %>% 
    dplyr::summarize(dplyr::across(c('EpCAM', 'EGFR', 'BetaCatenin', 'HLADR', 'CD138', 'Ki67', 'p53', 'PDL1', 'Caspase3cleaved'), .fns = mean)) %>% 
    tidyr::pivot_longer(cols = !ROI_ID, names_to = 'Channel', values_to = 'Value') %>% 
    dplyr::left_join(MT_CT %>% 
                       dplyr::distinct(ROI_ID, BNR) %>% 
                       dplyr::mutate(logBNR = log10(BNR)),
                     by = 'ROI_ID') %>% 
    dplyr::group_by(Channel) %>%
    rstatix::cor_test(vars = 'Value', vars2 = 'BNR', method = 'spearman') %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(cor) %>% 
    dplyr::mutate(Channel = Channel %>% forcats::fct_inorder()) %>% 
    dplyr::mutate(Label = paste0('cor = ', cor, '\np = ', ifelse(p > 0.001, sprintf('%.4f', p), sprintf('%.2e', p)))) %>% 
    ggplot(aes(y = Channel)) +
    geom_rect(aes(ymin = as.integer(Channel) - 0.5, ymax = as.integer(Channel) + 0.5, 
                  xmin = 0, xmax = 1, fill = cor), color = 'black') + 
    ggprism::theme_prism(base_fontface = 'plain') +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_blank(), 
          legend.title = element_text(size = 10),
          legend.position = 'left') +
    scale_fill_gradient2(name = 'Spearman\'s rho', high = '#CB2314', mid = 'white', low = '#046C9A', midpoint = 0) +
    labs(y = NULL, x = NULL) +
    coord_fixed(ratio = 1, clip = 'off') + 
    geom_text(aes(x = 1.2, label = Label), hjust = 0),
  
  MT_CT %>%
    dplyr::filter(CT0 == 'Epi') %>%
    dplyr::group_by(ROI_ID) %>% 
    dplyr::summarize(dplyr::across(c('BetaCatenin', 'EGFR', 'EpCAM', 'CD138'), .fns = mean)) %>% 
    tidyr::pivot_longer(cols = !ROI_ID, names_to = 'Channel', values_to = 'Value') %>% 
    dplyr::mutate(Channel = Channel %>% forcats::fct_relevel(c('BetaCatenin', 'EGFR', 'EpCAM', 'CD138'))) %>% 
    dplyr::left_join(MT_CT %>% 
                       dplyr::distinct(ROI_ID, BNR) %>% 
                       dplyr::mutate(logBNR = log10(BNR)),
                     by = 'ROI_ID') %>% 
    ggpubr::ggscatter(x = 'logBNR', y = "Value", color = "#6a8fad", size = 2,
                      add = "reg.line", conf.int = T, add.params = list(color = "black", fill = "lightgrey")
    ) +
    facet_wrap(.~Channel, ncol = 2, scales = 'free_y', axes = 'all') +
    ggprism::theme_prism(base_fontface = 'plain', border = T, base_size = 10) +
    theme(aspect.ratio = 1, strip.text = element_text(size = 12)) +
    labs(y = 'Mean Intensity (Z-score)'),
  nrow = 1
)
dev.off()




# Figure 2c
pdf('Fig2c.pdf', width = 9, height = 5)
MT_CT %>% 
  dplyr::filter(CT0 == 'TME') %>%
  dplyr::mutate(CT3 = CT3 %>% forcats::fct_drop()) %>%
  Helper1(UNIT_ID = 'ROI_ID', TERM_A = 'CT3') %>% 
  dplyr::filter(Term_A != 'TME_unknown') %>% 
  dplyr::left_join(MT_CT %>% 
                     Helper1('ROI_ID', 'SA') %>% 
                     dplyr::filter(Feature == 'LE') %>% 
                     dplyr::transmute(ROI_ID, LE = Value),
                   by = 'ROI_ID') %>% 
  dplyr::left_join(MT_CT %>% 
                     dplyr::distinct(ROI_ID, BNR) %>% 
                     dplyr::mutate(logBNR = log10(BNR)),
                   by = 'ROI_ID') %>% 
  dplyr::group_by(Term_A) %>% 
  rstatix::cor_test(vars = 'Value', vars2 = c('LE', 'BNR'), method = 'spearman') %>% 
  tidyr::pivot_wider(id_cols = 'Term_A', names_from = 'var2', values_from = c('cor', 'p')) %>% 
  dplyr::mutate(p_LE = ifelse(p_LE < 0.05, p_LE, NA),
                p_BNR = ifelse(p_BNR < 0.05, p_BNR, NA)) %>% 
  ggplot(aes(x = !!as.name(paste0('cor_', 'BNR')), y = !!as.name(paste0('cor_', 'LE')))) + 
  geom_hline(yintercept = c(-0.1, 0.1), linetype = 'dotted', linewidth = 0.2) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = 'dotted', linewidth = 0.2) +
  geom_point(aes(fill = - sign(!!as.name(paste0('cor_', 'BNR'))) * pmax(log10(!!as.name(paste0('p_', 'BNR'))), -3), 
                 color = - sign(!!as.name(paste0('cor_', 'LE'))) * pmax(log10(!!as.name(paste0('p_', 'LE'))), -3)), 
             shape = 21, size = 4, stroke = 3) + 
  ggrepel::geom_text_repel(aes(label = Term_A), 
                           max.overlaps = Inf, min.segment.length = 0, size = 4) +
  scale_color_gradientn(colours = c('#046C9A', 'white', '#CB2314'),
                        breaks = c(0.001, 0.01, 0.05, 1/0.05, 1/0.01, 1/0.001) %>% log10(),
                        labels = c('< 0.001', '0.01', '0.05', '0.05', '0.01', '0.001'),
                        name = 'P-value\nInner circle (BNR)\nOuter circle(LE)', na.value = 'grey50') +
  scale_fill_gradientn(guide = 'none', colours = c('#046C9A', 'white', '#CB2314'),
                       breaks = c(0.001, 0.01, 0.05, 1/0.05, 1/0.01, 1/0.001) %>% log10(), 
                       na.value = 'grey50') +
  scale_x_continuous(breaks = -4:4/10) +
  scale_y_continuous(breaks = -4:4/10) +
  ggprism::theme_prism(border = T, base_size = 12, base_line_size = 0.5, base_fontface = 'plain') +
  theme(legend.title = element_text(), aspect.ratio = 0.75)
dev.off()






# Figure 2d
CancerCellDensity <- MT_CT %>% 
  dplyr::filter(CT0 %in% c('Epi')) %>%
  dplyr::filter(SA %in% c('LE', 'Tumor')) %>% 
  dplyr::filter(!is.na(Distance2Boundary)) %>% 
  dplyr::nest_by(CT4) %>% 
  dplyr::mutate(CT4 = CT4 %>% forcats::fct_relabel(stringr::str_remove, 'Epi_')) %>%
  dplyr::mutate(data_bg = .$data %>% dplyr::bind_rows() %>% list()) %>%
  dplyr::mutate(ks_test = ks.test(data$Distance2Boundary, data_bg$Distance2Boundary) %>%
                  broom::glance()) %>% 
  dplyr::mutate(dens = density(data$Distance2Boundary) %>%
                  list(),
                x = dens$x %>% list(),
                y = dens$y %>% list(),
                peak_x = dens$x[which.max(dens$y)],
                peak_y = max(dens$y)
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(CT4 %in% c('p53', 'HLADR', 'PDL1', 'ki67', 'BetaCatenin', 'CD138')) %>%
  dplyr::mutate(CT4 = CT4 %>% forcats::fct_relevel(c('p53', 'HLADR', 'PDL1', 'ki67', 'BetaCatenin', 'CD138'))) %>%
  tidyr::unnest(ks_test)


pdf('Fig2d.pdf', height = 5, width = 6)
ggplot() +
  annotate(geom = 'rect', 
           xmin = -Inf, xmax = 20, 
           ymin = -Inf, ymax = Inf, 
           fill = '#1AE4B630') +
  annotate(geom = 'rect', 
           xmin = 20, xmax = Inf, 
           ymin = -Inf, ymax = Inf,
           fill = '#1965B030') +
  geom_segment(data = CancerCellDensity %>% dplyr::distinct(CT4, peak_x, peak_y),
               aes(x = -20, xend = peak_x, y = peak_y, yend = peak_y, color = CT4),
               linetype = 'dashed') +
  geom_segment(data = CancerCellDensity %>% dplyr::distinct(CT4, peak_x, peak_y),
               aes(x = peak_x, xend = peak_x, y = 0, yend = peak_y, color = CT4),
               linetype = 'dashed') +
  geom_line(data = CancerCellDensity %>% 
              dplyr::select(CT4, x, y) %>% 
              tidyr::unnest(c(x, y)),
            aes(x = x, y = y, color = CT4),
            linewidth = 1) +
  geom_point(data = CancerCellDensity %>% dplyr::distinct(CT4, peak_x, peak_y),
             aes(x = peak_x, y = peak_y, color = CT4),
             size = 3, show.legend = FALSE) +
  ggprism::theme_prism(border = T, base_fontface = 'plain', base_size = 12) +
  scale_x_continuous(breaks = (-1:5) * 20,
                     expand = expansion(c(0, 0.02))) +
  scale_y_continuous(expand = expansion(c(0, 0.1))) +
  scale_color_manual(
    name = 'Kolmogorov-Smirnov statistic',
    values = c('#6351A0', '#88775F', '#6B8993',
               '#DAA520', '#E7695D', '#D070B9') %>%
      setNames(CancerCellDensity %>%
                 dplyr::distinct(CT4, peak_x, peak_y, statistic) %>%
                 dplyr::pull(CT4)),
    labels = CancerCellDensity %>%
      dplyr::arrange(CT4) %>% 
      dplyr::distinct(CT4, statistic) %>%
      dplyr::mutate(label = paste0(CT4, ': ', round(statistic, 3))) %>%
      dplyr::pull(label)
  ) +
  theme(legend.title = element_text(),
        legend.position = c(1, 1),
        legend.justification = c(1, 1)) + 
  coord_cartesian(xlim = c(-20, 100)) +
  labs(x = "Distance to Boundary (um)", y = "Cancer cell density")
dev.off()




# Figure 2e
temp <- MT_CT %>% 
  dplyr::filter(CT0 == 'TME') %>%
  dplyr::filter(!is.na(Layer)) %>% 
  dplyr::transmute(Layer = Layer %>% forcats::as_factor(), 
                   CT3 = CT3 %>% forcats::fct_drop(), 
                   CT4 = CT4 %>% forcats::fct_drop())

OR_result <- as.factor(-10:10) %>% 
  purrr::map(\(X) {
    temp %>% 
      dplyr::mutate(Layer = ifelse(Layer == X, 'Y', 'N') %>% 
                      forcats::fct_expand(c('Y', 'N')) %>% 
                      forcats::fct_relevel(c('Y', 'N'))) %>% 
      dplyr::select(CT4, Layer) %>% 
      table() %>% 
      rstatix::row_wise_fisher_test(detailed = T) %>% 
      dplyr::transmute(CT = group,
                       Layer = X,
                       OR = estimate)
  })

pdf('Fig2e.pdf', height = 10, width = 8)
OR_result %>% 
  dplyr::bind_rows() %>% 
  dplyr::group_by(CT) %>% 
  dplyr::mutate(OR_lag = OR %>% lag()) %>% 
  dplyr::mutate(OR_lead = OR %>% lead()) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(OR = mean(c(OR, OR_lag, OR_lead), na.rm = T)) %>% 
  dplyr::group_by(CT) %>% 
  dplyr::mutate(Label = ifelse(OR == max(OR), sprintf('%.2f', OR), '')) %>%
  dplyr::ungroup() %>%
  dplyr::filter(CT != 'TME_unknown') %>% 
  dplyr::mutate(CT = CT %>%
                  forcats::fct_relevel(c('PC', 'CD8_T_Tim3', 'vCAF', 'Granulo', 'CAF_col_low', 'NK', 'CD8_T_other',
                                         'DC', 'TME_unknown', 'CD4_Treg1', 'Mono', 'M1', 'CAF_col_hi', 'M2', 'Endo',
                                         'DNT', 'CD8_T_LAG3', 'CD4_Treg2_ICOS', 'CD8_T_PD1', 'CD4_T_ICOS', 'CD4_T_Tim3',
                                         'CD8_T_GranzymeB_2', 'CD8_T_GranzymeB_1', 'BC', 'CD8_T_CD45RO', 
                                         'CD4_Tn', 'CD4_T_PD1', 'CD4_T_CD45RO'))) %>% 
  ggplot(aes(y = CT, x = Layer)) + 
  geom_tile(aes(fill = OR)) +
  geom_text(aes(label = Label), size = 3) + 
  scale_fill_gradient2(low = '#5B8DD8', mid = "#FFFFFF", high = '#CB2314',
                       midpoint = 1, transform = 'log10') +
  scale_x_discrete(expand = expansion()) + 
  scale_y_discrete(expand = expansion()) + 
  ggprism::theme_prism(border = T, base_fontface = 'plain', base_size = 10) +
  theme(legend.title = element_text()) +
  labs(y = NULL) +
  coord_fixed(ratio = 1)
dev.off()
