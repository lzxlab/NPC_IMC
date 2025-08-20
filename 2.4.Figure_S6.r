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



# Figure S6a
pdf('FigS6a.pdf', width = 10, height = 5)
patchwork::wrap_plots(
  MT_CT %>% 
    dplyr::filter(CT0 == 'Epi') %>% 
    dplyr::group_by(ROI_ID) %>% 
    dplyr::summarize(dplyr::across(c('EpCAM', 'EGFR', 'BetaCatenin', 'HLADR', 'CD138', 'Ki67', 'p53', 'PDL1', 'Caspase3cleaved'), 
                                   .fns = mean)) %>% 
    tidyr::pivot_longer(cols = !ROI_ID, names_to = 'Channel', values_to = 'Value') %>% 
    dplyr::left_join(MT_CT %>% 
                       Helper1('ROI_ID', 'SA') %>% 
                       dplyr::filter(Feature == 'LE') %>% 
                       dplyr::transmute(ROI_ID, LE = Value),
                     by = 'ROI_ID') %>% 
    dplyr::group_by(Channel) %>%
    rstatix::cor_test(vars = 'Value', vars2 = 'LE', method = 'spearman') %>% 
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
    dplyr::summarize(dplyr::across(c('EpCAM', 'CD138', 'BetaCatenin', 'EGFR'), .fns = mean)) %>% 
    tidyr::pivot_longer(cols = !ROI_ID, names_to = 'Channel', values_to = 'Value') %>% 
    dplyr::mutate(Channel = Channel %>% forcats::fct_relevel(c('EpCAM', 'CD138', 'BetaCatenin', 'EGFR'))) %>% 
    dplyr::left_join(MT_CT %>% 
                       Helper1('ROI_ID', 'SA') %>% 
                       dplyr::filter(Feature == 'LE') %>% 
                       dplyr::transmute(ROI_ID, LE = Value),
                     by = 'ROI_ID') %>% 
    ggpubr::ggscatter(x = 'LE', y = "Value", color = "#6a8fad", size = 2,
                      add = "reg.line", conf.int = T, add.params = list(color = "black", fill = "lightgrey")
    ) +
    facet_wrap(.~Channel, ncol = 2, scales = 'free_y', axes = 'all') +
    ggprism::theme_prism(base_fontface = 'plain', border = T, base_size = 10) +
    theme(aspect.ratio = 1, strip.text = element_text(size = 12)) +
    labs(y = 'Mean Intensity (Z-score)'),
  nrow = 1
)
dev.off()


# Figure S6b
pdf('FigS6b.pdf', width = 9, height = 5)
MT_CT %>% 
  dplyr::filter(CT0 == 'TME') %>%
  dplyr::mutate(CT4 = CT4 %>% forcats::fct_drop()) %>%
  Helper1(UNIT_ID = 'ROI_ID', TERM_A = 'CT4') %>% 
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


# Figure S6c
pdf('FigS6c.pdf', height = 4.5, width = 6)
MT_CT %>% 
  dplyr::filter(!is.na(Layer)) %>% 
  dplyr::mutate(Layer = Layer %>% as.factor()) %>% 
  dplyr::count(Layer, CT3) %>% 
  dplyr::group_by(Layer) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = Layer, y = prop, fill = CT3), color = 'white', 
           linewidth = 0.05, stat = 'identity', position = 'stack') +
  ggprism::theme_prism(border = T, base_size = 12, base_fontface = 'plain') +
  theme(aspect.ratio = 1) +
  scale_x_discrete(expand = expansion(mult = c(0.01, 0.01)),
                   breaks = c('-10', '-5', '0', '5', '10')) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)), 
                     breaks = c(0, 0.5, 1), 
                     labels = scales::percent_format()) +
  scale_fill_manual(values = c('#3B597AFF', 
                               paletteer::paletteer_d("Redmonder::sPBIGn", n = 3)[-1],
                               paletteer::paletteer_d("Redmonder::sPBIYl", n = 6)[-1],
                               paletteer::paletteer_d("Redmonder::sPBIRd", n = 8)[-1],
                               'grey70') %>% 
                      setNames(LEVELS$CT3)) +
  labs(y = NULL)
dev.off()


# Figure S6d
pdf('FigS6d.pdf', width = 8, height = 4)
patchwork::wrap_plots(nrow = 1, 
                      MT_CT %>% 
                        Helper2(UNIT_ID = 'Pat_ID', TERM_A = 'CT4', TERM_B = 'SA', REL2A = F, REL2B = T) %>% 
                        dplyr::filter(Term_B %in% c('LE', 'Tumor')) %>% 
                        dplyr::filter(Term_A == c('Epi_HLADR')) %>% 
                        dplyr::mutate(Value = ifelse(Prop4Filter >= 0.01, Value, 0)) %>% 
                        dplyr::left_join(ClinicalInfo_Pat_Binary) %>% 
                        dplyr::filter(OS %in% c(0, 1)) %>% 
                        dplyr::mutate(OS = ifelse(OS == 0, 'Alive', 'Dead') %>% as.factor()) %>% 
                        dplyr::mutate(Term_B = ifelse(Term_B == 'LE', paste0(Term_A, ' on LE'), paste0(Term_A, ' on Tumor')) %>% as.factor()) %>% 
                        ggplot(aes(x = OS, y = Value)) + 
                        facet_wrap(. ~ Term_B, scales = 'free_y') +
                        geom_boxplot(aes(fill = Term_B)) +
                        ggprism::theme_prism(border = T, base_size = 10, base_fontface = 'plain') +
                        scale_fill_manual(values = c('#6388b8', '#a5c561'), guide = 'none') + 
                        theme(aspect.ratio = 2.5, strip.text = element_text(size = 8)) +
                        ggpubr::stat_compare_means(method = 'wilcox', size = 3) +
                        labs(y = NULL), 
                      MT_CT %>% 
                        Helper2(UNIT_ID = 'Pat_ID', TERM_A = 'CT4', TERM_B = 'SA', REL2A = F, REL2B = T) %>% 
                        dplyr::filter(Term_B %in% c('LE', 'Tumor')) %>% 
                        dplyr::filter(Term_A == 'Epi_p53') %>% 
                        dplyr::mutate(Value = ifelse(Prop4Filter >= 0.01, Value, 0)) %>% 
                        dplyr::left_join(ClinicalInfo_Pat_Binary) %>% 
                        dplyr::mutate(T7th = ifelse(T7th == 0, 'T2-T3', 'T4') %>% as.factor()) %>% 
                        dplyr::mutate(Term_B = ifelse(Term_B == 'LE', paste0(Term_A, ' on LE'), paste0(Term_A, ' on Tumor')) %>% as.factor()) %>% 
                        ggplot(aes(x = T7th, y = Value)) + 
                        facet_wrap(. ~ Term_B, scales = 'free_y') +
                        geom_boxplot(aes(fill = Term_B)) +
                        ggprism::theme_prism(border = T, base_size = 10, base_fontface = 'plain') +
                        scale_fill_manual(values = c('#6388b8', '#a5c561'), guide = 'none') + 
                        theme(aspect.ratio = 2.5, strip.text = element_text(size = 8)) +
                        ggpubr::stat_compare_means(method = 'wilcox', size = 3) + 
                        labs(y = NULL)
                      )
dev.off()


# Figure S6e
pdf('FigS6e.pdf', height = 4.5, width = 6)
MT_CT %>% 
  dplyr::filter(CT0 == 'TME') %>% 
  dplyr::filter(!is.na(Layer)) %>% 
  dplyr::mutate(Layer = Layer %>% as.factor()) %>% 
  dplyr::count(Layer, CT3) %>% 
  dplyr::group_by(Layer) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::ungroup() %>% 
  ggplot() +
  geom_bar(aes(x = Layer, y = prop, fill = CT3), color = 'white', 
           linewidth = 0.05,
           stat = 'identity', position = 'stack') +
  ggprism::theme_prism(border = T, base_size = 12, base_fontface = 'plain') +
  theme(aspect.ratio = 1) +
  scale_x_discrete(expand = expansion(mult = c(0.01, 0.01)),
                   breaks = c('-10', '-5', '0', '5', '10')) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01)), 
                     breaks = c(0, 0.5, 1), 
                     labels = scales::percent_format()) +
  scale_fill_manual(values = c('#3B597AFF', 
                               paletteer::paletteer_d("Redmonder::sPBIGn", n = 3)[-1],
                               paletteer::paletteer_d("Redmonder::sPBIYl", n = 6)[-1],
                               paletteer::paletteer_d("Redmonder::sPBIRd", n = 8)[-1],
                               'grey70') %>% 
                      setNames(LEVELS$CT3)) +
  labs(y = NULL)
dev.off()


# Figure S6f
temp <- MT_CT %>% 
  dplyr::filter(CT0 == 'TME') %>%
  dplyr::filter(!is.na(Layer)) %>% 
  dplyr::transmute(Layer = Layer %>% forcats::as_factor(), 
                   CT3 = CT3 %>% forcats::fct_drop(),
                   CT4 = CT4 %>% forcats::fct_drop())
OR_result2 <- as.factor(-10:10) %>% 
  purrr::map(\(X) {
    temp %>% 
      dplyr::mutate(Layer = ifelse(Layer == X, 'Y', 'N') %>% 
                      forcats::fct_expand(c('Y', 'N')) %>% 
                      forcats::fct_relevel(c('Y', 'N'))) %>% 
      dplyr::select(CT3, Layer) %>% 
      table() %>% 
      rstatix::row_wise_fisher_test(detailed = T) %>% 
      dplyr::transmute(CT = group,
                       Layer = X,
                       OR = estimate)
  })

pdf('FigS6f.pdf', height = 6, width = 8)
OR_result2 %>% 
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
                  forcats::fct_relevel(c('PC', 'Granulo', 'NK', 'DC', 'TME_unknown', 'Mono', 'M1', 'CAF', 
                                         'M2', 'Endo', 'DNT', 'CD4_Treg', 'CD8_T', 'BC', 'CD4_Tconv'))) %>% 
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
