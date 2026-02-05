setwd('~/projects/NPC_SpatialProteomics/')
library(tidyverse)
library(data.table)
library(tidygraph)
library(ggraph)

MT_CT <- qs::qread('SingleCellDataFrame.qs')
LEVELS <- qs::qread('AnnotationLevels.qs')
ClinicalInfo_ROI <- data.table::fread('ClinicalInfo.csv')
ClinicalInfo_Pat <- ClinicalInfo_ROI %>% 
  dplyr::select(-ROI_ID, -TLS_ROI) %>% 
  dplyr::distinct()

TLSProfile_ROI <- ClinicalInfo_ROI %>% 
  dplyr::select(Pat_ID, ROI_ID, TLS_Pat, TLS_ROI) %>% 
  dplyr::mutate(TLS_group = dplyr::case_when(
    TLS_ROI == 0 & TLS_Pat == 0 ~ 'TLS_negative',
    TLS_ROI == 0 & TLS_Pat > 0 ~ 'TLS_adjacent',
    T ~ 'TLS_positive'
  ) %>% forcats::fct_relevel('TLS_negative', 'TLS_adjacent', 'TLS_positive'))

Net_Result <- qs::qread('Net_Result.qs')


# Figure S10a
pdf('FigS10a.pdf', height = 3, width = 5)
Net_Result %>% 
  dplyr::transmute(ROI_ID,
                   Lymphoid_Metrics = GRAPH_METRICS %>% 
                     tidygraph::activate(nodes) %>% 
                     dplyr::as_tibble() %>% 
                     list()) %>% 
  tidyr::unnest() %>% 
  dplyr::filter(CT2 == 'Lymphoid')  %>% 
  dplyr::mutate(Reticulation = ifelse(n_group_size >= 10, 'Reticulated', 'Sparse') %>% forcats::fct_relevel(c('Sparse', 'Reticulated'))) %>% 
  dplyr::mutate(CT = ifelse(CT3 %in% LEVELS$CT3[9:12], 'TC', as.character(CT3))) %>% 
  dplyr::count(CT, Reticulation) %>% 
  dplyr::group_by(CT) %>% 
  dplyr::mutate(n_sum = sum(n), 
                prop = n / sum(n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(Reticulation), prop) %>% 
  dplyr::mutate(CT = CT %>% forcats::fct_inorder()) %>% 
  {list(
    ggplot(., aes(y = CT, fill = Reticulation, x = prop)) + 
      geom_bar(stat = 'identity', position = 'stack') + 
      geom_text(aes(label = ifelse(Reticulation == 'Reticulated', paste0(round(prop * 100), '%'), NA)), 
                size = 2.5, hjust = -0.2) + 
      ggprism::theme_prism(base_size = 10, base_fontface = 'plain', border = T) +
      theme(strip.text = element_text(size = 12), 
            panel.spacing = unit(1, 'lines')) +
      scale_fill_manual(values = c('#DDDDDD', '#885093')) + 
      scale_x_continuous(breaks = c(0, 0.5, 1),
                         labels = c('0', '50%', '100%')
      ) +
      labs(x = 'Proportion of Reticulation', y = NULL) +
      coord_cartesian(clip = 'off'),
    ggplot(data = dplyr::distinct(., CT, Reticulation, n),
           aes(y = CT, fill = Reticulation, x = n)) + 
      geom_bar(stat = 'identity', position = 'stack') + 
      ggprism::theme_prism(base_size = 10, base_fontface = 'plain') +
      theme(strip.text = element_text(size = 12), 
            panel.spacing = unit(1, 'lines'), 
            legend.position = 'none', 
            axis.text.y = element_blank(), 
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.x = element_line(color = '#00000070', linewidth = 0.25)) +
      scale_fill_manual(values = c('#DDDDDD', '#885093'), guide = 'none') + 
      scale_x_continuous(breaks = c(0, 1.5e5, 3e5),
                         labels = ~ format(.x, big.mark = ",", scientific = F),
                         expand = expansion(c(0, 0))
      ) +
      labs(x = 'Cell Count', y = NULL) +
      coord_cartesian(clip = 'off')
  )} %>% 
  patchwork::wrap_plots(nrow = 1, widths = c(1, 1), guides = 'collect')
dev.off()


# Figure S10b
pdf('FigS10b.pdf', width = 12, height = 7)
Net_Result %>%
  dplyr::mutate(GRAPH_METRICS = GRAPH_METRICS %>%
                  tidygraph::activate(edges) %>% 
                  dplyr::mutate(n_group_id = tidygraph::.N()$n_group_id[from]) %>%
                  dplyr::as_tibble() %>% 
                  list(),
                NROW = nrow(GRAPH_METRICS)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(NROW > 0) %>% 
  dplyr::select(ROI_ID, GRAPH_METRICS) %>% 
  tidyr::unnest(GRAPH_METRICS) %>% 
  dplyr::group_by(ROI_ID, n_group_id, n_group_size, e_pair_CT4) %>% 
  dplyr::summarize(MeanBetweenness = exp(mean(log(e_c_betweenness))),
                   nEdge = dplyr::n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(ROI_ID, n_group_id, n_group_size) %>% 
  dplyr::mutate(FreqEdge = nEdge / sum(nEdge),
                SumEdge = sum(nEdge)) %>% 
  dplyr::filter(FreqEdge >= 0.1) %>%
  dplyr::slice_max(order_by = MeanBetweenness, n = 1, with_ties = F) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(e_pair_CT4) %>% 
  dplyr::summarize(MeanSize = mean(n_group_size),
                   NetCount = dplyr::n()) %>%
  dplyr::ungroup() %>% 
  tidyr::separate_wider_delim(cols = 'e_pair_CT4', delim = ';', names =  c('pair_1', 'pair_2')) %>% 
  dplyr::mutate(pair_1 = pair_1 %>% forcats::fct_relevel(LEVELS$CT4[17:34]),
                pair_2 = pair_2 %>% forcats::fct_relevel(LEVELS$CT4[17:34]) %>% forcats::fct_rev()) %>% 
  ggplot(aes(x = pair_1, y = pair_2)) + 
  geom_tile(aes(fill = NetCount)) +
  geom_text(aes(label = NetCount)) +
  scale_fill_gradientn(colours = c('#5D74A5FF', '#B0CBE7FF', '#FEF7C7FF', '#EBA07EFF', '#A8554EFF'),
                       limits = c(0, 2000)) +
  coord_fixed(ratio = 0.6) +
  ggprism::theme_prism(base_fontface = 'plain') + 
  theme(legend.title = element_text(), 
        legend.position = c(0.8, 0.8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = NULL, y = NULL)
dev.off()


# Figure S10c
pdf('FigS10c.pdf', width = 12, height = 7)
Net_Result %>%
  dplyr::mutate(GRAPH_METRICS = GRAPH_METRICS %>%
                  tidygraph::activate(edges) %>% 
                  dplyr::mutate(n_group_id = tidygraph::.N()$n_group_id[from]) %>%
                  dplyr::as_tibble() %>% 
                  list(),
                NROW = nrow(GRAPH_METRICS)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(NROW > 0) %>% 
  dplyr::select(ROI_ID, GRAPH_METRICS) %>% 
  tidyr::unnest(GRAPH_METRICS) %>% 
  dplyr::group_by(ROI_ID, n_group_id, n_group_size, e_pair_CT4) %>% 
  dplyr::summarize(MeanBetweenness = exp(mean(log(e_c_betweenness))),
                   nEdge = dplyr::n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(ROI_ID, n_group_id, n_group_size) %>% 
  dplyr::mutate(FreqEdge = nEdge / sum(nEdge),
                SumEdge = sum(nEdge)) %>% 
  dplyr::filter(FreqEdge >= 0.1) %>%
  dplyr::slice_max(order_by = MeanBetweenness, n = 1, with_ties = F) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(e_pair_CT4) %>% 
  dplyr::summarize(MeanSize = mean(n_group_size),
                   NetCount = dplyr::n()) %>%
  dplyr::ungroup() %>% 
  tidyr::separate_wider_delim(cols = 'e_pair_CT4', delim = ';', names =  c('pair_1', 'pair_2')) %>% 
  dplyr::mutate(pair_1 = pair_1 %>% forcats::fct_relevel(LEVELS$CT4[17:34]),
                pair_2 = pair_2 %>% forcats::fct_relevel(LEVELS$CT4[17:34]) %>% forcats::fct_rev()) %>% 
  ggplot(aes(x = pair_1, y = pair_2)) + 
  geom_tile(aes(fill = MeanSize)) +
  geom_text(aes(label = sprintf('%.1f', MeanSize))) +
  scale_fill_gradientn(colours = c('#5D74A5FF', '#B0CBE7FF', '#FEF7C7FF', '#EBA07EFF', '#A8554EFF'),
                       limits = c(10, 200)) +
  coord_fixed(ratio = 0.6) +
  ggprism::theme_prism(base_fontface = 'plain') + 
  theme(legend.title = element_text(), 
        legend.position = c(0.8, 0.8),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = NULL, y = NULL)
dev.off()


# Figure S10d
pdf('FigS10d.pdf', width = 12, height = 6)
Net_Result %>% 
  dplyr::filter(ROI_ID %in% c('T052_ROI4', 'T161_ROI1')) %>% 
  dplyr::as_tibble() %>% 
  dplyr::rowwise() %>% 
  dplyr::transmute(ROI_ID, 
                   NODES = GRAPH_METRICS %>% 
                     tidygraph::activate(nodes) %>%
                     data.frame() %>% 
                     list(),
                   EDGES = GRAPH_METRICS %>% 
                     tidygraph::activate(edges) %>%
                     dplyr::mutate(n_group_c_id = tidygraph::.N()$n_group_c_id[from]) %>% 
                     data.frame() %>% 
                     list()) %>% 
  dplyr::mutate(NODES = NODES %>% 
                  dplyr::mutate(CT = dplyr::case_when(
                    CT4 %in% c('CD4_T_CD45RO', 'CD8_T_GranzymeB_2') ~ as.character(CT4),
                    T ~ 'Others')) %>% 
                  list(),
                EDGES = EDGES %>% 
                  dplyr::filter(n_group_c_id == 1) %>% 
                  dplyr::arrange(e_c_betweenness) %>% 
                  list()) %>% 
  dplyr::mutate(PLOT = {ggplot() + 
      geom_segment(aes(x = x1, xend = x2, 
                       y = -y1, yend = -y2, 
                       color = log10(e_c_betweenness)), 
                   linewidth = 3,
                   linejoin = "round", lineend = "round",
                   data = EDGES) + 
      scale_color_gradientn(name = 'Edge Centrality', 
                            colors = c('#5D74A5FF', '#B0CBE7FF', '#FEF7C7FF', '#EBA07EFF', '#A8554EFF'),
                            labels = function(x) parse(text = gsub("e[+]", "%*%10^", scales::scientific_format()(10^x)))) +
      geom_segment(aes(x = x1, xend = x2, y = -y1, yend = -y2),
                   color = '#AAAAAA', linewidth = 0.2,
                   data = EDGES) +
      geom_point(aes(x = X_position, y = -Y_position, fill = !!as.name('CT')), 
                 size = 1.5, shape = 21, stroke = 0, color = '#000000', 
                 data = NODES) +
      scale_fill_manual(values = c('#1270AE', '#22ADA5', 'grey90') %>% 
                          setNames(c('CD4_T_CD45RO', 'CD8_T_GranzymeB_2', 'Others')), 
                        na.value = 'grey90') +
      scale_x_continuous(expand = expansion(mult = c(0, 0), add = c(0, 0))) +
      scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0, 0))) +
      scale_linewidth_identity() +
      theme_void() +
      labs(y = paste0(ROI_ID, '\n')) +
      coord_equal(clip = 'off') +
      theme(axis.title.y = element_text(angle = 90), 
            plot.margin = unit(c(1,1,1,1), 'line'))} %>% 
        list()) %>% 
  dplyr::pull(PLOT) %>% 
  purrr::walk(print)
dev.off()


# Figure S10e
Net_Result_Nodes <- Net_Result %>%
  dplyr::filter((GRAPH_METRICS %>% tidygraph::activate(edges) %>% dplyr::as_tibble() %>% nrow()) > 0) %>% 
  dplyr::mutate(DATA = GRAPH_METRICS %>% 
                  tidygraph::activate(nodes) %>%
                  dplyr::as_tibble() %>% 
                  dplyr::filter(n_group_size >= 10) %>%
                  list()) %>% 
  dplyr::select(ROI_ID, DATA) %>% 
  tidyr::unnest()

pdf('FigS10e.pdf', width = 10, height = 5)
Net_Result_Nodes %>% 
  dplyr::mutate(Net_ID = paste0(ROI_ID, '_', sprintf('%02.f', n_group_id))) %>% 
  dplyr::relocate(Pat_ID, ROI_ID, Net_ID, n_group_size) %>% 
  dplyr::count(Net_ID, CT4) %>% 
  dplyr::mutate(CT4 = CT4 %>% forcats::fct_drop()) %>% 
  tidyr::complete(Net_ID, CT4, fill = list(n = 0)) %>% 
  dplyr::group_by(Net_ID) %>% 
  dplyr::mutate(Size = sum(n),
                logSize = log10(sum(n)), 
                Prop = n / sum(n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::nest_by(CT4) %>%
  tibble::rowid_to_column('Row_ID') %T>% {temp_nrow <<- nrow(.)} %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Model = mgcv::gam(data = data, 
                                  formula = Prop ~ s(logSize, bs = "cs", k = 3), 
                                  family = stats::quasipoisson(),
                                  fitted = T) %>% 
                  list()) %>%
  dplyr::mutate(Prediction = Model %>% predict(type = "response") %>% list()) %>% 
  dplyr::select(-Model) %>% 
  tidyr::unnest(cols = c(data, Prediction)) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(CT4 %in% LEVELS$CT4[17:23]) %>% 
  {ggplot(., aes(x = logSize, y = Prediction, group = CT4, color = CT4)) + 
      geom_line(linewidth = 1) +
      geom_text(aes(x = logSize.x,
                    y = Prediction.x,
                    label = paste0(CT4, ': ', sprintf('%.2f', FC))),
                hjust = 0, 
                data = 
                  dplyr::left_join(dplyr::slice_max(dplyr::group_by(., CT4), order_by = logSize, n = 1, with_ties = F),
                                   dplyr::slice_min(dplyr::group_by(., CT4), order_by = logSize, n = 1, with_ties = F),
                                   by = 'CT4') %>% 
                  dplyr::ungroup() %>% 
                  dplyr::mutate(FC = Prediction.x / Prediction.y)) +
      ggprism::theme_prism(base_size = 10, base_fontface = 'plain', border = T) +
      theme(legend.position = 'none',
            aspect.ratio = 1.2,
            panel.grid = element_line(color = 'grey70', linetype = 'dashed', linewidth = 0.2)) +
      scale_x_continuous(expand = expansion(), breaks = c(1, 2, 3), labels = c(10, 100, 1000)) +
      scale_y_continuous(expand = expansion(), limits = c(0, 0.1), breaks = 0:4 / 40) +
      scale_color_manual(name = 'CT',
                         values = c('#43978FFF', '#56A09FFF', '#68AAAFFF', '#79B3BFFF', '#8ABDCFFF', '#9AC6E0FF',
                                    '#ABD0F1FF', '#C5B9CBFF', '#D5A1A5FF', '#E08981FF', '#E56F5EFF', '#EA875DFF', 
                                    '#EF9D5CFF', '#F3B35AFF') %>% 
                           setNames(LEVELS$CT4[17:30]),
                         na.value = 'grey80', drop = F) +
      labs(x = 'Net Size (in log scale)', y = 'Estimated Prop')+ 
      coord_fixed(clip = 'off')}
Net_Result_Nodes %>% 
  dplyr::mutate(Net_ID = paste0(ROI_ID, '_', sprintf('%02.f', n_group_id))) %>% 
  dplyr::relocate(Pat_ID, ROI_ID, Net_ID, n_group_size) %>% 
  dplyr::count(Net_ID, CT4) %>% 
  dplyr::mutate(CT4 = CT4 %>% forcats::fct_drop()) %>% 
  tidyr::complete(Net_ID, CT4, fill = list(n = 0)) %>% 
  dplyr::group_by(Net_ID) %>% 
  dplyr::mutate(Size = sum(n),
                logSize = log10(sum(n)), 
                Prop = n / sum(n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::nest_by(CT4) %>%
  tibble::rowid_to_column('Row_ID') %T>% {temp_nrow <<- nrow(.)} %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Model = mgcv::gam(data = data, 
                                  formula = Prop ~ s(logSize, bs = "cs", k = 3), 
                                  family = stats::quasipoisson(),
                                  fitted = T) %>% 
                  list()) %>%
  dplyr::mutate(Prediction = Model %>% predict(type = "response") %>% list()) %>% 
  dplyr::select(-Model) %>% 
  tidyr::unnest(cols = c(data, Prediction)) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(CT4 %in% LEVELS$CT4[24:30]) %>%
  {ggplot(., aes(x = logSize, y = Prediction, group = CT4, color = CT4)) + 
      geom_line(linewidth = 1) +
      geom_text(aes(x = logSize.x,
                    y = Prediction.x,
                    label = paste0(CT4, ': ', sprintf('%.2f', FC))),
                hjust = 0, 
                data = 
                  dplyr::left_join(dplyr::slice_max(dplyr::group_by(., CT4), order_by = logSize, n = 1, with_ties = F),
                                   dplyr::slice_min(dplyr::group_by(., CT4), order_by = logSize, n = 1, with_ties = F),
                                   by = 'CT4') %>% 
                  dplyr::ungroup() %>% 
                  dplyr::mutate(FC = Prediction.x / Prediction.y)) +
      ggprism::theme_prism(base_size = 10, base_fontface = 'plain', border = T) +
      theme(legend.position = 'none',
            aspect.ratio = 1.2,
            panel.grid = element_line(color = 'grey70', linetype = 'dashed', linewidth = 0.2)) +
      scale_x_continuous(expand = expansion(), breaks = c(1, 2, 3), labels = c(10, 100, 1000)) +
      scale_y_continuous(expand = expansion(), limits = c(0, 0.1), breaks = 0:4 / 40) +
      scale_color_manual(name = 'CT',
                         values = c('#43978FFF', '#56A09FFF', '#68AAAFFF', '#79B3BFFF', '#8ABDCFFF', '#9AC6E0FF',
                                    '#ABD0F1FF', '#C5B9CBFF', '#D5A1A5FF', '#E08981FF', '#E56F5EFF', '#EA875DFF', 
                                    '#EF9D5CFF', '#F3B35AFF') %>% 
                           setNames(LEVELS$CT4[17:30]),
                         na.value = 'grey80', drop = F) +
      labs(x = 'Net Size (in log scale)', y = 'Estimated Prop')+ 
      coord_fixed(clip = 'off')}
dev.off()


# Figure S10f
NetSynora_Result <- qs::qread('NetSynora_Result.qs')
pdf('FigS10f.pdf', width = 8, height = 7)
MT_CT %>% 
  dplyr::left_join(NetSynora_Result) %>%
  dplyr::filter(CT2 == 'Lymphoid') %>% 
  dplyr::mutate(NetLayer = ifelse(Distance2NetBoundary >= 0, 
                                  pmin(Distance2NetBoundary %/% 10, 10), 
                                  pmax(Distance2NetBoundary %/% 10 + 1, -10))) %>% 
  dplyr::mutate(CT4 = CT4 %>% forcats::fct_drop()) %>% 
  dplyr::count(CT4, NetLayer) %>% 
  dplyr::group_by(NetLayer) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = NetLayer, y = prop, fill = CT4)) + 
  geom_bar(stat = 'identity') + 
  ggprism::theme_prism(border = T, base_size = 12, base_fontface = 'plain') +
  theme(aspect.ratio = 1) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_manual(name = 'CT',
                    values = c('#43978FFF', '#56A09FFF', '#68AAAFFF', '#79B3BFFF', '#8ABDCFFF', '#9AC6E0FF',
                               '#ABD0F1FF', '#C5B9CBFF', '#D5A1A5FF', '#E08981FF', '#E56F5EFF', '#EA875DFF', 
                               '#EF9D5CFF', '#F3B35AFF', '#F6C957', '#F19685', '#FFB77F', '#C59D94') %>% 
                      setNames(LEVELS$CT4[17:34]),
                    na.value = 'grey80', drop = F)
dev.off()
