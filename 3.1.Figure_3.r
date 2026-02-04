setwd('~/projects/NPC_SpatialProteomics/')
library(tidyverse)
library(data.table)
library(tidygraph)
library(ggraph)


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

TLSProfile_ROI <- ClinicalInfo_ROI %>% 
  dplyr::select(Pat_ID, ROI_ID, TLS_Pat, TLS_ROI) %>% 
  dplyr::mutate(TLS_group = dplyr::case_when(
    TLS_ROI == 0 & TLS_Pat == 0 ~ 'TLS_negative',
    TLS_ROI == 0 & TLS_Pat > 0 ~ 'TLS_adjacent',
    T ~ 'TLS_positive'
  ) %>% forcats::fct_relevel('TLS_negative', 'TLS_adjacent', 'TLS_positive'))

Net_Result <- qs::qread('Net_Result.qs')





# Figure 3a
pdf('Fig3a.pdf', width = 12, height = 6)
Net_Result %>% 
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
  dplyr::filter(ROI_ID %in% c("T057_ROI3", "T045_ROI4")) %>% 
  dplyr::mutate(PLOT = {ggplot() + 
      geom_segment(aes(x = x1, xend = x2, 
                       y = -y1, yend = -y2, 
                       color = as.factor(n_group_c_id), 
                       linewidth = log10(e_c_betweenness) + 1), 
                   linejoin = "round", lineend = "round",
                   data = EDGES) + 
      scale_color_manual(values = paletteer::paletteer_d("colorBlindness::SteppedSequential5Steps") %>% .[-(1:4)] %>% rep(2), guide = 'none') +
      geom_segment(aes(x = x1, xend = x2, y = -y1, yend = -y2),
                   color = '#222222', linewidth = 0.2,
                   data = EDGES) +
      geom_point(aes(x = X_position, y = -Y_position, fill = !!as.name('CT3')), 
                 size = 1.5, shape = 21, stroke = 0, color = '#000000', 
                 data = NODES) +
      scale_fill_manual(values = c('#43978F', '#ABD0F1', '#E56F5E', '#F6C957', '#F19685', '#FFB77F', '#C59D94') %>% 
                          setNames(LEVELS$CT3[9:15]), 
                        na.value = 'grey90') +
      scale_x_continuous(expand = expansion(mult = c(0, 0), add = c(0, 0))) +
      scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0, 0))) +
      scale_linewidth_identity() +
      theme_void() +
      coord_equal(clip = 'off') +
      labs(y = paste0(ROI_ID, '\n')) +
      theme(axis.title.y = element_text(angle = 90), 
            plot.margin = unit(c(1,1,1,1), 'line'))} %>% 
        list()) %>% 
  dplyr::pull(PLOT) %>% 
  purrr::walk(print)
dev.off()


# Figure 3b
AverageLikelihood <- Net_Result %>% 
  dplyr::transmute(ROI_ID,
                   Lymphoid_Metrics = GRAPH_METRICS %>% 
                     tidygraph::activate(nodes) %>% 
                     dplyr::as_tibble() %>% 
                     list()) %>% 
  tidyr::unnest() %>% 
  dplyr::filter(CT2 == 'Lymphoid') %>% 
  dplyr::summarize(AverageLikelihood = mean(n_group_size >= 10)) %>% 
  dplyr::pull(AverageLikelihood)

pdf('Fig3b.pdf', height = 6, width = 7)
Net_Result %>% 
  dplyr::transmute(ROI_ID,
                   Lymphoid_Metrics = GRAPH_METRICS %>% 
                     tidygraph::activate(nodes) %>% 
                     dplyr::as_tibble() %>% 
                     list()) %>% 
  tidyr::unnest() %>% 
  dplyr::filter(CT2 == 'Lymphoid')  %>% 
  dplyr::mutate(Reticulation = ifelse(n_group_size >= 10, 'Reticulated', 'Sparse') %>% forcats::fct_relevel(c('Sparse', 'Reticulated'))) %>% 
  dplyr::count(CT4, Reticulation) %>% 
  dplyr::group_by(CT4) %>% 
  dplyr::mutate(n_sum = sum(n), 
                prop = n / sum(n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(Reticulation), prop) %>% 
  dplyr::mutate(CT4 = CT4 %>% forcats::fct_inorder()) %>% 
  {list(
    ggplot(., aes(y = CT4, fill = Reticulation, x = prop)) + 
      geom_bar(stat = 'identity', position = 'stack') + 
      geom_vline(xintercept = AverageLikelihood, linetype = 'dotted', color = 'red', linewidth = 1) + 
      geom_text(aes(label = ifelse(Reticulation == 'Reticulated', paste0(round(prop * 100), '%'), NA)), 
                size = 3, hjust = -0.2) + 
      ggprism::theme_prism(base_size = 10, base_fontface = 'plain', border = T) +
      theme(strip.text = element_text(size = 12), 
            panel.spacing = unit(1, 'lines'), 
            plot.margin = unit(c(1,1,0,1), "cm")
      ) +
      scale_fill_manual(values = c('#DDDDDD', '#885093')) + 
      scale_x_continuous(breaks = c(0, 0.5, AverageLikelihood, 1),
                         labels = c('0', '50%', '', '100%')
      ) +
      annotate(geom = 'text', label = paste0('Average Proportion ', round(100 * AverageLikelihood), '%'), 
               x = 0.65, y = Inf, size = 4, hjust = 0.5, vjust = -1, color = 'red') + 
      labs(x = 'Proportion of Reticulation', y = NULL) +
      coord_cartesian(clip = 'off'),
    ggplot(data = dplyr::distinct(., CT4, Reticulation, n),
             aes(y = CT4, fill = Reticulation, x = n)) + 
      geom_bar(stat = 'identity', position = 'stack') + 
      ggprism::theme_prism(base_size = 10, base_fontface = 'plain', border = F) +
      theme(strip.text = element_text(size = 12), 
            panel.spacing = unit(1, 'lines'), 
            legend.position = 'none', 
            axis.text.y = element_blank(), 
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.x = element_line(color = '#00000070', linewidth = 0.25)) +
      scale_fill_manual(values = c('#DDDDDD', '#885093')) + 
      scale_x_continuous(
        breaks = c(0, 100000, 200000),
        labels = ~ format(.x, big.mark = ",", scientific = F),
        expand = expansion(c(0, 0))
      ) +
      labs(x = 'Cell Count', y = NULL) +
      coord_cartesian(clip = 'off')
  )} %>% 
  patchwork::wrap_plots(nrow = 1, widths = c(1, 1.5), guides = 'collect')
dev.off()





# Figure 3c
pdf('Fig3c.pdf', width = 12, height = 7)
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
  dplyr::group_by(ROI_ID, n_group_id, n_group_size, e_pair_CT3) %>% 
  dplyr::summarize(MeanBetweenness = exp(mean(log(e_c_betweenness))),
                   nEdge = dplyr::n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(ROI_ID, n_group_id, n_group_size) %>% 
  dplyr::mutate(FreqEdge = nEdge / sum(nEdge),
                SumEdge = sum(nEdge)) %>% 
  dplyr::filter(FreqEdge >= 0.1) %>%
  dplyr::slice_max(order_by = MeanBetweenness, n = 1, with_ties = F) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(e_pair_CT3) %>% 
  dplyr::summarize(NetCount = dplyr::n()) %>%
  dplyr::ungroup() %>% 
  tidyr::separate_wider_delim(cols = 'e_pair_CT3', delim = ';', names =  c('pair_1', 'pair_2')) %>% 
  dplyr::mutate(pair_1 = pair_1 %>% forcats::fct_relevel(LEVELS$CT3[9:15]),
                pair_2 = pair_2 %>% forcats::fct_relevel(LEVELS$CT3[9:15])) %>% 
  tidyr::pivot_wider(names_from = 'pair_1', values_from = 'NetCount', names_sort = T) %>% 
  dplyr::arrange(pair_2) %>% 
  tibble::column_to_rownames('pair_2') %>% 
  as.matrix() %>% 
  tidygraph::as_tbl_graph() %>% 
  tidygraph::activate(nodes) %>% 
  dplyr::mutate(name = name %>% forcats::fct_inorder()) %>% 
  tidygraph::activate(edges) %>%
  dplyr::arrange(weight) %>%
  ggraph(layout = 'circle') +
  geom_edge_fan2(aes(width = weight, alpha = weight, edge_colour = weight), 
                 label_dodge = T, alpha = 1, lineend = 'round') +
  geom_edge_loop(aes(width = weight, alpha = weight, edge_colour = weight, direction = (from - 1) * 360 / length(LEVELS$CT3[9:15])), 
                 label_dodge = T, alpha = 1, lineend = 'round') +
  geom_node_point(aes(color = name), size = 6) +
  theme_graph(base_family = 'Helvetica') +
  coord_equal(clip = 'off') +
  scale_color_manual(name = 'CT',
                     values = c('#43978F', '#ABD0F1', '#E56F5E', '#F6C957', '#F19685', '#FFB77F', '#C59D94') %>% 
                       setNames(LEVELS$CT3[9:15]),
                     na.value = 'grey80', drop = F) +
  scale_edge_width_continuous(name = 'NetCount', 
                              trans = 'log10', range = c(0.1, 3)) +
  scale_edge_alpha_continuous(trans = 'log10') +
  scale_edge_color_gradientn(name = 'NetCount', 
                             trans = 'log10', 
                             colors = paletteer::paletteer_c('grDevices::OrRd', n = 100)[10:100] %>% rev())
dev.off()


# Figure 3d
pdf('Fig3d.pdf', width = 12, height = 7)
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
  dplyr::group_by(ROI_ID, n_group_id, n_group_size, e_pair_CT3) %>% 
  dplyr::summarize(MeanBetweenness = exp(mean(log(e_c_betweenness))),
                   nEdge = dplyr::n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(ROI_ID, n_group_id, n_group_size) %>% 
  dplyr::mutate(FreqEdge = nEdge / sum(nEdge),
                SumEdge = sum(nEdge)) %>% 
  dplyr::filter(FreqEdge >= 0.1) %>%
  dplyr::slice_max(order_by = MeanBetweenness, n = 1, with_ties = F) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(e_pair_CT3) %>% 
  dplyr::summarize(AverageNetSize = mean(n_group_size)) %>%
  dplyr::ungroup() %>% 
  tidyr::separate_wider_delim(cols = 'e_pair_CT3', delim = ';', names =  c('pair_1', 'pair_2')) %>% 
  dplyr::mutate(pair_1 = pair_1 %>% forcats::fct_relevel(LEVELS$CT3[9:15]),
                pair_2 = pair_2 %>% forcats::fct_relevel(LEVELS$CT3[9:15])) %>% 
  tidyr::pivot_wider(names_from = 'pair_1', values_from = 'AverageNetSize', names_sort = T) %>% 
  dplyr::arrange(pair_2) %>% 
  tibble::column_to_rownames('pair_2') %>% 
  as.matrix() %>% 
  tidygraph::as_tbl_graph() %>% 
  tidygraph::activate(nodes) %>% 
  dplyr::mutate(name = name %>% forcats::fct_inorder()) %>% 
  tidygraph::activate(edges) %>%
  dplyr::arrange(weight) %>%
  ggraph(layout = 'circle') +
  geom_edge_fan2(aes(width = weight, alpha = weight, edge_colour = weight), 
                 label_dodge = T, alpha = 1, lineend = 'round') +
  geom_edge_loop(aes(width = weight, alpha = weight, edge_colour = weight, direction = (from - 1) * 360 / length(LEVELS$CT3[9:15])), 
                 label_dodge = T, alpha = 1, lineend = 'round') +
  geom_node_point(aes(color = name), size = 6) +
  theme_graph(base_family = 'Helvetica') +
  coord_equal(clip = 'off') +
  scale_color_manual(name = 'CT',
                     values = c('#43978F', '#ABD0F1', '#E56F5E', '#F6C957', '#F19685', '#FFB77F', '#C59D94') %>% 
                       setNames(LEVELS$CT3[9:15]),
                     na.value = 'grey80', drop = F) +
  scale_edge_width_continuous(name = 'AverageNetSize', 
                              trans = 'log10', 
                              range = c(0.1, 3)) +
  scale_edge_alpha_continuous(trans = 'log10') +
  scale_edge_color_gradientn(name = 'AverageNetSize', 
                             trans = 'log10', 
                             colors = paletteer::paletteer_c('grDevices::PuBu', n = 100)[10:100] %>% rev())
dev.off()


# Figure 3e
Net_Result_Nodes <- Net_Result %>%
  dplyr::filter((GRAPH_METRICS %>% tidygraph::activate(edges) %>% dplyr::as_tibble() %>% nrow()) > 0) %>% 
  dplyr::mutate(DATA = GRAPH_METRICS %>% 
                  tidygraph::activate(nodes) %>%
                  dplyr::as_tibble() %>% 
                  dplyr::filter(n_group_size >= 10) %>%
                  list()) %>% 
  dplyr::select(ROI_ID, DATA) %>% 
  tidyr::unnest()

pdf('Fig3e.pdf', width = 6, height = 5)
Net_Result_Nodes %>% 
  dplyr::mutate(Net_ID = paste0(ROI_ID, '_', sprintf('%02.f', n_group_id))) %>% 
  dplyr::relocate(Pat_ID, ROI_ID, Net_ID, n_group_size) %>% 
  dplyr::count(Net_ID, CT3) %>% 
  dplyr::mutate(CT3 = CT3 %>% forcats::fct_drop()) %>% 
  tidyr::complete(Net_ID, CT3, fill = list(n = 0)) %>% 
  dplyr::group_by(Net_ID) %>% 
  dplyr::mutate(Size = sum(n),
                logSize = log10(sum(n)), 
                Prop = n / sum(n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::nest_by(CT3) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Model = mgcv::gam(data = data, 
                                  formula = Prop ~ s(logSize, bs = "cs", k = 3), 
                                  family = stats::quasipoisson(),
                                  fitted = T) %>% 
                  list()) %>%
  dplyr::mutate(Prediction = Model %>% predict(type = "response") %>% list()) %>% 
  dplyr::select(CT3, data, Prediction) %>% 
  tidyr::unnest(cols = c(data, Prediction)) %>%
  dplyr::ungroup() %>% 
  {ggplot(., aes(x = logSize, y = Prediction, color = CT3)) + 
      geom_line(linewidth = 1) +
      geom_text(aes(x = logSize.x,
                    y = Prediction.x,
                    label = paste0(CT3, ': ', sprintf('%.2f', FC))),
                hjust = 0, 
                data = 
                  dplyr::left_join(dplyr::slice_max(dplyr::group_by(., CT3), order_by = logSize, n = 1, with_ties = F),
                                   dplyr::slice_min(dplyr::group_by(., CT3), order_by = logSize, n = 1, with_ties = F),
                                   by = 'CT3') %>% 
                  dplyr::ungroup() %>% 
                  dplyr::mutate(FC = Prediction.x / Prediction.y)) +
      ggprism::theme_prism(base_size = 10, base_fontface = 'plain', border = T) +
      theme(legend.position = 'none', 
            aspect.ratio = 1.2,
            panel.grid = element_line(color = 'grey70', linetype = 'dashed', linewidth = 0.2)) +
      scale_x_continuous(expand = expansion(), breaks = c(1, 2, 3), labels = c(10, 100, 1000)) +
      scale_y_continuous(expand = expansion(), limits = c(0, 0.55), breaks = 1:5 / 10) +
      scale_color_manual(values = c('#43978F', '#ABD0F1', '#E56F5E', '#F6C957', '#F19685', '#FFB77F', '#C59D94') %>% 
                           setNames(LEVELS$CT3[9:15])) + 
      labs(x = 'Net Size (in log scale)', y = 'Estimated Prop', title = '') + 
      coord_cartesian(clip = 'off')}

Net_Result_Nodes %>% 
  dplyr::mutate(Net_ID = paste0(ROI_ID, '_', sprintf('%02.f', n_group_id))) %>% 
  dplyr::relocate(Pat_ID, ROI_ID, Net_ID, n_group_size) %>% 
  dplyr::count(Net_ID, CT3) %>% 
  dplyr::mutate(CT3 = CT3 %>% forcats::fct_drop()) %>% 
  tidyr::complete(Net_ID, CT3, fill = list(n = 0)) %>% 
  dplyr::group_by(Net_ID) %>% 
  dplyr::mutate(Size = sum(n),
                logSize = log10(sum(n)), 
                Prop = n / sum(n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::nest_by(CT3) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Model = mgcv::gam(data = data, 
                                  formula = n ~ s(logSize, bs = "cs", k = 3), 
                                  family = stats::quasipoisson(),
                                  fitted = T) %>% 
                  list()) %>%
  dplyr::mutate(Prediction = Model %>% predict(type = "response") %>% list()) %>% 
  dplyr::select(CT3, data, Prediction) %>% 
  tidyr::unnest(cols = c(data, Prediction)) %>%
  dplyr::ungroup() %>% 
  {ggplot(., aes(x = logSize, y = log10(Prediction), color = CT3)) +
      geom_abline(slope = 1, intercept = log10(1), linetype = 'dashed', color = 'grey30') +
      geom_line(linewidth = 1) +
      annotate(geom = 'text', label = 'Prop = 1', x = Inf, y = Inf, hjust = 2, vjust = 2.5) +
      geom_text(aes(x = logSize,
                    y = log10_Prediction,
                    label = CT3),
                hjust = 0, 
                data = dplyr::transmute(., CT3, logSize, log10_Prediction = log10(Prediction)) %>% 
                  dplyr::group_by(CT3) %>%
                  dplyr::slice_max(order_by = logSize, n = 1, with_ties = F)) +
      ggprism::theme_prism(base_size = 10, base_fontface = 'plain', border = T) +
      theme(legend.position = 'none',
            aspect.ratio = 1.2,
            panel.grid = element_line(color = 'grey70', linetype = 'dashed', linewidth = 0.2)) +
      scale_x_continuous(expand = expansion(),
                         breaks = c(1, 2, 3),
                         labels = c(10, 100, 1000)
      ) +
      scale_y_continuous(expand = expansion(),
                         breaks = c(0, 1, 2, 3),
                         labels = c(1, 10, 100, 1000),
                         limits = c(0, 3.6)) +
      scale_color_manual(values = c('#43978F', '#ABD0F1', '#E56F5E', '#F6C957', '#F19685', '#FFB77F', '#C59D94') %>%
                           setNames(LEVELS$CT3[9:15])) +
      labs(x = 'Net Size (in log scale)', y = 'Estimated Count', title = '') +
      coord_fixed(clip = 'off')}
dev.off()






# Figure 3f
NetSynora_Result <- qs::qread('NetSynora_Result.qs')
pdf('Fig3f.pdf', width = 7, height = 7)
MT_CT %>% 
  dplyr::left_join(NetSynora_Result) %>%
  dplyr::filter(CT2 == 'Lymphoid') %>% 
  dplyr::mutate(NetLayer = ifelse(Distance2NetBoundary >= 0, 
                                  pmin(Distance2NetBoundary %/% 10, 10), 
                                  pmax(Distance2NetBoundary %/% 10 + 1, -10))) %>% 
  dplyr::mutate(CT3 = CT3 %>% forcats::fct_drop()) %>% 
  dplyr::count(CT3, NetLayer) %>% 
  dplyr::group_by(NetLayer) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::ungroup() %>% 
  ggplot(aes(x = NetLayer, y = prop, fill = CT3)) + 
  geom_bar(stat = 'identity') + 
  ggprism::theme_prism(border = T, base_size = 12, base_fontface = 'plain') +
  theme(aspect.ratio = 1) +
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_manual(name = 'none',
                    values = c('#43978F', '#ABD0F1', '#E56F5E', '#F6C957', '#F19685', '#FFB77F', '#C59D94') %>% 
                      setNames(LEVELS$CT3[9:15]),
                    na.value = 'grey80', drop = F)
dev.off()


# Figure 3g
pdf('Fig3g.pdf', height = 6, width = 10)
Net_Result_Nodes %>% 
  dplyr::transmute(ROI_ID, 
                   Net_ID = paste0(ROI_ID, '_', sprintf('%02.f', n_group_id)), 
                   n_group_size) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(ROI_ID) %>% 
  dplyr::mutate(MeanSize = mean(n_group_size),
                MaxSize = max(n_group_size), 
                N = dplyr::n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(desc(MaxSize), desc(MeanSize), desc(N), ROI_ID) %>% 
  dplyr::left_join(TLSProfile_ROI %>%
                     dplyr::transmute(ROI_ID, 
                                      TLS_group)) %>% 
  dplyr::mutate(ROI_ID = ROI_ID %>% forcats::fct_inorder()) %>% 
  dplyr::mutate(log10Size = log10(n_group_size)) %>% 
  {list(
    dplyr::distinct(., TLS_group, Net_ID, n_group_size, log10Size) %T>% 
      {temp_label <<- dplyr::group_by(., TLS_group) %>% 
        dplyr::slice_max(order_by = n_group_size) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(Label = paste0('Max Size:\n', n_group_size))} %>% 
      dplyr::mutate(Bin = log10Size %>% Hmisc::cut2(cuts = 20:71/20 - 0.025)) %>% 
      dplyr::count(TLS_group, Bin) %>%
      tidyr::complete(TLS_group, Bin, fill = list(n = 0)) %>% 
      ggplot() + 
      facet_wrap(. ~ TLS_group, scales = 'free_x', labeller = label_parsed) +
      geom_bar(aes(y = Bin, x = n, fill = TLS_group), 
               linewidth = 0.2, stat = 'identity', 
               position = 'stack',
               color = 'black') +
      geom_text(aes(label = Label), 
                x = -Inf, y = Inf, size = 3,
                hjust = 0,
                vjust = 1.5, 
                data = temp_label) +
      ggprism::theme_prism(base_size = 8, base_fontface = 'plain', border = T) +
      scale_y_discrete(expand = expansion(add = c(1, 1)),
                       breaks = c('[0.975,1.025)', '[1.975,2.025)', '[2.975,3.025)'),
                       labels = c('10', '100', '1000')) +
      scale_x_continuous(transform = 'reverse') +
      scale_fill_manual(values = c('#D6C7C2', '#7DA9B8', '#A94D42'), guide = 'none') + 
      labs(y = 'Size (in log scale)', x = NULL) + 
      coord_cartesian(clip = 'off'),
    dplyr::arrange(. ,TLS_group) %>% 
      ggplot() +
      geom_point(aes(x = ROI_ID, y = log10Size, color = TLS_group), size = 0.5) + 
      ggprism::theme_prism(base_size = 10, base_fontface = 'plain', border = T) + 
      scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
      scale_y_continuous(expand = expansion(add = c(0.05, 0.05))) +
      scale_color_manual(values = c('#D6C7C2', '#7DA9B8', '#A94D42')) + 
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(), 
            legend.position = c(0.8, 0.8)),
    patchwork::plot_spacer(),
    ggplot(., aes(x = ROI_ID)) +
      geom_tile(aes(y = 1, fill = TLS_group)) +
      ggprism::theme_prism(base_size = 10, base_fontface = 'plain', border = T) + 
      scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
      # scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()) +
      scale_fill_manual(values = c('#D6C7C2', '#7DA9B8', '#A94D42'), guide = 'none')
  )} %>% 
  patchwork::wrap_plots(nrow = 2,
                        heights = c(1, 0.05),
                        widths = c(0.3, 1))
dev.off()

# Figure 3h
NetProfile_ROI <- fread('NetProfile_ROI.csv')
pdf('Fig3h.pdf', width = 5, height = 6)
NetProfile_ROI %>% 
  tidyr::pivot_longer(cols = !ROI_ID, names_to = 'Feature', values_to = 'Value') %>% 
  dplyr::filter(Feature == 'LSS') %>% 
  dplyr::left_join(TLSProfile_ROI %>%
                     dplyr::transmute(ROI_ID, 
                                      TLS_group)) %>%
  ggplot(aes(x = TLS_group, y = Value)) +
  facet_wrap(. ~ Feature, scales = 'free_y') +
  geom_violin(aes(fill = TLS_group)) +
  geom_boxplot(width = 0.1, outliers = F) +
  ggpubr::stat_compare_means(method = 'wilcox',
                             comparisons = list(
                               c('TLS_negative', 'TLS_adjacent'),
                               c('TLS_adjacent', 'TLS_positive'),
                               c('TLS_negative', 'TLS_positive')), 
  ) +
  scale_fill_manual(values = c('#D6C7C2', '#7DA9B8', '#A94D42'), guide = 'none') + 
  scale_y_continuous(expand = expansion(c(0.01, 0.05))) +
  ggprism::theme_prism(base_size = 12, base_fontface = 'plain', border = T) +
  labs(x = NULL, y = NULL)
dev.off()

# Figure 3i
NetProfile <- fread('NetProfile_Pat.csv')
pdf('Fig3i.pdf', width = 7, height = 6)
NetProfile %>% 
  tidyr::pivot_longer(cols = !Pat_ID, names_to = 'Feature', values_to = 'Value') %>% 
  dplyr::filter(Feature == 'LSS') %>% 
  dplyr::group_by(Feature) %>% 
  dplyr::mutate(Group = ifelse(Value >= median(Value), 1, 0)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(ClinicalInfo_Pat %>%
                     dplyr::transmute(Pat_ID, EFStime, EFS, TLS_Pat = as.integer(TLS_Pat > 0))) %>% 
  dplyr::nest_by(Feature) %>%
  dplyr::mutate(data = data %>% 
                  dplyr::mutate(Group = dplyr::case_when(
                    TLS_Pat == 1 ~ paste0('TLS+'),
                    TLS_Pat == 0 & Group == 0 ~ paste0(Feature, '_low'),
                    TLS_Pat == 0 & Group == 1 ~ paste0(Feature, '_hi')) %>% 
                      forcats::fct_relevel(paste0(Feature, '_low'), 'TLS+', paste0(Feature, '_hi'))) %>% 
                  list()) %>% 
  dplyr::mutate(Result = data %>%
                  survminer::ggsurvplot(
                    data = .,
                    fit = survival::survfit(survival::Surv(EFStime, EFS) ~ Group, data = .),
                    ylab = 'EFS', xlab = 'Months', break.time.by = 24,
                    palette = c('#D6C7C2', '#7DA9B8', '#A94D42'),
                    ggtheme = ggprism::theme_prism(base_size = 12, base_fontface = 'plain'),
                    conf.int = F, pval = F, risk.table = F, tables.height = 0.3,
                  ) %>% 
                  .$plot %>%
                  list()) %>% 
  dplyr::mutate(Result2 = data %>% 
                  dplyr::nest_by() %>% 
                  tidyr::expand_grid(Group_Subset = list(c(1, 2), c(1, 3), c(2, 3))) %>% 
                  dplyr::rowwise() %>% 
                  dplyr::mutate(data_subset = data %>% 
                                  dplyr::filter(as.numeric(Group) %in% Group_Subset) %>% 
                                  survival::survdiff(survival::Surv(EFStime, EFS) ~ Group, data = .) %>% 
                                  broom::glance() %>% 
                                  list()) %>% 
                  dplyr::mutate(Group_Subset = paste0(levels(data$Group)[Group_Subset], collapse = ' vs ')) %>% 
                  tidyr::unnest(data_subset) %>%
                  dplyr::mutate(Group_Subset = 
                                  paste0(Group_Subset, ' = ', sprintf('%.3f', p.value))) %>% 
                  dplyr::pull(Group_Subset) %>% 
                  paste0(collapse = '\n') %>% 
                  list()) %>% 
  dplyr::mutate(Result = {Result +
      annotate(geom = 'text', x = 0, hjust = 0, vjust = 0, y = 0, label = Result2)} %>% 
                  list()) %>% 
  dplyr::pull(Result)
dev.off()
