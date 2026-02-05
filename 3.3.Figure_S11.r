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



# Figure S11a
Net_Result_Nodes <- Net_Result %>%
  dplyr::filter((GRAPH_METRICS %>% tidygraph::activate(edges) %>% dplyr::as_tibble() %>% nrow()) > 0) %>% 
  dplyr::mutate(DATA = GRAPH_METRICS %>% 
                  tidygraph::activate(nodes) %>%
                  dplyr::as_tibble() %>% 
                  dplyr::filter(n_group_size >= 10) %>%
                  list()) %>% 
  dplyr::select(ROI_ID, DATA) %>% 
  tidyr::unnest()

pdf('FigS11a.pdf', height = 6, width = 8)
Net_Result_Nodes %>% 
  dplyr::transmute(ROI_ID, 
                   Net_ID = paste0(ROI_ID, '_', sprintf('%02.f', n_group_id)), 
                   n_group_size) %>% 
  dplyr::distinct() %>% 
  dplyr::count(ROI_ID, name = 'Net Count') %>% 
  dplyr::mutate(ROI_ID = ROI_ID %>% forcats::fct_expand(TLSProfile_ROI$ROI_ID)) %>% 
  tidyr::complete(ROI_ID, fill = list(`Net Count` = 0)) %>% 
  dplyr::left_join(TLSProfile_ROI %>%
                     dplyr::transmute(ROI_ID, 
                                      TLS_group)) %>% 
  dplyr::count(TLS_group, `Net Count`, name = 'ROI Count') %>% 
  ggplot() + 
  geom_bar(aes(x = `Net Count`, y = `ROI Count`, fill = TLS_group), 
           color = 'black', width = 1, position = 'stack', stat = 'identity') +
  geom_vline(xintercept = 41, alpha = 0.5, linetype = 'dotted') +
  annotate(geom = 'text', x = 41, y = Inf, size = 6,
           hjust = 1.1, vjust = 2, label = 'Max Count: 41') +
  ggprism::theme_prism(base_size = 16, base_fontface = 'plain', border = T) + 
  theme(legend.position = c(0.3, 0.8)) +
  scale_x_continuous(expand = expansion(c(0.02, 0.02))) + 
  scale_y_continuous(expand = expansion(c(0, 0.02))) + 
  scale_fill_manual(values = c('#D6C7C2', '#7DA9B8', '#A94D42')) + 
  coord_fixed(ratio = 0.6)
dev.off()


# Figure S11b
NetProfile_ROI <- fread('NetProfile_ROI.csv')

pdf('FigS11b.pdf', height = 4, width = 4)
NetProfile_ROI %>% 
  dplyr::left_join(TLSProfile_ROI %>% dplyr::transmute(ROI_ID, TLS_group)) %>%
  ggplot(aes(x = TLS_group, y = LymAbund)) +
  geom_violin(aes(fill = TLS_group)) +
  geom_boxplot(width = 0.1, outliers = F) +
  ggpubr::stat_compare_means(method = 'wilcox',
                             comparisons = list(
                               c('TLS_negative', 'TLS_adjacent'),
                               c('TLS_adjacent', 'TLS_positive'),
                               c('TLS_negative', 'TLS_positive'))
  ) +
  scale_fill_manual(values = c('#D6C7C2', '#7DA9B8', '#A94D42'), guide = 'none') + 
  scale_y_continuous(expand = expansion(c(0.05, 0.1))) +
  ggprism::theme_prism(base_size = 12, base_fontface = 'plain', border = T) +
  labs(x = NULL, y = 'LymAbund')
dev.off()


# Figure S11c
Net_Assort <- Net_Result %>%
  dplyr::filter((GRAPH_METRICS %>% tidygraph::activate(edges) %>% dplyr::as_tibble() %>% nrow()) > 0) %>% 
  dplyr::mutate(GRAPH_SPLIT = GRAPH_METRICS %>% 
                  tidygraph::activate(nodes) %>%
                  dplyr::filter(n_group_size >= 10) %>% 
                  tidygraph::to_split(n_group_id, split_by = 'nodes') %>% 
                  list()) %>% 
  dplyr::mutate(GRAPH_MEASUREMENTS = GRAPH_SPLIT %>%
                  purrr::map(.progress = T, \(X) X %>% 
                               dplyr::mutate(graph_assortativity = tidygraph::graph_assortativity(CT3, directed = F)) %>% 
                               dplyr::as_tibble() %>% 
                               dplyr::select(tidyr::starts_with('graph_')) %>% 
                               dplyr::distinct()
                  ) %>% 
                  tibble::enframe(name = 'n_group_id', value = 'Metrics') %>% 
                  tidyr::unnest() %>% 
                  list())

pdf('FigS11c.pdf', height = 4, width = 4)
Net_Assort %>% 
  dplyr::select(ROI_ID, GRAPH_MEASUREMENTS) %>% 
  tidyr::unnest() %>% 
  dplyr::left_join(TLSProfile_ROI %>%
                     dplyr::transmute(ROI_ID, TLS_group)) %>%
  ggplot(aes(x = TLS_group, y = graph_assortativity)) +
  geom_violin(aes(fill = TLS_group)) +
  geom_boxplot(width = 0.1, outliers = F) +
  ggpubr::stat_compare_means(method = 'wilcox',
                             comparisons = list(
                               c('TLS_negative', 'TLS_adjacent'),
                               c('TLS_adjacent', 'TLS_positive'),
                               c('TLS_negative', 'TLS_positive'))
  ) +
  scale_fill_manual(values = c('#D6C7C2', '#7DA9B8', '#A94D42'), guide = 'none') + 
  scale_y_continuous(expand = expansion(c(0.05, 0.1))) +
  ggprism::theme_prism(base_size = 12, base_fontface = 'plain', border = T) +
  labs(x = NULL, y = 'graph_assortativity')
dev.off()


# Figure S11d
pdf('FigS11d.pdf', height = 4, width = 4)
NetProfile_ROI %>% 
  dplyr::left_join(TLSProfile_ROI %>% dplyr::transmute(ROI_ID, TLS_group)) %>%
  ggplot(aes(x = TLS_group, y = LRI)) +
  geom_violin(aes(fill = TLS_group)) +
  geom_boxplot(width = 0.1, outliers = F) +
  ggpubr::stat_compare_means(method = 'wilcox',
                             comparisons = list(
                               c('TLS_negative', 'TLS_adjacent'),
                               c('TLS_adjacent', 'TLS_positive'),
                               c('TLS_negative', 'TLS_positive'))
  ) +
  scale_fill_manual(values = c('#D6C7C2', '#7DA9B8', '#A94D42'), guide = 'none') + 
  scale_y_continuous(expand = expansion(c(0.05, 0.1))) +
  ggprism::theme_prism(base_size = 12, base_fontface = 'plain', border = T) +
  labs(x = NULL, y = 'LRI')
dev.off()


# Figure S11e
NetProfile_Pat <- fread('NetProfile_Pat.csv')
NetProfile_Surv <- NetProfile_Pat %>%
  tidyr::pivot_longer(cols = !Pat_ID, names_to = 'Feature', values_to = 'Value') %>% 
  dplyr::left_join(ClinicalInfo_Pat) %>% 
  dplyr::group_by(Feature) %>% 
  dplyr::mutate(Group = ifelse(Value >= median(Value), 1, 0)) %>% 
  dplyr::ungroup() %>% 
  dplyr::nest_by(Feature) %>%
  dplyr::mutate(data = data %>% 
                  dplyr::mutate(Group = dplyr::case_when(
                    TLS_Pat > 0 ~ paste0('TLS+'),
                    TLS_Pat == 0 & Group == 0 ~ paste0(Feature, '_low'),
                    TLS_Pat == 0 & Group == 1 ~ paste0(Feature, '_hi')) %>% 
                      forcats::fct_relevel(paste0(Feature, '_low'), 'TLS+', paste0(Feature, '_hi'))) %>% 
                  list()) %>% 
  dplyr::mutate(plot = data %>%
                  survminer::ggsurvplot(
                    data = .,
                    fit = survival::survfit(survival::Surv(EFStime, EFS) ~ Group, data = .),
                    ylab = 'Event-free Survival', xlab = 'Months', break.time.by = 24, 
                    # legend.title = Subset,
                    # legend.labs = paste0(Feature, c(' Low', ' High')),
                    palette = c('#D6C7C2', '#7DA9B8', '#A94D42'),
                    ggtheme = ggprism::theme_prism(base_size = 12, base_fontface = 'plain'),
                    conf.int = F, pval = F,
                  ) %>% 
                  .[['plot']] %>%
                  list()) %>% 
  dplyr::mutate(test_result = data %>% 
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
                                  paste0(Group_Subset, ' = ', sprintf('%.2g', p.value))) %>% 
                  dplyr::pull(Group_Subset) %>% 
                  paste0(collapse = '\n') %>% 
                  list()) %>% 
  dplyr::mutate(plot = {plot  +
      annotate(geom = 'text', x = 0, hjust = 0, vjust = 0, y = 0, label = test_result)} %>% 
        list())

pdf('FigS11e.pdf', height = 4, width = 5)
NetProfile_Surv %>% 
  dplyr::filter(Feature == 'LRI') %>% 
  .[['plot']] %>% 
  print()
dev.off()


# Figure S11f
pdf('FigS11f.pdf', height = 4, width = 5)
NetProfile_Surv %>% 
  dplyr::filter(Feature == 'NetCountPerROI') %>% 
  .[['plot']] %>% 
  print()
dev.off()


# Figure S11g
pdf('FigS11g.pdf', height = 4, width = 5)
NetProfile_Surv %>% 
  dplyr::filter(Feature == 'MaxNetSize') %>% 
  .[['plot']] %>% 
  print()
dev.off()


# Figure S11h
pdf('FigS11h.pdf', height = 4, width = 5)
NetProfile_Surv %>% 
  dplyr::filter(Feature == 'MeanNetSize') %>% 
  .[['plot']] %>% 
  print()
dev.off()


# Figure S11i
pdf('FigS11i.pdf', height = 4, width = 4)
NetProfile_Pat %>%
  tidyr::pivot_longer(cols = !Pat_ID, names_to = 'Feature', values_to = 'Value') %>% 
  dplyr::group_by(Feature) %>% 
  dplyr::mutate(Group = ifelse(Value >= median(Value), 1, 0)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(ClinicalInfo_Pat %>% dplyr::mutate(TLS_Presence = as.integer(TLS_Pat > 0))) %>% 
  dplyr::nest_by(Feature) %>%
  dplyr::mutate(Result = data %>% 
                  survival::coxph(survival::Surv(EFStime, EFS) ~ Group + TLS_Presence, data = .) %>% 
                  summary() %>% 
                  list()) %>% 
  dplyr::mutate(Result = dplyr::left_join(Result$conf.int %>% as.data.frame() %>% tibble::rownames_to_column('Var'),
                                          Result$coefficients %>% as.data.frame() %>% tibble::rownames_to_column('Var'),
                                          by = c('Var', 'exp(coef)')) %>% 
                  list()) %>% 
  dplyr::select(Feature, Result) %>%
  tidyr::unnest() %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(Var = ifelse(Var == 'TLS_Presence', Var, Feature) %>%
                  forcats::fct_relevel('TLS_Presence')) %>%
  dplyr::mutate(Feature = Feature %>% forcats::fct_relevel('LymAbund', 'LRI', 'LSS', 'NetCountPerROI', 'MaxNetSize', 'MeanNetSize')) %>%
  ggplot(aes(x = `exp(coef)`, y = Var)) +
  facet_wrap(.~ Feature, ncol = 1, scales = 'free_y') +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'black') +
  geom_pointrange(aes(xmin = `lower .95`, xmax = `upper .95`)) +
  geom_text(aes(x = Inf, label = paste0('  ', sprintf('%.2g', `Pr(>|z|)`))), hjust = 0) +
  theme_bw() + 
  theme(strip.text = element_blank(), 
        axis.text.y = element_text(color = 'black'),
        plot.margin = unit(c(1,3,1,1), 'line')) +
  scale_x_continuous(transform = 'log2') +
  coord_cartesian(clip = 'off') +
  labs(x = 'Hazard Ratio for EFS', y = NULL)
dev.off()

