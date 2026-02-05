setwd("~/projects/NPC_SpatialProteomics/")
library(tidyverse)
library(data.table)
source('0.Helpers.R')

MT_CT <- qs::qread('SingleCellDataFrame.qs')
LEVELS <- qs::qread('AnnotationLevels.qs')
ClinicalInfo_ROI <- data.table::fread('ClinicalInfo.csv')
ClinicalInfo_Pat <- ClinicalInfo_ROI %>% 
  dplyr::select(-ROI_ID, -TLS_ROI) %>% 
  dplyr::distinct()

Interaction_Results <- qs::qread('Interaction_Results.qs')


# Figure S12a
pdf('FigS12a.pdf', height = 4, width = 5)
Interaction_Results %>% 
  dplyr::filter(SA %in% c("Global") & LABEL_USED == 'CT3') %>% 
  dplyr::mutate(PLOT = {INTERACTION %>% 
      dplyr::mutate(rowname = rowname %>% forcats::fct_rev()) %>% 
      dplyr::group_by(rowname, colname) %>% 
      dplyr::summarize(Intensity = mean(Value, na.rm = T)) %>% 
      dplyr::ungroup() %>% 
      ggplot(aes(x = colname, y = rowname)) + 
      geom_tile(aes(fill = Intensity), color = 'black', linewidth = 0.5) +
      scale_x_discrete(expand = expansion()) +
      scale_y_discrete(expand = expansion()) +
      scale_fill_gradient2(name = 'Interaction', limits = c(NA, 1),
                           low = '#1965B0', mid = '#EEEEEE', high = '#CB2314') +
      ggprism::theme_prism(base_size = 8, base_fontface = 'plain') +
      theme(legend.title = element_text(size = 8), 
            legend.key.height = rel(0.5),
            legend.key.width = rel(0.5),
            axis.text.x = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0),
                                       angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
            axis.line = element_blank(), 
            axis.ticks = element_blank()) +
      labs(x = NULL, y = NULL) +
      coord_fixed(ratio = 1)} %>% 
        list()) %>% 
  dplyr::pull(PLOT) %>% 
  purrr::walk(print)
dev.off()


# Figure S12b
pdf('FigS12b.pdf', height = 4, width = 8)
Interaction_Results %>% 
  tidyr::unnest(INTERACTION) %>% 
  dplyr::ungroup() %>% 
  dplyr::transmute(Pat_ID = group_by,
                   SA,
                   rowname, 
                   colname, 
                   Group = Value_Binary %>% as.factor()
  ) %>% 
  dplyr::filter(colname == 'DC' & rowname == 'Epi') %>% 
  dplyr::filter(SA %in% c('Tumor', 'Stroma')) %>% 
  dplyr::mutate(SA = SA %>% forcats::fct_relevel(c('Tumor', 'Stroma'))) %>% 
  dplyr::nest_by(SA) %>% 
  dplyr::mutate(SurvData = data %>%
                  dplyr::left_join(ClinicalInfo_Pat %>% 
                                     dplyr::select(Pat_ID, DMFS, DMFStime), 
                                   by = 'Pat_ID') %>%
                  list()) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(PLOT = {
    survminer::ggsurvplot(
      fit = survival::survfit(survival::Surv(DMFStime, DMFS) ~ Group, 
                              data = SurvData), 
      data = SurvData, xlab = 'Months', ylab = 'DMFS', break.time.by = 48,
      legend.title = 'Factor',
      pval = T, tables.height = 0.3, conf.int = F,
      palette = c('grey70', '#EC7723'),
      risk.table = T)[['plot']] +
      theme(legend.direction = 'vertical', 
            plot.title = element_text(size = 16, hjust = 0.5)) +
      labs(title = paste0('At ', SA)) } %>%
      list()) %>% 
  dplyr::pull(PLOT) %>% 
  patchwork::wrap_plots(nrow = 1, byrow = F, guides = 'collect')
dev.off()


# Figure S12c
pdf('FigS12c.pdf', height = 4, width = 7)
list(Interaction_Results %>%
       dplyr::filter(SA %in% c('Stroma', 'LE', 'Tumor'), LABEL_USED == 'CT3') %>% 
       tidyr::unnest(INTERACTION) %>% 
       dplyr::filter(rowname %in% LEVELS$CT3[9:15] & colname %in% LEVELS$CT3[9:15]) %>% 
       dplyr::filter(rowname == colname) %>% 
       dplyr::group_by(SA, group_by) %>% 
       dplyr::summarize(Value = mean(Value)) %>% 
       dplyr::ungroup() %>% 
       ggplot(aes(x = SA, y = Value)) +
       labs(x = NULL, y = 'Average Intensity',
            title = 'Homotypic lymphocyte-related\ninteractions intensity'),
     
     Interaction_Results %>%
       dplyr::filter(SA %in% c('Stroma', 'LE', 'Tumor'), LABEL_USED == 'CT3') %>% 
       tidyr::unnest(INTERACTION) %>% 
       dplyr::filter(rowname %in% LEVELS$CT3[9:15] & colname %in% LEVELS$CT3[9:15]) %>% 
       dplyr::filter(rowname != colname) %>% 
       dplyr::group_by(SA, group_by) %>% 
       dplyr::summarize(Value = mean(Value)) %>% 
       dplyr::ungroup() %>% 
       ggplot(aes(x = SA, y = Value)) +
       labs(x = NULL, y = 'Average Intensity',
            title = 'Heterotypic lymphocyte-related\ninteractions intensity')
) %>% 
  purrr::map(\(X) X +
               geom_boxplot(aes(fill = SA), outliers = F) +
               geom_jitter(aes(fill = SA), stroke = 0.2, shape = 21, size = 0.5, width = 0.15) +
               ggprism::theme_prism(base_size = 12, base_fontface = 'plain', border = T) +
               ggpubr::stat_compare_means(comparisons = list(c('0', '1')),
                                          label = "p.format",
                                          method = 'wilcox') +
               scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
               scale_fill_manual(values = c("#762A83", "#1AE4B6", "#1965B0") %>% 
                                   setNames(LEVELS$SA[1:3]), guide = 'none') +
               ggpubr::stat_compare_means(method = 'wilcox',
                                          comparisons = list(c('Stroma', 'LE'),
                                                             c('LE', 'Tumor'),
                                                             c('Stroma', 'Tumor'))) +
               theme(legend.title = element_text())
  ) %>% 
  patchwork::wrap_plots(nrow = 1, guides = 'collect')
dev.off()



# Figure S12d
pdf('FigS12d.pdf', height = 7, width = 4)
list(Interaction_Results %>%
       dplyr::filter(SA %in% c('Global'), LABEL_USED == 'CT3') %>% 
       tidyr::unnest(INTERACTION) %>% 
       dplyr::filter(rowname == 'CAF' | colname == 'CAF') %>% 
       dplyr::group_by(rowname, colname) %>% 
       dplyr::summarize(Value = mean(Value), .groups = 'drop') %>% 
       dplyr::arrange(Value) %>% 
       dplyr::mutate(Pair = paste0(rowname, '-', colname) %>% forcats::fct_inorder()) %>% 
       ggplot(aes(y = Pair, x = Value, fill = Value)) + 
       geom_bar(stat = 'identity') + 
       geom_text(aes(label = ifelse(Value > 0.4, sprintf('%.2f', Value), NA)), size = 3, hjust = 1) +
       scale_fill_gradient2(name = 'Interaction', limits = c(-1, 1),
                            low = '#1965B0', mid = '#EEEEEE', high = '#CB2314', guide = 'none') +
       scale_x_continuous(limits = c(-0.5, 1),
                          breaks = c(-0.4, 0, 0.4, 0.8)) + 
       ggprism::theme_prism(base_size = 12, base_fontface = 'plain', border = T) +
       theme(panel.grid.major = element_line(color = '#CCCCCC', linewidth = 0.2, linetype = 'dashed'),
             legend.title = element_text()) +
       labs(x = NULL, y = NULL),
     Interaction_Results %>%
       dplyr::filter(SA %in% c('LE'), LABEL_USED == 'CT4') %>% 
       tidyr::unnest(INTERACTION) %>% 
       dplyr::filter(rowname %in%  LEVELS[['CT4']][9:11] | colname %in%  LEVELS[['CT4']][9:11]) %>% 
       dplyr::group_by(rowname, colname) %>% 
       dplyr::summarize(Value = mean(Value), .groups = 'drop') %>% 
       dplyr::arrange(Value) %>% 
       dplyr::slice_tail(n = 5) %>% 
       dplyr::mutate(Pair = paste0(rowname, '-', colname) %>% forcats::fct_inorder()) %>% 
       ggplot(aes(y = Pair, x = Value, fill = Value)) + 
       geom_bar(stat = 'identity') + 
       geom_text(aes(label = ifelse(Value > 0.4, sprintf('%.2f', Value), NA)), size = 3, hjust = 1) +
       scale_fill_gradient2(name = 'Interaction', limits = c(-1, 1),
                            low = '#1965B0', mid = '#EEEEEE', high = '#CB2314', guide = 'none') +
       scale_x_continuous(limits = c(-0.5, 1),
                          breaks = c(-0.4, 0, 0.4, 0.8)) + 
       ggprism::theme_prism(base_size = 12, base_fontface = 'plain', border = T) +
       theme(panel.grid.major = element_line(color = '#CCCCCC', linewidth = 0.2, linetype = 'dashed'),
             legend.title = element_text()) +
       labs(x = NULL, y = NULL)
) %>% 
  patchwork::wrap_plots(ncol = 1, heights = c(29, 5))
dev.off()


# Figure S12e
pdf('FigS12e.pdf', height = 6, width = 8)
Interaction_Results %>% 
  dplyr::filter(SA %in% c('Stroma', 'LE', 'Tumor'), LABEL_USED == 'CT3') %>% 
  dplyr::select(SA, INTERACTION) %>% 
  tidyr::unnest(INTERACTION) %>% 
  dplyr::mutate(SA = SA %>% forcats::fct_relevel('LE')) %>% 
  dplyr::mutate(rowname = rowname %>% forcats::fct_drop() %>% forcats::fct_rev()) %>% 
  dplyr::mutate(colname = colname %>% forcats::fct_drop()) %>% 
  dplyr::group_by(rowname, colname, SA) %>% 
  dplyr::summarize(Intensity = mean(Value, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(rowname, colname) %>% 
  dplyr::slice_max(order_by = Intensity, n = 1, with_ties = F) %>% 
  dplyr::ungroup() %T>% {temp_count <<- dplyr::count(., SA); temp_label <<- temp_count %>% 
    dplyr::mutate(SA2 = paste0(SA, ': ', n, ' (', round(100 * n / sum(n)), '%)')) %>% 
    {stats::setNames(.$SA2, .$SA)}} %>% 
  dplyr::mutate(SA = SA %>% dplyr::recode(!!!temp_label))%>%
  dplyr::mutate(xmin = as.integer(colname) - 0.5,
                xmax = as.integer(colname) - 0.5 + 1,
                ymin = as.integer(rowname) - 0.5,
                ymax = as.integer(rowname) - 0.5 + 1) %>% 
  ggplot(aes(x = colname, y = rowname)) +
  geom_point(aes(size = Intensity, color = SA)) +
  ggprism::theme_prism(base_size = 12, base_fontface = 'plain', border = T) +
  scale_color_manual(name = 'Spatial annotation with\nstrongest interaction', 
                     values = c("#1AE4B6", "#762A83", "#1965B0")) +
  scale_size_continuous(range = c(0, 10), limits = c(0, 1)) +
  theme(legend.title = element_text(), 
        axis.text.x = element_text(#margin = margin(t = 0, r = 0, b = 0, l = 0),
          angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_blank(), 
        axis.ticks = element_blank()) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 1, clip = 'off')
dev.off()


# Figure S12f
pdf('FigS12f.pdf', height = 4, width = 8)
Interaction_Results %>% 
  dplyr::filter(LABEL_USED == 'CT3') %>% 
  dplyr::select(SA, INTERACTION) %>% 
  tidyr::unnest(INTERACTION) %>% 
  base::split(.$SA == 'Global') %>% 
  with(dplyr::left_join(.[[1]], 
                        .[[2]] %>% dplyr::transmute(rowname, colname, group_by, Value_Global = Value),
                        by = c('rowname', 'colname', 'group_by'))) %>% 
  dplyr::group_by(rowname, colname, SA) %>% 
  dplyr::summarize(Accuracy = mean(Value_Global == Value), 
                   Intensity = mean(Value)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(Cancer = rowname %in% LEVELS$CT3[[1]] | colname %in% LEVELS$CT3[[1]]) %>% 
  dplyr::mutate(Stromal = rowname %in% LEVELS$CT3[2:3] | colname %in% LEVELS$CT3[2:3]) %>% 
  dplyr::mutate(Myeloid = rowname %in% LEVELS$CT3[4:8] | colname %in% LEVELS$CT3[4:8]) %>% 
  dplyr::mutate(Lymphoid = rowname %in% LEVELS$CT3[9:15] | colname %in% LEVELS$CT3[9:15]) %>% 
  tidyr::pivot_longer(cols = c(Cancer, Stromal, Myeloid, Lymphoid), names_to = 'Lineage') %>% 
  dplyr::mutate(Lineage = Lineage %>% paste0('-related') %>% forcats::fct_inorder()) %>% 
  dplyr::filter(value) %>% 
  ggplot(aes(x = SA, y = Accuracy)) +
  facet_wrap(~ Lineage, nrow = 1)+
  geom_boxplot(aes(fill = SA)) + 
  ggprism::theme_prism(base_fontface = 'plain', border = T) + 
  scale_fill_manual(values = c("#762A83", "#1AE4B6", "#1965B0") %>% 
                      setNames(LEVELS$SA[1:3]), na.value = 'grey70', guide = 'none') +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)),
                     breaks = seq(0, 1, by = 0.25)) + 
  theme(axis.text.x = element_text(vjust = 1, hjust = 0.5),
        legend.position = 'none',
        strip.text = element_text(size = 12)) + 
  ggpubr::stat_compare_means(comparisons = list(c('Stroma', 'LE'), c('LE', 'Tumor')),
                             tip.length = 0.01, step.increase = 0, size = 3, 
                             method = 'wilcox', paired = T) +
  labs(x = NULL, y = NULL, title = 'Accuracy of predicting global interaction intensities')
dev.off()

