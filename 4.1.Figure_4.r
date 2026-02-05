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

Interaction_Results <- qs::qread('Interaction_Results.qs')


# Figure 4a
pdf('Fig4a.pdf', height = 4, width = 5)
Interaction_Results %>% 
  dplyr::filter(SA %in% c("Stroma", "LE", "Tumor") & LABEL_USED == 'CT3') %>% 
  tidyr::unnest(INTERACTION) %>% 
  dplyr::group_by(SA, rowname, colname) %>% 
  dplyr::summarize(Intensity = mean(Value, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(rev(rowname), SA) %>% 
  dplyr::mutate(rowname2 = paste0(rowname, '_', SA) %>% forcats::fct_inorder()) %>% 
  {ggplot(data = .) + 
      geom_tile(aes(x = colname, y = rowname2, fill = Intensity), color = NA, linewidth = 0.5) +
      geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                fill = NA, color = 'black',
                data = dplyr::distinct(., rowname, colname) %>% 
                  dplyr::mutate(ymin = as.integer(rowname) * 3 - 2.5,
                                ymax = as.integer(rowname) * 3 + 3 - 2.5,
                                xmin = as.integer(colname) - 0.5,
                                xmax = as.integer(colname) + 1 - 0.5)) +
      scale_x_discrete(expand = expansion()) +
      scale_y_discrete(expand = expansion(), labels = as.vector(rbind("", rev(LEVELS$CT3)[-1], ""))) +
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
      coord_fixed(ratio = 1/3)} %>% 
  print()
dev.off()


# Figure 4b
pdf('Fig4b.pdf', height = 4, width = 8)
Interaction_Results %>% 
  dplyr::filter(SA %in% c("Global") & LABEL_USED == 'CT3') %>% 
  tidyr::unnest(INTERACTION) %>% 
  dplyr::mutate(rowname = rowname %>% forcats::fct_rev()) %>%
  dplyr::mutate(colname = colname %>% forcats::fct_rev()) %>%
  dplyr::group_by(rowname, colname) %>% 
  dplyr::summarize(Intensity = mean(Value, na.rm = T)) %>% 
  dplyr::ungroup() %>% 
  {dplyr::bind_rows(dplyr::filter(., colname == 'Epi' & rowname != colname) %>% 
                      dplyr::transmute(TME = rowname, 
                                       TME_as = 'rowname',
                                       Intensity),
                    dplyr::filter(., rowname == 'Epi' & rowname != colname) %>% 
                      dplyr::transmute(TME = colname, 
                                       TME_as = 'colname',
                                       Intensity))} %>% 
  dplyr::mutate(TME_as = TME_as %>% forcats::fct_relevel(c('rowname', 'colname'))) %>% 
  dplyr::arrange(Intensity) %>% 
  ggplot(aes(x = TME_as, y = TME))+
  geom_segment(aes(color = Intensity, linewidth = abs(Intensity), xend = 1.5, yend = 7.5)) +
  geom_point(aes(fill = TME), shape = 21, size = 6) +
  geom_text(aes(x = as.integer(TME_as) + 0.05 * (as.integer(TME_as) * 2 - 3), 
                label = paste0(TME, ': ', sprintf('%.2f', Intensity)), hjust = 2 - 1 * as.integer(TME_as))) + 
  annotate(geom = 'point', x = 1.5, y = 7.5, fill = '#3B597AFF', shape = 21, size = 10) +
  annotate(geom = 'text', x = 1.5, y = 9, label = 'Cancer') +
  scale_linewidth(limits = c(0, 1), range = c(1, 5)) + 
  scale_color_gradient2(name = 'Interaction', limits = c(-1, 1),
                        low = '#1965B0', mid = '#EEEEEE', high = '#CB2314') +
  scale_fill_manual(values = c('#CCE5E2', '#AFD1CE',
                               '#ECE6CE', '#DBD4B2', '#CBC196', '#BBAF7A', '#AA9C5D',
                               '#EEDAD9', '#DEC2C1', '#CEAAA8', '#BF9290', '#AF7978', '#9F6160', '#8F4947') %>%
                      rev(), guide = 'none') +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 0.075)
dev.off()


# Figure 4c
ExpProfile_MyeloidCell <- MT_CT %>% 
  dplyr::filter(CT2 == 'Myeloid') %>% 
  dplyr::rename(CC3 = Caspase3cleaved) %>% 
  dplyr::group_by(., Pat_ID, CT3) %>% 
  dplyr::summarize(dplyr::across(c(Ki67, HLADR, PDL1, CC3), .fns = mean)) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_longer(cols = c(Ki67, HLADR, PDL1, CC3), names_to = 'Channel', values_to = 'MeanExp') %>% 
  dplyr::mutate(Channel = Channel %>% forcats::fct_relevel(c('Ki67', 'HLADR', 'PDL1', 'CC3'))) %>% 
  dplyr::rename(TME = CT3) %>% 
  dplyr::nest_by(TME, Channel, .key = 'ExpData') %>% 
  dplyr::ungroup()

pdf('Fig4c.pdf', height = 4, width = 6)
Interaction_Results %>% 
  dplyr::filter(SA %in% c("Global") & LABEL_USED == 'CT3') %>% 
  tidyr::unnest(INTERACTION) %>% 
  dplyr::transmute(Pat_ID = group_by,
                   rowname, 
                   colname, 
                   Value_Binary
  ) %>% 
  dplyr::filter(colname == 'Epi' & rowname != colname) %>% 
  dplyr::rename(TME = rowname) %>% 
  dplyr::mutate(TME_as = 'rowname') %>% 
  dplyr::select(-colname) %>% 
  dplyr::mutate(TME_as = TME_as %>% forcats::fct_relevel(c('colname', 'rowname'))) %>% 
  dplyr::nest_by(TME_as, TME) %>% 
  dplyr::ungroup() %>% 
  dplyr::right_join(ExpProfile_MyeloidCell, by = c('TME')) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(data = data %>% 
                  dplyr::left_join(ExpData, by = c('Pat_ID')) %>% 
                  list()) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(result = 
                  purrr::map(.$data, .progress = T, \(X) {
                    purrr::safely(dplyr::left_join)(rstatix::wilcox_test(X, MeanExp ~ Value_Binary),
                                                    rstatix::cohens_d(X, MeanExp ~ Value_Binary),
                                                    by = join_by(.y., group1, group2, n1, n2))$result
                  })) %>% 
  tidyr::unnest(result) %>% 
  dplyr::mutate(Pair = ifelse(TME_as == 'rowname', 
                              paste0(TME, '-Epi'), 
                              paste0('Epi-', TME)) %>% 
                  forcats::fct_inorder() %>%
                  forcats::fct_rev()) %>% 
  dplyr::mutate(log10(p)) %>% 
  ggplot(aes(y = Pair, x = Channel)) +
  geom_point(aes(alpha = ifelse(p < 0.05, 'Sig', 'Non_Sig'),
                 fill = - sign(effsize) * log10(p), 
                 size = -log10(p)), shape = 21) + 
  scale_fill_gradient(name = 'log10(p)', high = '#EEEEEE', low = '#CB2314', na.value = '#EEEEEE') + 
  scale_alpha_manual(name = 'Sig', values = c(0, 1), guide = 'none') + 
  scale_size_continuous(range = c(1, 10)) + 
  ggprism::theme_prism(base_size = 12, border = T, base_fontface = 'plain') +
  theme(legend.title = element_text(), 
        panel.grid = element_line(colour = "grey92", linewidth = rel(0.5))) +
  labs(x = NULL, y = NULL, title = 'Expression on myeloid cells')
dev.off()


# Figure 4d
ExpProfile_DC_SA <- MT_CT %>% 
  dplyr::filter(CT3 == 'DC') %>% 
  dplyr::group_by(Pat_ID, SA) %>% 
  dplyr::summarize(HLADR = mean(HLADR)) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(SA != 'Noise') %>% 
  dplyr::nest_by(SA, .key = 'ExpData') %>% 
  dplyr::ungroup()

pdf('Fig4d.pdf', height = 4, width = 6)
Interaction_Results %>% 
  dplyr::filter(SA %in% c("Stroma", "LE", "Tumor") & LABEL_USED == 'CT3') %>% 
  dplyr::mutate(SA = SA %>% forcats::fct_relevel(c("LE", "Stroma", "Tumor"))) %>% 
  tidyr::unnest(INTERACTION) %>% 
  dplyr::ungroup() %>% 
  dplyr::transmute(Pat_ID = group_by,
                   SA,
                   rowname, 
                   colname, 
                   Group = Value_Binary %>% as.factor()
  ) %>% 
  dplyr::filter(rowname == 'Epi' & colname == 'DC') %>% 
  dplyr::nest_by(SA) %>% 
  dplyr::left_join(ExpProfile_DC_SA, by = 'SA') %>% 
  dplyr::mutate(ExpData = data %>%
                  dplyr::left_join(ExpData, by = 'Pat_ID') %>%
                  list()) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(PLOT = {ExpData %>% 
      ggplot(aes(x = Group, y = HLADR)) + 
      geom_boxplot(aes(fill = Group), outliers = F) +
      geom_jitter(aes(fill = Group), stroke = 0.2, shape = 21, size = 0.5, width = 0.15) +
      ggpubr::stat_compare_means(comparisons = list(c('0', '1')),
                                 label = "p.format",
                                 method = 'wilcox') +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      scale_fill_manual(values = c('grey70', '#EC7723') %>% setNames(c('0', '1'))) +
      ggprism::theme_prism(base_size = 10, base_fontface = 'plain', border = T) +
      labs(x = paste0('At ', SA),
           y = paste0('HLA-DR expression of DC'))
      } %>%
        list()) %>% 
  dplyr::pull(PLOT) %>% 
  patchwork::wrap_plots(nrow = 1, guides = 'collect') &
  theme(legend.title = element_blank(), 
        legend.position = 'top',
        plot.title = element_text(size = 10))
dev.off()


# Figure 4e
pdf('Fig4e.pdf', height = 4, width = 8)
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
  dplyr::filter(SA %in% c('LE', 'Global')) %>% 
  dplyr::mutate(SA = SA %>% forcats::fct_relevel(c('LE', 'Global'))) %>% 
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


# Figure 4f
pdf('Fig4f.pdf', height = 4, width = 7)
list(MT_CT %>%
       Helper2(UNIT_ID = 'Pat_ID', TERM_A = 'CT3', TERM_B = 'SA',
               REL2A = F, REL2B = T) %>% 
       dplyr::filter(Term_A == 'CAF') %>% 
       dplyr::filter(Term_B %in% c('Stroma', 'LE', 'Tumor')) %>%
       dplyr::mutate(SA = Term_B,
                     Value = ifelse(Prop4Filter < 0.01, 0, Value)) %>% 
       ggplot(aes(x = SA, y = Value)) + 
       labs(x = NULL, y = 'CAF Proportion'),
     
     Interaction_Results %>%
       dplyr::filter(SA %in% c('Stroma', 'LE', 'Tumor'), LABEL_USED == 'CT3') %>% 
       tidyr::unnest(INTERACTION) %>% 
       dplyr::filter(rowname == 'CAF' | colname == 'CAF') %>% 
       dplyr::group_by(SA, group_by) %>% 
       dplyr::summarize(Value = mean(Value)) %>% 
       dplyr::ungroup() %>% 
       ggplot(aes(x = SA, y = Value)) +
       labs(x = NULL, y = 'CAF-related\ninteractions intensity')) %>% 
  purrr::map(\(X) X +
               geom_boxplot(aes(fill = SA), outliers = F) +
               geom_jitter(aes(fill = SA), stroke = 0.2, shape = 21, size = 0.5, width = 0.15) +
               ggprism::theme_prism(base_size = 12, base_fontface = 'plain', border = T) +
               ggpubr::stat_compare_means(comparisons = list(c('0', '1')),
                                          label = "p.format",
                                          method = 'wilcox') +
               scale_y_continuous(expand = expansion(mult = c(0.01, 0.08))) +
               scale_fill_manual(values = c("#762A83", "#1AE4B6", "#1965B0") %>% 
                                   setNames(LEVELS$SA[1:3]), guide = 'none') +
               ggpubr::stat_compare_means(method = 'wilcox',
                                          comparisons = list(c('Stroma', 'LE'),
                                                             c('LE', 'Tumor'),
                                                             c('Stroma', 'Tumor'))) +
               theme(legend.title = element_text(), plot.title = element_text(size = 10))
  ) %>% 
  patchwork::wrap_plots(guides = 'collect')
dev.off()


# Figure 4g
AllFeatures <- qs::qread('AllFeatures.qs')
TherapeuticBenefits_Results <- AllFeatures %>%
  dplyr::group_by(SA_group, Feature, Label) %>% 
  dplyr::mutate(Value = as.integer(Value > median(Value))) %>% 
  dplyr::group_by(SA_group, Feature, Label) %>%
  dplyr::mutate(ndis = dplyr::n_distinct(Value)) %>%
  dplyr::filter(ndis != 1) %>%
  dplyr::ungroup() %>% 
  dplyr::left_join(ClinicalInfo_Pat, by = 'Pat_ID') %>%
  dplyr::mutate(Group = Value) %>% 
  dplyr::nest_by(SA_group, Feature, Label) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(p_coxph_Interaction = survival::coxph(data = data, formula = survival::Surv(EFStime, EFS) ~ Group + InductChemo + Group * InductChemo) %>% 
                  summary() %>%
                  .$coef %>% 
                  .['Group:InductChemo', 'Pr(>|z|)']) %>% 
  dplyr::mutate(Details = data %>%
                  split(.$Group) %>% 
                  purrr::imap_dfr(\(X, Y) {
                    result_coxph <- survival::coxph(data = X, formula = survival::Surv(EFStime, EFS) ~ InductChemo)
                    result_survdiff <-  survival::survdiff(data = X, formula = survival::Surv(EFStime, EFS) ~ InductChemo)
                    return(list(
                      Group = Y,
                      p_logrank = result_survdiff %>% broom::glance() %>% dplyr::pull(p.value),
                      p_coxph = result_coxph %>% summary() %>% .$coef %>% .[, 'Pr(>|z|)'],
                      HR = result_coxph %>% summary() %>% .$coef %>% .[, 'exp(coef)'],
                      HR_CI_l = result_coxph %>% summary() %>% .$conf.int %>% .[, 'lower .95'],
                      HR_CI_u = result_coxph %>% summary() %>% .$conf.int %>% .[, 'upper .95']
                    ))
                  }) %>%
                  list()) %>% 
  dplyr::ungroup()

pdf('Fig4g.pdf', width = 6, height = 8)
TherapeuticBenefits_Results %>% 
  dplyr::filter(p_coxph_Interaction < 0.05) %>% 
  dplyr::rowwise() %>%
  dplyr::filter((Details %>% dplyr::filter(p_coxph < 0.05) %>% dplyr::pull(Group)) == 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(p_coxph_Interaction)) %>% 
  tidyr::unnest(Details, names_sep = '_') %>% 
  dplyr::mutate(Details_Group = Details_Group %>% as.factor()) %>% 
  dplyr::mutate(Label = Label %>% forcats::fct_inorder()) %>% 
  dplyr::mutate(Alpha = ifelse(Details_p_logrank < 0.05, 1, 0.2)) %>% 
  dplyr::mutate(SA_group = SA_group %>% dplyr::recode(`Spatially_Ignorant` = 'Spatially\nIgnorant\n')) %>% 
  {ggplot(., aes(y = Label)) +
      geom_vline(xintercept = 1, alpha = 0.5, linetype = 'dotted') +
      geom_pointrange(aes(
        x = Details_HR, 
        xmin = Details_HR_CI_l, 
        xmax = Details_HR_CI_u, 
        color = Details_Group, 
        alpha = Alpha), 
        size = 0.5,
        position = position_dodge(width = 1)) +
      geom_text(aes(x = 4, 
                    label = sprintf('%.3f', p_coxph_Interaction)), 
                size = 3, 
                position = position_dodge(width = 1)) +
      facet_grid(SA_group ~ ., scales = 'free_y', space = 'free', switch = 'y') + 
      scale_alpha_identity(guide = "none") +
      scale_color_manual(values = c('#b58463', '#6A3D9A')) +
      scale_x_continuous(trans = 'log10', 
                         labels = scales::label_number(drop0trailing = T),
                         expand = expansion(mult = c(0.05, 0.15))) +
      scale_y_discrete(expand = expansion(mult = c(0, 0)), 
                       position = 'right', labels = scales::parse_format()) +
      ggprism::theme_prism(border = T, base_size = 10, base_fontface = 'plain') +
      theme(axis.ticks.y = element_line(linewidth = 0.2),
            strip.text.y.left = element_text(size = 10)) + 
      labs(y = NULL, x = 'Hazard Ratios')
  }
dev.off()

