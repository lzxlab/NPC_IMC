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


# Figure S13a
AllFeatures <- qs::qread('AllFeatures.qs')

Prognosis_Results <- AllFeatures %>%
  dplyr::group_by(SA_group, Feature, Label) %>% 
  dplyr::mutate(Group = as.integer(Value > median(Value))) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(ClinicalInfo_Pat, by = 'Pat_ID') %>%
  dplyr::nest_by(SA_group, Feature, Label) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(coxph = survival::coxph(data = data, formula = survival::Surv(EFStime, EFS) ~ Group) %>% list(),
                survdiff = survival::survdiff(data = data, formula = survival::Surv(EFStime, EFS) ~ Group) %>% list()) %>% 
  dplyr::mutate(p_logrank = survdiff %>% broom::glance() %>% dplyr::pull(p.value),
                p_coxph = coxph %>% summary() %>% .$coef %>% .['Group', 'Pr(>|z|)'],
                HR = coxph %>% summary() %>% .$coef %>% .[, 'exp(coef)'],
                HR_CI_l = coxph %>% summary() %>% .$conf.int %>% .[, 'lower .95'],
                HR_CI_u = coxph %>% summary() %>% .$conf.int %>% .[, 'upper .95']) %>% 
  dplyr::ungroup()


pdf('FigS13a.pdf', height = 8, width = 12)
Prognosis_Results %>% 
  dplyr::mutate(Label = ifelse(p_coxph < 0.05, Label, NA)) %>%
  dplyr::arrange(desc(p_coxph)) %>% 
  dplyr::mutate(nudge = dplyr::case_when(
    HR > 1 & SA_group == 'Stroma' ~ 1,
    HR > 1 & SA_group == 'LE' ~ 1,
    HR > 1 & SA_group == 'Tumor' ~ 0.5,
    HR > 1 & SA_group == 'Lymphoid_Net' ~ 0,
    HR > 1 & SA_group == 'Cell_Abundance' ~ 1,
    HR < 1 & SA_group == 'Stroma' ~ -0.5,
    HR < 1 & SA_group == 'LE' ~ -1,
    HR < 1 & SA_group == 'Tumor' ~ -1,
    HR < 1 & SA_group == 'Lymphoid_Net' ~ -0.5,
    HR < 1 & SA_group == 'Cell_Abundance' ~ -1
  )) %>% 
  ggplot(aes(y = -log10(p_coxph), x = log10(HR))) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = log10(HR_CI_l), xmax = log10(HR_CI_u), color = SA_group), 
                  alpha = 0.25, 
                  size = 0.3) +
  geom_point(aes(x = log10(HR), color = SA_group, alpha = (p_coxph < 0.05)), size = 2) +
  ggrepel::geom_label_repel(aes(label = Label, 
                                fill = SA_group, 
                                hjust = as.integer(HR < 1), 
                                nudge_x = nudge), 
                            direction = 'y', force = 1, 
                            box.padding = 0.1,
                            seed = 12, size = 2, 
                            min.segment.length = 0, max.overlaps = Inf, parse = T) +
  scale_alpha_manual(values = c(0.1, 1), guide = "none") +
  scale_color_manual(values = c("#7876B1", "#20854E", "#003C67FF", "#BC3C29", "#767676FF")) +
  scale_fill_manual(values = c("#7876B1", "#20854E", "#003C67FF", "#BC3C29", "#767676FF") %>%
                      colorspace::lighten(0.5),
                    guide = 'none') +
  scale_x_continuous(expand = expansion(add = c(1, 1))) +
  scale_y_continuous(breaks = -log10(c(1e-3, 5e-2, 1e-2, 1e-1, 1)),
                     labels = c('0.001', '0.05', '0.01', '0.1', '1')) +
  ggprism::theme_prism(border = F, base_size = 12, base_fontface = 'plain') +
  theme(aspect.ratio = 0.8,
        panel.grid.major = element_line(color = 'grey30', linewidth = 0.25, linetype = 'dotted'),
        legend.title = element_text()) +
  labs(y = 'p_coxph') +
  coord_cartesian(clip = 'off')
dev.off()


# Figure S13b
pdf('FigS13b.pdf', height = 4, width = 5)
Prognosis_Results %>% 
  dplyr::mutate(Sig = ifelse(p_coxph < 0.05, 'Sig', 'Not_Sig') %>% 
                  forcats::fct_relevel('Sig', 'Not_Sig')) %>% 
  dplyr::nest_by() %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(RESULT_rowwise = table(data$SA_group, data$Sig) %>% 
                  rstatix::row_wise_fisher_test(detailed = T) %>% 
                  list(),
                RESULT_global = table(data$SA_group, data$Sig) %T>% 
                  {set.seed(42)} %>% 
                  rstatix::fisher_test(simulate.p.value = T, detailed = T) %>%
                  list()) %>% 
  dplyr::mutate(Plot = {data %>% 
      dplyr::count(Sig, SA_group) %>%
      dplyr::group_by(SA_group) %>% 
      dplyr::mutate(x_pos = cumsum(0.5 * cumsum(n))) %>% 
      ggplot() + 
      geom_bar(aes(y = SA_group, x = n, fill = SA_group, alpha = Sig %>% forcats::fct_rev()), 
               stat = 'identity', position = 'stack') +
      geom_text(aes(y = SA_group, x = x_pos, label = n, color = ifelse(Sig == 'Sig', 'white', 'black'))) +
      geom_text(aes(y = group, 
                    label = paste0('P = ', sprintf('%.2g', p.adj)), 
                    color = dplyr::case_when(p.adj < 0.05 & estimate > 1 ~ 'red', 
                                             p.adj < 0.05 & estimate < 1 ~ 'blue', 
                                             T ~ 'black')),
                x = 90, size = 3, hjust = 0, vjust = 3,
                data = RESULT_rowwise) +
      
      scale_x_continuous(expand = expansion(c(0, 0.02))) +
      scale_y_discrete(limits = rev, position = 'right') +
      scale_fill_manual(values = c(
        "#7876B1", "#20854E", "#003C67FF", "#BC3C29", "#767676FF"), 
        guide = 'none') +
      scale_alpha_manual(values = c(0.2, 1), guide = 'none') +
      scale_color_identity() +
      ggprism::theme_prism(border = T, base_fontface = 'plain', base_size = 12) + 
      theme(title = element_text(size = 12),
            axis.title = element_text(size = 12)) +
      coord_fixed(ratio = 20, clip = 'off') +
      
      labs(x = 'Number of prognostic features\n', y = NULL, 
           title = paste0(RESULT_global$method, ': p = ', RESULT_global$p))
  } %>% 
    list()) %>% 
  dplyr::pull(Plot)
dev.off()


# Figure S13c
pdf('FigS13c.pdf', height = 4, width = 5)
TherapeuticBenefits_Results %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(Group_Benefits = Details %>% 
                  dplyr::filter(p_logrank < 0.05) %>%
                  dplyr::pull(Group) %>% 
                  list()) %>% 
  dplyr::mutate(Group_Benefits = ifelse(length(Group_Benefits) == 0, NA, Group_Benefits)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(Sig = ifelse((p_coxph_Interaction < 0.05 & Group_Benefits == 1), 'Sig', 'Not_Sig') %>% 
                  forcats::fct_relevel('Sig', 'Not_Sig')) %>% 
  dplyr::nest_by() %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(RESULT_rowwise = table(data$SA_group, data$Sig) %>% 
                  rstatix::row_wise_fisher_test(detailed = T) %>% 
                  list(),
                RESULT_global = table(data$SA_group, data$Sig) %T>% 
                  {set.seed(42)} %>% 
                  rstatix::fisher_test(simulate.p.value = T, detailed = T) %>%
                  list()) %>% 
  dplyr::mutate(Plot = {data %>% 
      dplyr::count(Sig, SA_group) %>%
      dplyr::group_by(SA_group) %>% 
      dplyr::mutate(x_pos = cumsum(0.5 * cumsum(n))) %>% 
      ggplot() + 
      geom_bar(aes(y = SA_group, x = n, fill = SA_group, alpha = Sig %>% forcats::fct_rev()), 
               stat = 'identity', position = 'stack') +
      geom_text(aes(y = SA_group, x = x_pos, label = n, color = ifelse(Sig == 'Sig', 'white', 'black'))) +
      geom_text(aes(y = group, 
                    label = paste0('P = ', sprintf('%.2g', p.adj)), 
                    color = dplyr::case_when(p.adj < 0.05 & estimate > 1 ~ 'red', 
                                             p.adj < 0.05 & estimate < 1 ~ 'blue', 
                                             T ~ 'black')),
                x = 90, size = 3, hjust = 0, vjust = 3,
                data = RESULT_rowwise) +
      
      scale_x_continuous(expand = expansion(c(0, 0.02))) +
      scale_y_discrete(limits = rev, position = 'right') +
      scale_fill_manual(values = c("#7876B1", "#20854E", "#003C67FF", "#BC3C29", "#767676FF"), 
                        guide = 'none') +
      scale_alpha_manual(values = c(0.2, 1), guide = 'none') +
      scale_color_identity() +
      ggprism::theme_prism(border = T, base_fontface = 'plain', base_size = 12) + 
      theme(title = element_text(size = 12),
            axis.title = element_text(size = 12)) +
      coord_fixed(ratio = 20, clip = 'off') +
      
      labs(x = 'Number of features predictive of\ntherapeutic response', y = NULL, 
           title = paste0(RESULT_global$method, ': p = ', RESULT_global$p))
  } %>% 
    list()) %>% 
  dplyr::pull(Plot)
dev.off()

