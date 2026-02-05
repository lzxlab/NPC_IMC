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


# Figure S3a
pdf('FigS3a.pdf', width = 14, height = 12)
HM_CTxChannel <- MT_CT %>% 
  dplyr::group_by(CT2, CT4) %>% 
  dplyr::summarize(Count = dplyr::n(), 
                   across(LEVELS$Channel, .fns = ~ mean(.)))
HM_CTxChannel %>% 
  tibble::column_to_rownames('CT4') %>% 
  .[, c(-1,-2)] %>% 
  as.matrix() %>% 
  ComplexHeatmap::Heatmap(name = 'Channel\nExpression\n(Z-score)', cluster_rows = F, cluster_columns = F, 
                          circlize::colorRamp2(c(-1.5, 0, 1.5), c("#046C9A", "white", "#CB2314")),
                          width = ncol(.) * unit(6, 'mm'), 
                          height = nrow(.) * unit(6, 'mm'),
                          row_title = NULL, 
                          column_title = NULL, 
                          row_names_side = 'left',
                          border = T,
                          column_split = c(3, 4, 8, 10, 9) %>% rep(letters[1:length(.)], .),
                          row_split = HM_CTxChannel$CT2,
                          right_annotation = rowAnnotation(
                            Count = anno_barplot(HM_CTxChannel$Count, add_numbers = T, border = F, 
                                                 ylim = c(0, 300000), 
                                                 width = unit(2.4, 'cm'))))
dev.off()


# Figure S3b
pdf('FigS3b.pdf', height = 9, width = 26)
MT_CT %>%
    Helper2(UNIT_ID = 'Pat_ID', TERM_A = 'CT4', TERM_B = 'CT0',
            REL2A = F, REL2B = T) %>%
    dplyr::mutate(Value = ifelse(Prop4Filter < 0.01, 0, Value)) %>%
    dplyr::group_by(Term_A) %>%
    dplyr::mutate(ndis = dplyr::n_distinct(Value)) %>%
    dplyr::filter(ndis > 1) %>%
    dplyr::ungroup() %>% 
    dplyr::transmute(Pat_ID, Feature = Term_A, Value) %>% 
    dplyr::nest_by() %>% 
    tidyr::expand_grid(Group2Compare = c('Age', 'Sex', 
                                         'T7th', 'N7th', 'TNM', 
                                         'EBV', 'TLS_Pat', 
                                         'OS', 'EFS', 'DMFS')) %>% 
    dplyr::left_join(ClinicalInfo_Pat_Binary %>%
                         dplyr::select('Pat_ID', 'Age', 'Sex', 
                                       'T7th', 'N7th', 'TNM', 
                                       'EBV', 'TLS_Pat', 
                                       'OS', 'EFS', 'DMFS') %>% 
                         tidyr::pivot_longer(cols = !Pat_ID, 
                                             names_to = 'Group2Compare', 
                                             values_to = 'Group') %>% 
                         dplyr::filter(!Group %>% is.na()) %>% 
                         dplyr::nest_by(Group2Compare, .key = 'Clinical')) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(Result = data %>% 
                      dplyr::left_join(Clinical, by = 'Pat_ID') %>% 
                      dplyr::filter(Group %in% c(0, 1)) %>% 
                      dplyr::group_by(Feature) %>% 
                      {dplyr::left_join(rstatix::cohens_d(., Value ~ Group),
                                        rstatix::wilcox_test(., Value ~ Group),
                                        by = dplyr::join_by(.y., group1, group2, Feature, n1, n2))} %>% 
                      list()) %>% 
    dplyr::select(Group2Compare, Result) %>% 
    tidyr::unnest(Result) %>% 
    dplyr::nest_by(.key = 'Result') %>% 
    dplyr::mutate(
        Result_2 = Result %>% 
            dplyr::filter(!Feature %>% str_detect('unknown'))%>% 
            dplyr::mutate(Feature = Feature %>% forcats::fct_relevel(LEVELS[['CT4']])) %>% 
            dplyr::mutate(Group2Compare = Group2Compare %>% 
                              forcats::fct_relevel('Age', 'Sex', 
                                                   'T7th', 'N7th', 'TNM', 
                                                   'EBV', 'TLS_Pat', 
                                                   'OS', 'EFS', 'DMFS') %>% 
                              forcats::fct_rev()) %>% 
            dplyr::mutate(color = Group2Compare %>% 
                              as.character() %>% dplyr::recode(
                                  !!!setNames(
                                      c("#be0032", '#f99379', 
                                        '#f38400', "#f6a600", '#f3c300', 
                                        '#8db600', '#008856', 
                                        '#0067a5', '#604e97', '#8d6ad2'
                                      ) %>% colorspace::lighten(amount = 0.4),
                                      c('Age', 'Sex', 
                                        'T7th', 'N7th', 'TNM', 
                                        'EBV', 'TLS_Pat', 
                                        'OS', 'EFS', 'DMFS')
                                  ))) %>% 
            dplyr::mutate(color1 = ifelse(effsize > 0,
                                          color,
                                          color %>% 
                                              colorspace::darken(amount = 0.6) %>% 
                                              colorspace::desaturate(0.3)))  %>%
            dplyr::mutate(color2 = ifelse(p < 0.05,
                                          color1,
                                          color1 %>% 
                                              colorspace::adjust_transparency(alpha = 0.5))) %>%
            dplyr::mutate(color2 = color2 %>% fct_inorder()) %>% 
            dplyr::mutate(border = ifelse(p  < 0.05, 'black', '#FFFFFF')) %>% 
            dplyr::arrange(Feature) %>% 
            list()
    ) %>% 
    dplyr::mutate(
        temp_legend = {Result_2 %>% 
                dplyr::mutate(Group = ifelse(effsize > 0, 0, 1)) %>% 
                dplyr::select(color1, Group, Group2Compare) %>% 
                dplyr::distinct() %>% 
                ggplot(aes(x = Group, y = Group2Compare, fill = color1)) + 
                geom_point(size = 10, shape = 21, stroke = 2) + 
                geom_text(aes(x = Group + 0.3, y = Group2Compare, label = Group), size = 6) +
                theme_void() +
                scale_x_continuous(expand = expansion(0, 0.2)) +
                scale_fill_identity() +
                coord_equal(ratio = 0.5)} %>% 
            list()
    ) %>% 
    dplyr::mutate(
        temp_plot = {Result_2 %>% 
                ggplot(aes(x = Feature, y = Group2Compare)) +
                geom_point(aes(color = border, fill = color2, size = p), shape = 21, stroke = 2) + 
                ggprism::theme_prism(base_fontface = 'plain', base_size = 16) +
                theme(axis.title.y = element_blank(), 
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                      legend.position = 'right',
                      panel.border = element_rect(colour = "black", fill = NA, linewidth = 2)
                ) +
                scale_y_discrete(position = 'right')+
                scale_color_identity() +
                scale_fill_identity() +
                scale_size('p', trans='log10', range = c(20, 3), breaks=c(0.001, 0.01, 0.05)) + 
                labs(x = NULL, y = NULL) +
                coord_equal()} %>% 
            list()
    ) %>% 
    dplyr::mutate(Plot = list(temp_plot + theme(legend.position = 'none'), 
                              temp_legend, 
                              temp_plot %>% ggpubr::get_legend()) %>% 
                      patchwork::wrap_plots(widths = c(10, 1, 1)) %>%
                      list()) %>% 
    dplyr::pull(Plot) %>% 
    dplyr::first()
dev.off()


