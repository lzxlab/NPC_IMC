setwd("~/projects/NPC_SpatialProteomics/")
library(tidyverse)
library(data.table)
library(ComplexHeatmap)
source('0.Helpers.R')

MT_CT <- qs::qread('SingleCellDataFrame.qs')
LEVELS <- qs::qread('AnnotationLevels.qs')
ClinicalInfo_ROI <- data.table::fread('ClinicalInfo.csv')
ClinicalInfo_Pat <- ClinicalInfo_ROI %>% 
  dplyr::select(-ROI_ID, -TLS_ROI) %>% 
  dplyr::distinct()


PALETTE1 <- c('#005F73', 
              '#04707F', '#0A9396',
              '#9B9B7A', '#baa587', '#d9ae94', '#d08c60', '#b58463',
              '#E9D8A6', '#ECCE90', '#EFC479', '#F0B962', '#F0AF4A', '#EFA52F', '#EE9B00',
              'grey70') %>% 
  setNames(LEVELS$CT3)
PALETTE2 <- c('#6c757d', '#6A3D9A') %>% setNames(c('0', '1'))
PALETTE3 <- c('#046C9A', '#CB2314') %>% setNames(c('0', '1'))
HEIGHT1 <- unit(0.4, "cm")
HEIGHT2 <- unit(10, "cm")
WIDTH1 <- unit(0.15, "cm")



# Figure 1b
pdf('Fig1b.pdf', width = 14, height = 7)
HM_CTxChannel <- MT_CT %>% 
  dplyr::group_by(CT2, CT3) %>% 
  dplyr::summarize(Count = dplyr::n(), 
                   across(LEVELS$Channel[-c(2:3, 31:34, 19, 22, 26:30)], .fns = ~ mean(.)))

HM_CTxChannel %>% 
  tibble::column_to_rownames('CT3') %>% 
  .[, c(-1,-2)] %>% 
  as.matrix() %>% 
  ComplexHeatmap::Heatmap(name = 'Channel\nExpression\n(Z-score)',
                          cluster_rows = F, cluster_columns = F, 
                          circlize::colorRamp2(c(-1.5, 0, 1.5), c("#046C9A", "white", "#CB2314")),
                          width = ncol(.) * unit(6, 'mm'), 
                          height = nrow(.) * unit(6, 'mm'),
                          row_title = NULL, 
                          column_title = NULL, 
                          row_names_side = 'left',
                          border = T,
                          row_split = HM_CTxChannel$CT2,
                          right_annotation = rowAnnotation(
                            Count = anno_barplot(HM_CTxChannel$Count, add_numbers = T, border = F, 
                                                 ylim = c(0, 200000), 
                                                 width = unit(2.4, 'cm'))))
dev.off()


# Figure 1c
pdf('Fig1c.pdf', width = 15, height = 6)
MT_CT %>% 
  dplyr::filter(ROI_ID %in% c('T109_ROI1', 'T109_ROI2')) %>% 
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
    PLOT = {points_final %>% 
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
                          drop = F) +
        theme_void() +
        guides(fill = guide_legend(ncol = 1)) +
        theme(plot.margin = unit(c(0,0,0,0), "mm"))
    } %>% list()
  ) %>% 
  dplyr::pull(PLOT) %>% 
  purrr::walk(print)
dev.off()


# Figure 1d
pdf('Fig1d.pdf', height = 5, width = 27)
MT_CT %>%
  dplyr::filter(CT0 != 'Epi') %>% 
  Helper1('ROI_ID', 'CT3') %>% 
  tidyr::pivot_wider(id_cols = 'ROI_ID', names_from = 'Feature', values_from = 'Value', values_fill = 0) %>% 
  dplyr::left_join(ClinicalInfo_ROI[,.(ROI_ID, EFS, InductChemo)], by = 'ROI_ID') %>% 
  dplyr::arrange(EFS, InductChemo, desc(Endo + CAF + Granulo + Mono + M1 + M2 + DC)) %>% 
  dplyr::nest_by(EFS, InductChemo, .keep = T) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(Plot = data %>% 
                  {Heatmap(matrix = matrix(0, nrow = 0, ncol = nrow(.)), cluster_rows = F, cluster_columns = F,
                           row_title = NULL, height = 0, width = WIDTH1 * nrow(.), show_heatmap_legend = F,
                           top_annotation = HeatmapAnnotation(
                             gap = unit(1, "mm"),
                             '  ' = anno_simple(.$InductChemo, height = HEIGHT1, col = PALETTE2),
                             ' ' = anno_simple(.$EFS, height = HEIGHT1, col = PALETTE3),
                             '   ' = anno_barplot(dplyr::select(., LEVELS$CT3),
                                                  gp = gpar(fill = PALETTE1[LEVELS$CT3],
                                                            col = NA),
                                                  bar_width = 0.8, height = HEIGHT2, axis = F)
                           )
                  )} %>%
                  list()
  ) %>% 
  dplyr::select(InductChemo, Plot) %>%
  dplyr::nest_by(InductChemo) %>%
  dplyr::mutate(Plots = data$Plot %>%
                  purrr::reduce(`+`) %>%
                  draw(annotation_legend_list = list(
                    Legend(title = 'Event', labels =  c('0', '1'),
                           legend_gp = gpar(fill = PALETTE3)),
                    Legend(title = 'Induction\nChemotherapy', labels =  c('0', '1'),
                           legend_gp = gpar(fill = PALETTE2)),
                    Legend(title = 'Cell Type', labels = LEVELS$CT3,
                           legend_gp = gpar(fill = PALETTE1[LEVELS$CT3]))
                  )) %>%
                  list()) %>% 
  dplyr::pull(Plots)
dev.off()


# Figure 1e
ClinicalInfo_Pat_Binary <- ClinicalInfo_Pat %>% 
    dplyr::mutate(Age = ifelse(Age >= 45, 1, 0),
                  T7th = ifelse(T7th == '4', 1, 0),
                  N7th = ifelse(N7th == '3', 1, 0),
                  TNM = ifelse(TNM == '4a' | TNM == '4b', 1, 0),
                  TLS_Pat = ifelse(TLS_Pat > 0, 0, 1))


pdf('Fig1e.pdf', height = 8, width = 14)
MT_CT %>%
  Helper2(UNIT_ID = 'Pat_ID', TERM_A = 'CT3', TERM_B = 'CT0',
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
      dplyr::mutate(Feature = Feature %>% forcats::fct_relevel(LEVELS[['CT3']])) %>% 
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
                  patchwork::wrap_plots(widths = c(0.4, 0.1, 0.05)) %>%
                  list()) %>% 
  dplyr::pull(Plot) %>% 
  dplyr::first() %>% 
  print()
dev.off()
