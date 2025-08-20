setwd("~/projects/NPC_SpatialProteomics/")
library(tidyverse)
library(data.table)
library(ComplexHeatmap)
source('0.Helpers.R')

MT_CT <- qs::qread('SingleCellDataFrame.qs')

LEVELS <- c('CT0', 'CT1', 'CT2', 'CT3', 'CT4') %>% 
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


# Figure S4a
pdf('FigS4a.pdf', height = 9, width = 12)
MT_CT %>% 
  Helper1('ROI_ID', 'CT3') %>%
  tidyr::pivot_wider(id_cols = 'ROI_ID', names_from = 'Feature', names_sort = T,
                     values_from = 'Value', values_fill = 0) %>% 
  tibble::column_to_rownames('ROI_ID') %>%
  rstatix::cor_test(vars = LEVELS$CT3[1:15], 
                    vars2 = LEVELS$CT3[1:15], 
                    method = 'spearman') %>% 
  dplyr::mutate(var1 = var1 %>% forcats::fct_inorder() %>% forcats::fct_rev(),
                var2 = var2 %>% forcats::fct_inorder()) %>% 
  dplyr::filter(as.integer(var1) + as.integer(var2)  <= max(as.integer(var2)) + 1) %>% 
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", '')) %T>% 
  {temp_HM <<- tidyr::pivot_wider(., id_cols = c('var1'), 
                                  names_from = 'var2', values_from = 'cor') %>% 
    tibble::column_to_rownames('var1') %>% 
    as.matrix()} %T>% 
  {temp_TEXT <<- tidyr::pivot_wider(., id_cols = c('var1'), 
                                    names_from = 'var2', values_from = 'p.signif', values_fill = '') %>% 
    tibble::column_to_rownames('var1') %>% 
    as.matrix()} %>% 
  {ComplexHeatmap::Heatmap(matrix = temp_HM,
                           name = paste0('Spearman\'s rho'), 
                           cluster_rows = F, cluster_columns = F, 
                           row_title = NULL, column_title = NULL, 
                           cluster_row_slices = F, cluster_column_slices = F, 
                           row_names_side = 'left',
                           col = circlize::colorRamp2(c(-1, 0, 1), c('#1965B0', '#EEEEEE', '#CB2314')), 
                           na_col = '#00000000',
                           
                           row_split = rownames(temp_HM) %>%
                             dplyr::recode(!!!setNames(
                               c(rep(LEVELS$CT2[1:4], times = c(1, 2, 5, 7))),
                               c(rownames(temp_HM)))) %>%
                             forcats::fct_inorder(),
                           
                           column_split = colnames(temp_HM) %>%
                             dplyr::recode(!!!setNames(
                               c(rep(LEVELS$CT2[1:4], times = c(1, 2, 5, 7))),
                               c(rownames(temp_HM)))) %>%
                             forcats::fct_inorder(),
                           
                           row_gap = unit(3, 'mm'),
                           column_gap = unit(3, 'mm'),
                           
                           width = ncol(temp_HM) * unit(6, 'mm'), 
                           height = nrow(temp_HM) * unit(6, 'mm'),
                           
                           cell_fun = function(j, i, x, y, w, h, col) {
                             grid.text(temp_TEXT[i, j], x, y, gp = gpar(fontsize = 8))
                           })}

dev.off()

CellProportionPerROI <- dplyr::left_join(
  MT_CT %>% 
    Helper1('ROI_ID', 'CT2') %>%
    tidyr::pivot_wider(id_cols = 'ROI_ID', names_from = 'Feature', names_sort = T,
                       values_from = 'Value', values_fill = 0), 
  MT_CT %>% 
    dplyr::mutate(CT4 = CT4 %>% forcats::fct_drop()) %>%
    Helper1('ROI_ID', 'CT4') %>%
    tidyr::pivot_wider(id_cols = 'ROI_ID', names_from = 'Feature', names_sort = T,
                       values_from = 'Value', values_fill = 0)
)
CellProportionPerROI_cor <- list(
  CellProportionPerROI %>% 
    tibble::column_to_rownames('ROI_ID') %>% 
    rstatix::cor_test(vars = LEVELS$CT4[8:34], 
                      vars2 = LEVELS$CT2[2:4], 
                      method = 'spearman') %>% 
    dplyr::mutate(var1 = var1 %>% forcats::fct_inorder() %>% forcats::fct_rev(),
                  var2 = var2 %>% forcats::fct_inorder()), 
  
  CellProportionPerROI %>% 
    tibble::column_to_rownames('ROI_ID') %>% 
    rstatix::cor_test(vars = LEVELS$CT4[8:34], 
                      vars2 = LEVELS$CT4[8:34], 
                      method = 'spearman') %>% 
    dplyr::mutate(var1 = var1 %>% forcats::fct_inorder() %>% forcats::fct_rev(),
                  var2 = var2 %>% forcats::fct_inorder()) %>% 
    dplyr::filter(as.integer(var1) + as.integer(var2)  <= max(as.integer(var2)) + 1)
) %>% 
  dplyr::bind_rows()



# Figure S4b
pdf('FigS4b.pdf', height = 9, width = 12)
CellProportionPerROI_cor %>% 
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                            symbols = c("***", "**", "*", '')) %T>% 
  {temp_HM <<- tidyr::pivot_wider(., id_cols = c('var1'), 
                                  names_from = 'var2', values_from = 'cor') %>% 
    tibble::column_to_rownames('var1') %>% 
    as.matrix()} %T>% 
  {temp_TEXT <<- tidyr::pivot_wider(., id_cols = c('var1'), 
                                    names_from = 'var2', values_from = 'p.signif', values_fill = '') %>% 
    tibble::column_to_rownames('var1') %>% 
    as.matrix()} %>% 
  {ComplexHeatmap::Heatmap(matrix = temp_HM,
                           name = paste0('Spearman\'s rho'), 
                           cluster_rows = F, cluster_columns = F, 
                           row_title = NULL, column_title = NULL, 
                           cluster_row_slices = F, cluster_column_slices = F, 
                           row_names_side = 'left',
                           col = circlize::colorRamp2(c(-1, 0, 1), c('#1965B0', '#EEEEEE', '#CB2314')), 
                           na_col = '#00000000',
                           
                           row_split = rownames(temp_HM) %>%
                             dplyr::recode(!!!setNames(
                               c(rep(LEVELS$CT2[2:4], times = c(4, 5, 18))),
                               c(LEVELS$CT4[8:34]))) %>%
                             forcats::fct_inorder(),
                           
                           column_split = colnames(temp_HM) %>%
                             dplyr::recode(!!!setNames(
                               c(c(1, 1, 1, 1), rep(LEVELS$CT2[2:4], times = c(4, 5, 18))),
                               c(LEVELS$CT2[1:4], LEVELS$CT4[8:34]))) %>%
                             forcats::fct_inorder(),
                           
                           row_gap = unit(3, 'mm'),
                           column_gap = unit(3, 'mm'),
                           
                           width = ncol(temp_HM) * unit(6, 'mm'), 
                           height = nrow(temp_HM) * unit(6, 'mm'),
                           
                           cell_fun = function(j, i, x, y, w, h, col) {
                             grid.text(temp_TEXT[i, j], x, y, gp = gpar(fontsize = 8))
                           })}
dev.off()
