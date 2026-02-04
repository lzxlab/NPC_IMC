setwd('~/projects/NPC_SpatialProteomics/')
library(tidyverse)
library(data.table)


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



# Figure S8a
pdf('FigS8a.pdf', width = 16, height = 3)
IMC_CD8_T %>%
  dplyr::mutate(APT_Bin = APT %>% Hmisc::cut2(cuts = 0:20/ 20)) %>%
  dplyr::mutate(APT_Bin = dplyr::case_when(
    APT == 0 & APT_Bin == '[0.00,0.05)' ~ '0',
    APT > 0 & APT_Bin == '[0.00,0.05)' ~ '(0.00,0.05)',
    T ~ as.character(APT_Bin)) %>%
      forcats::fct_relevel('0')) %>%
  dplyr::filter(SA != 'Noise') %>%
  dplyr::select(ROI_ID, Cell_ID, Layer, APT, APT_Bin, c('GranzymeB', 'CD45RO', 'ICOS', 'LAG3', 'Tim3', 'Caspase3cleaved')) %>%
  tidyr::pivot_longer(cols = c('GranzymeB', 'CD45RO', 'ICOS', 'LAG3', 'Tim3', 'Caspase3cleaved'),
                      names_to = 'Channel', values_to = 'Expression') %>%
  dplyr::mutate(Channel = Channel %>% 
                  forcats::fct_relevel(c('GranzymeB', 'CD45RO', 'ICOS', 'LAG3', 'Tim3', 'Caspase3cleaved'))) %>%
  dplyr::nest_by(Channel) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Prediction =
                  speedglm::speedglm(data =  data,
                                     formula = Expression ~ splines::ns(APT, df = 4) + splines::ns(Layer, df = 4) + splines::ns(APT, df = 4) * splines::ns(Layer, df = 4),
                                     model = F, y = F, fitted = T) %>%
                  predict(type = "response") %>%
                  list()) %>%
  tidyr::unnest(cols = c(data, Prediction)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Channel, Layer, APT_Bin) %>%
  dplyr::summarize(MeanPrediction = mean(Prediction)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Channel) %>%
  dplyr::mutate(ScaledEstimate = MeanPrediction %>% scales::rescale(c(0, 1))) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(Channel) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = APT_Bin, y = Layer)) +
  facet_grid(. ~ Channel, axes = 'all', switch = 'y') +
  geom_tile(aes(fill = ScaledEstimate)) +
  scale_fill_gradientn(colors = circlize::colorRamp2(colors = c('#046C9A', '#FFFFFF', '#CB2314'),
                                                     breaks = 1:3)(seq(1, 3, length.out = 100))) +    #dc040c
  scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), expand = expansion(c(0, 0))) +
  scale_x_discrete(breaks = c('0', '[0.20,0.25)', '[0.45,0.50)', '[0.70,0.75)', '[0.95,1.00]'),
                   labels = c('', '[0.00, 0.25)', '[0.25, 0.50)', '[0.50, 0.75)', '[0.75, 1.00]')
  ) +
  ggprism::theme_prism(base_fontface = 'plain', base_line_size = 0.5) +
  theme(strip.placement = 'outside', axis.text.x = element_text(size = 6, hjust = 1.05, vjust = 1.15)) +
  coord_fixed(ratio = 1)
dev.off()



# Figure S8b
IMC_CD8_T <- MT_CT %>%
  dplyr::filter(CT3 == 'CD8_T') %>% 
  dplyr::mutate(APT = PD1 %>%
                  dplyr::dense_rank() %>%
                  scales::rescale(to = c(0, 1))) %>% 
  dplyr::filter(!is.na(Layer))

pdf('FigS8b.pdf', width = 16, height = 6)
IMC_CD8_T %>%
  dplyr::mutate(APT_Bin = APT %>% Hmisc::cut2(cuts = 0:20/ 20)) %>%
  dplyr::mutate(APT_Bin = dplyr::case_when(
    APT == 0 & APT_Bin == '[0.00,0.05)' ~ '0',
    APT > 0 & APT_Bin == '[0.00,0.05)' ~ '(0.00,0.05)',
    T ~ as.character(APT_Bin)) %>%
      forcats::fct_relevel('0')) %>%
  dplyr::filter(SA != 'Noise') %>%
  dplyr::select(ROI_ID, Cell_ID, Layer, APT, APT_Bin, c('GranzymeB', 'CD45RO', 'ICOS', 'LAG3', 'Tim3', 'Caspase3cleaved')) %>%
  tidyr::pivot_longer(cols = c('GranzymeB', 'CD45RO', 'ICOS', 'LAG3', 'Tim3', 'Caspase3cleaved'),
                      names_to = 'Channel', values_to = 'Expression') %>%
  dplyr::mutate(Channel = Channel %>% 
                  forcats::fct_relevel(c('GranzymeB', 'CD45RO', 'ICOS', 'LAG3', 'Tim3', 'Caspase3cleaved'))) %>%
  dplyr::left_join(ClinicalInfo_ROI, by = 'ROI_ID') %>%
  dplyr::mutate(DMFS = DMFS %>% as.factor()) %>%
  dplyr::filter(!is.na(DMFS)) %>%
  dplyr::nest_by(Channel, DMFS) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Prediction =
                  speedglm::speedglm(data =  data,
                                     formula = Expression ~ splines::ns(APT, df = 4) + splines::ns(Layer, df = 4) + splines::ns(APT, df = 4) * splines::ns(Layer, df = 4),
                                     model = F, y = F, fitted = T) %>%
                  predict(type = "response") %>%
                  list()) %>%
  tidyr::unnest(cols = c(data, Prediction)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Channel, DMFS, Layer, APT_Bin) %>%
  dplyr::summarize(MeanPrediction = mean(Prediction)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Channel) %>%
  dplyr::mutate(ScaledEstimate = MeanPrediction %>% scales::rescale(c(0, 1))) %>%
  dplyr::ungroup() %>% 
  dplyr::group_by(Channel, DMFS) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = APT_Bin, y = Layer)) +
  facet_grid(paste0('DM = ', DMFS) ~ Channel, axes = 'all', switch = 'y') +
  geom_tile(aes(fill = ScaledEstimate)) +
  scale_fill_gradientn(colors = circlize::colorRamp2(colors = c('#046C9A', '#FFFFFF', '#CB2314'),
                                                     breaks = 1:3)(seq(1, 3, length.out = 100))) +    #dc040c
  scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), expand = expansion(c(0, 0))) +
  scale_x_discrete(breaks = c('0', '[0.20,0.25)', '[0.45,0.50)', '[0.70,0.75)', '[0.95,1.00]'),
                   labels = c('', '[0.00, 0.25)', '[0.25, 0.50)', '[0.50, 0.75)', '[0.75, 1.00]')
  ) +
  ggprism::theme_prism(base_fontface = 'plain', base_line_size = 0.5) +
  theme(strip.placement = 'outside', axis.text.x = element_text(size = 6, hjust = 1.05, vjust = 1.15)) +
  coord_fixed(ratio = 1)
dev.off()


# Figure S8c
MT_unscaled_ICOS_Tim3 <- qs::qread('Unscaled_ICOS_Tim3.qs')
pdf('FigS8c.pdf', height = 4, width = 12)
IMC_CD8_T %>%
  dplyr::left_join(MT_unscaled_ICOS_Tim3, by = 'Cell_ID') %>% 
  dplyr::mutate(`ICOS/Tim3` = ifelse(ICOS_unscaled == 0, 0, ICOS_unscaled / Tim3_unscaled)) %>% 
  dplyr::select(Pat_ID, Layer, APT, c('Tim3', 'ICOS', 'ICOS/Tim3')) %>% 
  tidyr::pivot_longer(cols = c('Tim3', 'ICOS', 'ICOS/Tim3'), 
                      names_to = 'Channel', values_to = 'Expression') %>% 
  dplyr::mutate(Channel = Channel %>%
                  forcats::fct_relevel('Tim3', 'ICOS', 'ICOS/Tim3')) %>%
  dplyr::filter(APT >= 0.9 & Layer %in% c(0:4)) %>% 
  dplyr::group_by(Pat_ID, Channel) %>% 
  dplyr::summarize(Expression = mean(Expression)) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(ClinicalInfo_Pat, by = 'Pat_ID') %>% 
  dplyr::nest_by(Channel) %>% 
  dplyr::mutate(Cutpoint = data %>% 
                  survminer::surv_cutpoint(time = 'DMFStime', event = 'DMFS', variables = 'Expression', minprop = 0.3) %>%
                  .$cutpoint %>% 
                  .$cutpoint) %>% 
  dplyr::mutate(PLOT = data %>% 
                  dplyr::mutate(Group = ifelse(Expression > Cutpoint, 1, 0)) %>% 
                  {survminer::ggsurvplot(
                    data = ., 
                    fit = survival::survfit(survival::Surv(DMFStime, DMFS) ~ Group, data = .),
                    legend.title = Channel, xlab = 'Months', ylab = 'DMFS',
                    legend.labs = c('low', 'high'),
                    palette = c('#44978f', '#EE8026'),
                    ggtheme = ggprism::theme_prism(base_size = 10, base_fontface = 'plain') + 
                      theme(legend.position = 'right', 
                            legend.title = element_text()),
                    conf.int = F,  pval = T, risk.table = F, 
                  )} %>% 
                  .[[1]] %>% 
                  list()) %>% 
  dplyr::pull(PLOT) %>%
  patchwork::wrap_plots(nrow = 1)
dev.off()
