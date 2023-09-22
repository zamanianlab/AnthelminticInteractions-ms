# data wrangling/plotting
library(tidyverse)

# other plotting
library(ggbeeswarm)
library(ggtext)
library(cowplot)
library(ZamanianLabThemes)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(remotes)
# install.packages("remotes")
# remotes::install_github("ddsjoberg/bstfun")
# 
# install.packages("webshot2")
library(webshot2)

library(AICcmodavg)
# tables
library(gt)
library(gtExtras)

# stats
library(drc)
library(broom)

# misc
library(conflicted)
library(here)
library(DescTools)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

dr_data <- read_rds('/Users/elenagr/Desktop/dr_data_EGR.rds') %>% # when box : 'data_EGR/dr_data_EGR.rds'
  filter(case_when(
    metadata_date == '20211115' ~ FALSE,
    metadata_date == '20211118' ~ FALSE,
    metadata_date == '20211208' ~ FALSE,
    metadata_date == '20211216' ~ FALSE,
    metadata_date == '20220131' & metadata_plate %in% c('p01', 'p02', 'p03', 'p04') ~ TRUE,
    metadata_date == '20220203' & metadata_plate %in% c('p01', 'p02', 'p03', 'p04') ~ TRUE,
    metadata_date == '20220210' & metadata_plate %in% c('p01', 'p08', 'p15', 'p17', 'p23') ~ TRUE,
    metadata_date == '20220311' & metadata_plate %in% c('p01') ~ TRUE,
    metadata_date == '20220324' & metadata_plate %in% c('p05') ~ TRUE,
    metadata_date == '20220325' & metadata_plate %in% c('p01') ~ TRUE,
    metadata_date == '20220331' & metadata_plate %in% c('p01', 'p02', 'p05', 'p08')  ~ TRUE,
    #metadata_date == '20220331' & metadata_plate %in% c('p01', 'p02','p05', 'p08') & treatment %in% c('DMSO', 'Untreated', 'Ivermectin', 'Levamisole') ~ TRUE,
    metadata_date == '20220408' & metadata_plate %in% c('p01') ~ TRUE
  )) %>%
  # bad wells
  filter(case_when(
    
    TRUE ~ TRUE
  )) %>%
  mutate(
    assay_date = str_extract(plate, '202[1-2][0-9]{4}'),
    # fix conc metadata
    conc = case_when(
      treatment == 'Levamisole' & conc == '0.005uM' ~ '0.05uM',
      treatment == 'Levamisole' & conc == '6.9uM' ~ '6.95uM',
      treatment == 'AlbendazoleSulfoxide' & conc == '111.18uM' ~ '111.8uM',
      treatment == 'Ivermectin' & conc == '0.1865uM' ~ '0.1864uM',
      treatment == 'Ivermectin' & conc == '0.0013uM' ~ '0.00134uM',
      treatment == 'Ivermectin' & conc == '0.05uM' ~ '0.5uM',
      TRUE ~ conc),
    strains = case_when(
      assay_date == '20220408' ~ 'N2',
      TRUE ~ strains
    ),
    genotype = case_when(
      assay_date == '20220408' ~ 'N2',
      TRUE ~ genotype
    ),
  )

coef <- 1.5

dr_data <- dr_data %>% select(plate, well, row, col, conc, strains, treatment, metadata_date, metadata_plate,genotype, area_shape_major_axis_length)

outliers <- dr_data %>%
  group_by(treatment, strains, conc) %>%
  group_nest() %>%
  mutate(quant = map(data, ~ as.numeric(quantile(.x$area_shape_major_axis_length, c(0, 0.25, 0.5, 0.75, 1)))),
         iqr = unlist(map(quant, ~ diff(.x[c(2, 4)]))),
         top_cutoff = unlist(map(quant, ~ pluck(.x, 4))) + coef * iqr,
         bottom_cutoff = unlist(map(quant, ~ pluck(.x, 2))) - coef * iqr) %>% 
  select(treatment, strains, conc, contains('cutoff'))

trimmed_data <- dr_data %>% 
  left_join(outliers) %>% 
  mutate(outlier = case_when(
    area_shape_major_axis_length > top_cutoff | area_shape_major_axis_length < bottom_cutoff ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  filter(outlier == FALSE) %>% 
  select(-contains('cutoff'), -outlier) %>%
  rename('raw_length' = 'area_shape_major_axis_length') %>%
  mutate(
    sqrt_length = sqrt(raw_length)
  )

# DMSO_mean <- trimmed_data %>%
#   filter(treatment == 'DMSO') %>%
#   group_by(plate) %>% 
#   summarise(DMSO_mean = mean(raw_length, na.rm = FALSE)) %>%
#   ungroup()

# normalize by dividing by DMSO average
revised_data <- trimmed_data %>% 
  mutate(conc = case_when(
    treatment == 'DMSO' ~ 0.01,
    treatment == 'Untreated' ~ NA_real_,
    treatment != 'DMSO' ~ as.numeric(str_remove(conc, 'uM'))
  )) 
# %>%
#   left_join(DMSO_mean) %>%
#   mutate(raw = raw_length / DMSO_mean)


well_summary <- revised_data %>% 
  select(plate, well, conc, metadata_date, strains, treatment, genotype, raw_length) %>%
  group_by(metadata_date, plate, genotype, strains, well, treatment, conc) %>% 
  summarise(mean_raw_length = mean(raw_length)) %>% 
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  ))  %>%
  mutate(group = case_when(
    metadata_date %in% c('20220311', '20220324') ~ 1,
    metadata_date %in% c('20220131') & strains == 'PR672' ~ 4,
    metadata_date %in% c('20220131') & strains != 'PR672' ~ 1,
    metadata_date %in% c('20220203', '20220408', '20220325') ~ 2,
    metadata_date %in% c('20220210') ~ 3,
    metadata_date %in% c('20220331') ~ 4)) 


# plot all DMSO and Untreated for every strain and color by date -----------------
well_summary <- well_summary %>%
  mutate(drug_stock = case_when(
    metadata_date %in% c('20211115', '20211118', '20211208', '20211216') ~ 'before',
    metadata_date %in% c('20220131', '20220203', '20220210', '20220311', '20220324', '20220325', '20220331', '20220408') ~ 'after',
  ))

reps <- well_summary %>% 
  ungroup() %>% 
  filter(!treatment %in% c('DMSO', 'Untreated')) %>% 
  distinct(metadata_date, genotype, strains, treatment) %>% 
  arrange(desc(metadata_date)) %>% 
  group_by(genotype, strains, treatment) %>% 
  mutate(rep = row_number())
#5
well_summary <- well_summary %>% 
  left_join(reps)

### dose-responses

curves <- revised_data %>% 
  #left_join(reps) %>% 
  # use only the most recent 3 replicates (not just those that used the same drug preparation)
  filter(!treatment %in% c('DMSO', 'Untreated')) %>%
  mutate(group = case_when(
    metadata_date %in% c('20220311', '20220324') ~ 1,
    metadata_date %in% c('20220131') & strains == 'PR672' ~ 4,
    metadata_date %in% c('20220131') & strains != 'PR672' ~ 1,
    metadata_date %in% c('20220203', '20220408', '20220325') ~ 2,
    metadata_date %in% c('20220210') ~ 3,
    metadata_date %in% c('20220331') ~ 4))  %>%
  
  # uncomment below to fit to the summarized well data
  group_by(group, plate, well, treatment, strains, genotype, conc) %>%
  summarise(mean_raw_length = mean(raw_length)) %>% 
  ungroup() %>%
  group_nest(treatment, group) %>% 
  # rowwise() %>% 
  # fit the models using all the data but group by strain, give some parameters
  mutate(params = list(c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA))) %>% 
  mutate(drc_raw  = map2(data, params, ~ drm(.x$mean_raw_length ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_raw = map(drc_raw, glance)) %>%
  mutate(tidy_raw = map(drc_raw, tidy)) 

curves_all <- revised_data %>% 
  left_join(reps) %>% 
  # use only the most recent 3 replicates (not just those that used the same drug preparation)
  filter(!treatment %in% c('DMSO', 'Untreated')) %>%
  # filter(metadata_date %in% c('20220131', '20220311', '20220324')) %>%
  # filter(metadata_date %in% c('20220203', '20220408', '20220325')) %>%
  # filter(metadata_date %in% c('20220210')) %>%
  # filter(metadata_date %in% c('20220331')) %>%
  # uncomment below to fit to the summarized well data
  group_by(plate, treatment, strains, genotype, conc) %>%
  summarise(mean_raw_length = mean(raw_length)) %>% 
  ungroup() %>%
  group_nest(treatment) %>% 
  # rowwise() %>% 
  # fit the models using all the data but group by strain, give some parameters
  mutate(params = list(c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA))) %>% 
  mutate(drc_raw  = map2(data, params, ~ drm(.x$mean_raw_length ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_raw = map(drc_raw, glance)) %>%
  mutate(tidy_raw = map(drc_raw, tidy)) %>%
  mutate(group = '')


# get the ec50
ec50_raw <- curves %>%
  unnest(tidy_raw) %>%
  select(-drc_raw, -data, -glance_raw) %>%
  filter(term == 'e') %>%
  select(-term) %>%
  mutate(xmax = estimate + (0.434 * std.error / estimate),
         xmin = estimate - (0.434 * std.error / estimate)) %>%
  mutate(
    xmin = case_when(
      is.na(xmin) & treatment == 'Ivermectin' ~ 0.00025,
      is.na(xmin) & treatment == 'AlbendazoleSulfoxide' ~ 0.15,
      TRUE ~ xmin),
    xmax = case_when(
      is.na(xmax) & treatment == 'Ivermectin' ~ 1,
      is.na(xmax) & treatment == 'AlbendazoleSulfoxide' ~ 450,
      TRUE ~ xmax
    )
  ) %>%
  rename(curve_raw = curve, estimate_raw = estimate, std.error_raw = std.error,
         statistic_raw = statistic, p.value_raw = p.value, xmax_raw = xmax, xmin_raw = xmin) %>%
  group_by(curve_raw, treatment) %>%
  mutate(ec50_average = mean(estimate_raw)) %>%
  ungroup() %>%
  mutate(route = case_when(
    curve_raw == 'N2' ~ 'Wild type',
    curve_raw %in% c('PR672', 'SP1234') ~ 'Amphid',
    curve_raw %in% c('DC19', 'LC144') ~ 'Cuticle',
    curve_raw %in% c('DA453', 'AE501') ~ 'Digestive'
  )) %>%
  mutate( genotype = case_when(
    curve_raw == 'N2' ~ "Wild type",
    curve_raw == 'SP1234' ~ "*dyf-2(m160)*",
    curve_raw == 'DA453' ~ "*eat-2(ad453)*",
    curve_raw == 'LC144' ~ "*agmo-1(e3016)*",
    curve_raw == 'PR672' ~ "*che-1(p672)*",
    curve_raw == 'DC19' ~ "*bus-5(br19)*",
    curve_raw == 'AE501' ~ "*nhr-8(ok186)*"
  )) %>%
  mutate(route_group = paste(route, group)) %>%
  group_by(treatment, curve_raw) %>%
  mutate(p_val_avg = mean(p.value_raw), sd = sd(estimate_raw))

ec50_raw_all <- curves_all %>%
  unnest(tidy_raw) %>%
  select(-drc_raw, -data, -glance_raw) %>%
  filter(term == 'e') %>%
  select(-term) %>%
  mutate(xmax = estimate + (0.434 * std.error / estimate),
         xmin = estimate - (0.434 * std.error / estimate)) %>%
  mutate(
    xmin = case_when(
      is.na(xmin) & treatment == 'Ivermectin' ~ 0.00025,
      is.na(xmin) & treatment == 'AlbendazoleSulfoxide' ~ 0.15,
      TRUE ~ xmin),
    xmax = case_when(
      is.na(xmax) & treatment == 'Ivermectin' ~ 1,
      is.na(xmax) & treatment == 'AlbendazoleSulfoxide' ~ 450,
      TRUE ~ xmax
    )
  ) %>%
  rename(curve_raw = curve, estimate_raw = estimate, std.error_raw = std.error,
         statistic_raw = statistic, p.value_raw = p.value, xmax_raw = xmax, xmin_raw = xmin) %>%
  group_by(curve_raw, treatment) %>%
  mutate(ec50_average = mean(estimate_raw)) %>%
  ungroup() %>%
  mutate(route = case_when(
    curve_raw == 'N2' ~ 'Wild type',
    curve_raw %in% c('PR672', 'SP1234') ~ 'Amphid',
    curve_raw %in% c('DC19', 'LC144') ~ 'Cuticle',
    curve_raw %in% c('DA453', 'AE501') ~ 'Digestive'
  )) %>%
  mutate( genotype = case_when(
    curve_raw == 'N2' ~ "Wild type",
    curve_raw == 'SP1234' ~ "*dyf-2(m160)*",
    curve_raw == 'DA453' ~ "*eat-2(ad453)*",
    curve_raw == 'LC144' ~ "*agmo-1(e3016)*",
    curve_raw == 'PR672' ~ "*che-1(p672)*",
    curve_raw == 'DC19' ~ "*bus-5(br19)*",
    curve_raw == 'AE501' ~ "*nhr-8(ok186)*"
  )) %>%
  mutate(route_group = paste(route, group))
#7
terms_raw <- curves %>%
  unnest(tidy_raw) %>%
  select(treatment, curve, term, estimate, group) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(strains = curve, b_raw = b, c_raw = c, d_raw = d, e_raw = e)

terms_raw_all <- curves_all %>%
  unnest(tidy_raw) %>%
  select(treatment, curve, term, estimate, group) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(strains = curve, b_raw = b, c_raw = c, d_raw = d, e_raw = e)

# fit the model to draw a line
azs_conc <- tibble(treatment = 'AlbendazoleSulfoxide',
                   conc = exp(seq(log(min(filter(revised_data, treatment == 'AlbendazoleSulfoxide')$conc, na.rm = TRUE)), 
                                  log(max(filter(revised_data, treatment == 'AlbendazoleSulfoxide')$conc, na.rm = TRUE)), 
                                  length = 100)))

lev_conc <- tibble(treatment = 'Levamisole',
                   conc = exp(seq(log(min(filter(revised_data, treatment == 'Levamisole')$conc, na.rm = TRUE)), 
                                  log(max(filter(revised_data, treatment == 'Levamisole')$conc, na.rm = TRUE)), 
                                  length = 100)))

ivm_conc <- tibble(treatment = 'Ivermectin',
                   conc = exp(seq(log(min(filter(revised_data, treatment == 'Ivermectin')$conc, na.rm = TRUE)), 
                                  log(max(filter(revised_data, treatment == 'Ivermectin')$conc, na.rm = TRUE)), 
                                  length = 100)))

newdata <- expand.grid(
  strains = unique(dr_data$strains),
  treatment = c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin')) %>%
  left_join(bind_rows(azs_conc, lev_conc, ivm_conc)) %>% 
  select(treatment, conc, strains) %>% 
  left_join(distinct(select(dr_data, strains, genotype)))

predict <- newdata %>%
  left_join(terms_raw) %>%
  # mutate(
  #   c = case_when(
  #     treatment == 'AlbendazoleSulfoxide' ~ sqrt(0.6),
  #     treatment == 'Levamisole' ~ sqrt(0.5),
  #     TRUE ~ c
  #   )) %>%
  # manually predict using the LL.4 model and the tidied terms
  mutate(pred_raw = c_raw + ((d_raw - c_raw) / (1 + exp(b_raw*(log(conc) - log(e_raw)))))) %>%
  select(-b_raw, -c_raw, -d_raw) %>%
  select(treatment, strains, group, genotype, conc, pred_raw) %>%
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  )) %>%
  mutate(route_group = paste(route, group))

predict_all <- newdata %>%
  left_join(terms_raw_all) %>%
  # mutate(
  #   c = case_when(
  #     treatment == 'AlbendazoleSulfoxide' ~ sqrt(0.6),
  #     treatment == 'Levamisole' ~ sqrt(0.5),
  #     TRUE ~ c
  #   )) %>%
  # manually predict using the LL.4 model and the tidied terms
  mutate(pred_raw = c_raw + ((d_raw - c_raw) / (1 + exp(b_raw*(log(conc) - log(e_raw)))))) %>%
  select(-b_raw, -c_raw, -d_raw) %>%
  select(treatment, strains, genotype, conc, pred_raw) %>%
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  )) 

#make dataframes where the genotype for N2 is each mutation for the purpose of plotting vertical line reflecting N2 ec50s on graphs with mutant curves for easy comparison
N2_raw_ec50_dyf <- ec50_raw_all %>% filter(curve_raw == 'N2') %>% select(-genotype) %>% mutate(genotype = '*dyf-2(m160)*')
N2_raw_ec50_eat <- ec50_raw_all %>% filter(curve_raw == 'N2') %>% select(-genotype) %>% mutate(genotype = '*eat-2(ad453)*')
N2_raw_ec50_bus <- ec50_raw_all %>% filter(curve_raw == 'N2') %>% select(-genotype) %>% mutate(genotype = '*bus-5(br19)*')  
N2_raw_ec50_agm <- ec50_raw_all %>% filter(curve_raw == 'N2') %>% select(-genotype) %>% mutate(genotype = '*agmo-1(e3016)*')
N2_raw_ec50_nhr <- ec50_raw_all %>% filter(curve_raw == 'N2') %>% select(-genotype) %>% mutate(genotype = '*nhr-8(ok186)*')
N2_raw_ec50_che <- ec50_raw_all %>% filter(curve_raw == 'N2') %>% select(-genotype) %>% mutate(genotype = '*che-1(p672)*')  

# (curve_plot_WT <- revised_data %>%
#     filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin')) %>% 
#     filter(strains %in% c('N2')) %>% 
#     ggplot() +
#     geom_vline(data = ec50_raw %>% filter(curve_raw == 'N2'),
#                aes(xintercept = estimate_raw, color = route_group), linetype = 'dashed', size = 0.3, alpha = 0.75) +
#     geom_line(data = predict %>% filter(strains == 'N2', group == 1), aes(x = conc, y = pred_raw), color = '#622020', size = 0.3, alpha =0.75) +
#     geom_line(data = predict %>% filter(strains == 'N2', group == 2), aes(x = conc, y = pred_raw), color = '#b87333', size = 0.3, alpha =0.75) +
#     geom_line(data = predict %>% filter(strains == 'N2', group == 3), aes(x = conc, y = pred_raw), color = '#ebce75', size = 0.3, alpha =0.75) +
#     geom_line(data = predict %>% filter(strains == 'N2', group == 4), aes(x = conc, y = pred_raw), color = '#f3830d', size = 0.3, alpha =0.75) +
#     geom_quasirandom(data = well_summary %>% filter(strains == 'N2', treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
#                                                     plate != '20220210-p22-EJG_1187', group == 1), aes(x = conc, y = mean_raw_length), color = '#622020', size = 0.5, width = 0.1, alpha = 0.75) +
#     geom_quasirandom(data = well_summary %>% filter(strains == 'N2', treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
#                                                     plate != '20220210-p22-EJG_1187', group == 2), aes(x = conc, y = mean_raw_length), color = '#b87333', size = 0.5, width = 0.1, alpha = 0.75) +
#     geom_quasirandom(data = well_summary %>% filter(strains == 'N2', treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
#                                                     plate != '20220210-p22-EJG_1187', group == 3), aes(x = conc, y = mean_raw_length), color = '#ebce75', size = 0.5, width = 0.1, alpha = 0.75) +
#     geom_quasirandom(data = well_summary %>% filter(strains == 'N2', treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
#                                                     plate != '20220210-p22-EJG_1187', group == 4), aes(x = conc, y = mean_raw_length), color = '#f3830d', size = 0.5, width = 0.1, alpha = 0.75) +
#     geom_hline(yintercept = 0, color = 'black', size = 0.25) +
#     geom_vline(data = ec50_raw_all %>% filter(curve_raw == 'N2'), aes(xintercept = estimate_raw), linetype = 'dashed', color = 'black', size = 0.5) +
#     geom_line(data = predict_all %>% filter(strains == 'N2'), aes(x = conc, y = pred_raw), color = 'black', size = 0.75, alpha =1 ) +
#     
#     scale_x_log10(
#       breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
#       labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
#     annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
#                         short = unit(0.1, "cm"),
#                         mid = unit(0.15, "cm"),
#                         long = unit(0.2, "cm")) +
#     scale_y_continuous(limits = c(0, 300), expand = c(0, 0), breaks = seq(0, 300, 100)) +
#     #scale_color_manual(values = c('#94C58C',  '#89CFF0', '#AA98A9', 'orange',  'yellow', 'red') )+
#     scale_color_manual(values = c('#622020',  '#b87333', '#ebce75', '#f3830d')) +
#     scale_fill_manual(values = c('#622020',  '#b87333', '#ebce75', '#f3830d')) +
#     scale_shape(guide = 'none') +
#     coord_cartesian(clip = "off") +
#     
#     labs(x = 'Drug concentration (µM)', y = "Normalized length",
#          color = 'Route', fill = 'Route', shape = 'Replicate', tag = "D") +
#     facet_grid(rows = vars(strains), cols = vars(treatment), scales = 'free_x',
#                labeller = as_labeller(c(AlbendazoleSulfoxide = 'AZS',
#                                         Ivermectin = 'IVM',
#                                         Levamisole = 'LEV',
#                                         `N2` = 'Wild type'))) +
#     theme_nw2() +
#     ggtitle('Wild type') +
#     theme(
#       legend.position = 'none',
#       axis.text.x = element_markdown(angle = 0, hjust = 0.5),
#       axis.title.x = element_markdown(face = 'plain'),
#       axis.title.y = element_markdown(face = 'plain'),
#       plot.tag = element_text(size = 10),
#       strip.text = element_text(color = ('black'), size = 7),
#       plot.title = element_text(color = '#f3830d', hjust = 0.5, size = 15)) +
#     NULL)

(curve_plot_cuticle <- revised_data %>%
    filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin')) %>% 
    mutate(route = case_when(
      strains == 'N2' ~ 'Wild type',
      strains %in% c('PR672', 'SP1234') ~ 'Amphid',
      strains %in% c('DC19', 'LC144') ~ 'Cuticle',
      strains %in% c('DA453', 'AE501') ~ 'Digestive'
    )) %>% 
    # don't show N2 here
    filter(route %in% c('Cuticle')) %>% 
    ggplot() +
    geom_vline(data = ec50_raw %>% filter(curve_raw != 'N2', route == 'Cuticle'),
               aes(xintercept = estimate_raw, group = genotype, color = route_group), linetype = 'dashed', size = 0.3, alpha = 0.75) +
    geom_vline(data = N2_raw_ec50_bus, aes(xintercept = ec50_average), color = '#fac04e', linetype = 'solid', size = 0.3, alpha = 0.75) +
    geom_vline(data = N2_raw_ec50_agm, aes(xintercept = ec50_average), color = '#fac04e', linetype = 'solid', size = 0.3, alpha = 0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Cuticle', group == 1),
              aes(x = conc, group = genotype, y = pred_raw), color = '#94C58C', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Cuticle', group == 2),
              aes(x = conc, group = genotype, y = pred_raw), color = '#1A8828', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Cuticle', group == 3),
              aes(x = conc, group = genotype, y = pred_raw), color = '#094F26', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Cuticle', group == 4),
              aes(x = conc, group = genotype, y = pred_raw), color = '#2AAA8A', size = 0.3, alpha =0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Cuticle', group == 1), 
                     aes(x = conc, y = mean_raw_length), color = '#94C58C', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Cuticle', group == 2), 
                     aes(x = conc, y = mean_raw_length), color = '#1A8828', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Cuticle', group == 3), 
                     aes(x = conc, y = mean_raw_length), color = '#094F26', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Cuticle', group == 4), 
                     aes(x = conc, y = mean_raw_length), color = '#2AAA8A', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_hline(yintercept = 0, color = 'black', size = 0.25) +
    geom_vline(data = ec50_raw_all %>% filter(curve_raw != 'N2', route == 'Cuticle'),
               aes(xintercept = estimate_raw, group = genotype), linetype = 'dashed', color = 'black', size = 0.5) +
    geom_line(data = predict_all %>% filter(route == 'Cuticle'),
              aes(x = conc, group = genotype, y = pred_raw), color = 'black', size = 0.75, alpha =1 ) +
    
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
    annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scale_y_continuous(limits = c(0, 300), expand = c(0, 0), breaks = seq(0, 300, 50)) +
    #scale_color_manual(values = c('#94C58C',  '#89CFF0', '#AA98A9', 'orange',  'yellow', 'red') )+
    scale_color_manual(values = c('#94C58C',  '#1A8828', '#094F26', '#2AAA8A')) +
    scale_fill_manual(values = c('#94C58C',  '#1A8828', '#094F26', '#2AAA8A')) +
    scale_shape(guide = 'none') +
    coord_cartesian(clip = "off") +
    
    labs(x = 'Drug concentration (µM)', y = "Raw length",
         color = 'Route', fill = 'Route', shape = 'Replicate', tag = "B") +
    facet_grid(rows = vars(genotype), cols = vars(treatment), scales = 'free_x',
               labeller = as_labeller(c(AlbendazoleSulfoxide = 'AZS',
                                        Ivermectin = 'IVM',
                                        Levamisole = 'LEV',
                                        `*agmo-1(e3016)*` = '*agmo-1(e3016)*',
                                        `*bus-5(br19)*` = '*bus-5(br19)*'))) +
    theme_nw2() +
    ggtitle('Cuticle') +
    theme(
      legend.position = 'none',
      axis.text.x = element_markdown(angle = 0, hjust = 0.5),
      axis.title.x = element_markdown(face = 'plain'),
      axis.title.y = element_markdown(face = 'plain'),
      plot.tag = element_text(size = 10),
      strip.text = element_text(color = ('black'), size = 7),
      plot.title = element_text(color = '#1A8828', hjust = 0.5, size = 15)) +
    NULL)

(curve_plot_digestive <- revised_data %>%
    filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin')) %>% 
    mutate(route = case_when(
      strains == 'N2' ~ 'Wild type',
      strains %in% c('PR672', 'SP1234') ~ 'Amphid',
      strains %in% c('DC19', 'LC144') ~ 'Cuticle',
      strains %in% c('DA453', 'AE501') ~ 'Digestive'
    )) %>% 
    # don't show N2 here
    filter(route %in% c('Digestive')) %>% 
    ggplot() +
    geom_vline(data = ec50_raw %>% filter(curve_raw != 'N2', route == 'Digestive'),
               aes(xintercept = estimate_raw, group = genotype, color = route_group), linetype = 'dashed', size = 0.3, alpha = 0.75) +
    geom_vline(data = N2_raw_ec50_eat, aes(xintercept = ec50_average), color = '#fac04e', linetype = 'solid', size = 0.3, alpha = 0.75) +
    geom_vline(data = N2_raw_ec50_nhr, aes(xintercept = ec50_average), color = '#fac04e', linetype = 'solid', size = 0.3, alpha = 0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Digestive', group == 1),
              aes(x = conc, group = genotype, y = pred_raw), color = '#AA98A9', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Digestive', group == 2),
              aes(x = conc, group = genotype, y = pred_raw), color = '#800080', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Digestive', group == 3),
              aes(x = conc, group = genotype, y = pred_raw), color = '#702963', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Digestive', group == 4),
              aes(x = conc, group = genotype, y = pred_raw), color = '#CCCCFF', size = 0.3, alpha =0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Digestive', group == 1), 
                     aes(x = conc, y = mean_raw_length), color = '#AA98A9', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Digestive', group == 2), 
                     aes(x = conc, y = mean_raw_length), color = '#800080', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Digestive', group == 3), 
                     aes(x = conc, y = mean_raw_length), color = '#702963', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Digestive', group == 4), 
                     aes(x = conc, y = mean_raw_length), color = '#CCCCFF', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_hline(yintercept = 0, color = 'black', size = 0.25) +
    geom_line(data = predict_all %>% filter(route == 'Digestive'),
              aes(x = conc, group = genotype, y = pred_raw), color = 'black', size = 0.75, alpha =1 ) +
    geom_vline(data = ec50_raw_all %>% filter(curve_raw != 'N2', route == 'Digestive'),
               aes(xintercept = estimate_raw, group = genotype), linetype = 'dashed', color = 'black', size = 0.5, alpha = 1) +
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
    annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scale_y_continuous(limits = c(0, 300), expand = c(0, 0), breaks = seq(0, 300, 50)) +
    #scale_color_manual(values = c('#94C58C',  '#89CFF0', '#AA98A9', 'orange',  'yellow', 'red') )+
    scale_color_manual(values = c('#AA98A9', '#800080', '#702963', '#CCCCFF')) +
    scale_fill_manual(values = c('#AA98A9', '#800080', '#702963', '#CCCCFF')) +
    scale_shape(guide = 'none') +
    coord_cartesian(clip = "off") +
    
    labs(x = 'Drug concentration (µM)', y = "Raw length",
         color = 'Route', fill = 'Route', shape = 'Replicate', tag = "A") +
    facet_grid(rows = vars(genotype), cols = vars(treatment), scales = 'free_x',
               labeller = as_labeller(c(AlbendazoleSulfoxide = 'AZS',
                                        Ivermectin = 'IVM',
                                        Levamisole = 'LEV',
                                        `*eat-2(ad453)*` = '*eat-2(ad453)*',
                                        `*nhr-8(ok186)*` = '*nhr-8(ok186)*'))) +
    theme_nw2() +
    ggtitle('Digestive') +
    theme(
      legend.position = 'none',
      axis.text.x = element_markdown(angle = 0, hjust = 0.5),
      axis.title.x = element_markdown(face = 'plain'),
      axis.title.y = element_markdown(face = 'plain'),
      plot.tag = element_text(size = 10),
      strip.text = element_text(color = ('black'), size = 7),
      plot.title = element_text(color = '#702963', hjust = 0.5, size = 15)) +
    NULL)

(curve_plot_amphid <- revised_data %>%
    filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin')) %>% 
    mutate(route = case_when(
      strains == 'N2' ~ 'Wild type',
      strains %in% c('PR672', 'SP1234') ~ 'Amphid',
      strains %in% c('DC19', 'LC144') ~ 'Cuticle',
      strains %in% c('DA453', 'AE501') ~ 'Digestive'
    )) %>% 
    # don't show N2 here
    filter(route %in% c('Amphid')) %>% 
    ggplot() +
    geom_vline(data = ec50_raw %>% filter(curve_raw != 'N2', route == 'Amphid'),
               aes(xintercept = estimate_raw, group = genotype, color = route_group), linetype = 'dashed', size = 0.3, alpha = 0.75) +
    geom_vline(data = N2_raw_ec50_che, aes(xintercept = ec50_average), color = '#fac04e', linetype = 'solid', size = 0.3, alpha = 0.75) +
    geom_vline(data = N2_raw_ec50_dyf, aes(xintercept = ec50_average), color = '#fac04e', linetype = 'solid', size = 0.3, alpha = 0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Amphid', group == 2),
              aes(x = conc, group = genotype, y = pred_raw), color = '#4169E1', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Amphid', group == 3),
              aes(x = conc, group = genotype, y = pred_raw), color = '#191970', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Amphid', group == 4),
              aes(x = conc, group = genotype, y = pred_raw), color = '#0096FF', size = 0.3, alpha =0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Amphid', group == 2), 
                     aes(x = conc, y = mean_raw_length), color = '#4169E1', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Amphid', group == 3), 
                     aes(x = conc, y = mean_raw_length), color = '#191970', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Amphid', group == 4), 
                     aes(x = conc, y = mean_raw_length), color = '#0096FF', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_hline(yintercept = 0, color = 'black', size = 0.25) +
    geom_vline(data = ec50_raw_all %>% filter(curve_raw != 'N2', route == 'Amphid'),
               aes(xintercept = estimate_raw, group = genotype), linetype = 'dashed', color = 'black', size = 0.3, alpha = 1) +
    geom_line(data = predict_all %>% filter(route == 'Amphid'),
              aes(x = conc, group = genotype, y = pred_raw), color = 'black', size = 0.75, alpha =1 ) +
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
    annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scale_y_continuous(limits = c(0, 300), expand = c(0, 0), breaks = seq(0, 300, 50)) +
    #scale_color_manual(values = c('#94C58C',  '#89CFF0', '#AA98A9', 'orange',  'yellow', 'red') )+
    scale_color_manual(values = c('#4169E1', '#191970', '#0096FF')) +
    scale_fill_manual(values = c('#4169E1', '#191970', '#0096FF')) +
    scale_shape(guide = 'none') +
    coord_cartesian(clip = "off") +
    
    labs(x = 'Drug concentration (µM)', y = "Raw length",
         color = 'Route', fill = 'Route', shape = 'Replicate', tag = "C") +
    facet_grid(rows = vars(genotype), cols = vars(treatment), scales = 'free_x',
               labeller = as_labeller(c(AlbendazoleSulfoxide = 'AZS',
                                        Ivermectin = 'IVM',
                                        Levamisole = 'LEV',
                                        `*che-1(p672)*` = '*che-1(p672)*',
                                        `*dyf-2(m160)*` = '*dyf-2(m160)*'))) +
    theme_nw2() +
    ggtitle('Amphid') +
    theme(
      legend.position = 'none',
      axis.text.x = element_markdown(angle = 0, hjust = 0.5),
      axis.title.x = element_markdown(face = 'plain'),
      axis.title.y = element_markdown(face = 'plain'),
      plot.tag = element_text(size = 10),
      strip.text = element_text(color = ('black'), size = 7),
      plot.title = element_text(color = '#191970', hjust = 0.5, size = 15)) +
    NULL)

(curve_plots_mutants_raw <- plot_grid(curve_plot_digestive, curve_plot_cuticle, curve_plot_amphid, align = "h", rel_widths = c(1, 1, 1), nrow = 1))

(fig_2_raw_and_norm <- plot_grid(curve_plots_mutants_norm, curve_plots_mutants_raw, align = "v", rel_widths = c(1, 1), nrow = 2))

ggsave('/Users/elenagr/Desktop/Supp_Fig_3.tiff', curve_plots_mutants_raw, width = 300, height = 90, units = 'mm')
ggsave('/Users/elenagr/Desktop/Supp_Fig_3.png', curve_plots_mutants_raw, width = 300, height = 90, units = 'mm')

ggsave('/Users/elenagr/Desktop/Supp_Fig_Fig2Comp.pdf', fig_2_raw_and_norm, width = 550, height = 250, units = 'mm')
ggsave('/Users/elenagr/Desktop/Supp_Fig_Fig2Comp.png', fig_2_raw_and_norm, width = 550, height = 250, units = 'mm')



