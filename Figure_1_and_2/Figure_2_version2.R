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

dr_data <- read_rds(here('data_EGR/dr_data_EGR.rds')) %>% # when box : 'data_EGR/dr_data_EGR.rds'
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

DMSO_mean <- trimmed_data %>%
  filter(treatment == 'DMSO') %>%
  group_by(plate) %>% 
  summarise(DMSO_mean = mean(raw_length, na.rm = FALSE)) %>%
  ungroup()

# normalize by dividing by DMSO average
revised_data <- trimmed_data %>% 
  mutate(conc = case_when(
    treatment == 'DMSO' ~ 0.01,
    treatment == 'Untreated' ~ NA_real_,
    treatment != 'DMSO' ~ as.numeric(str_remove(conc, 'uM'))
  )) %>%
  left_join(DMSO_mean) %>%
  mutate(norm_DMSO = raw_length / DMSO_mean)
  

well_summary <- revised_data %>% 
  select(plate, well, conc, metadata_date, strains, treatment, genotype, raw_length, norm_DMSO) %>%
  group_by(metadata_date, plate, genotype, strains, well, treatment, conc) %>% 
  summarise(mean_raw_length = mean(raw_length),
            mean_norm_DMSO = mean(norm_DMSO)) %>% 
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
  summarise(mean_raw_length = mean(raw_length),
            mean_norm_DMSO = mean(norm_DMSO)) %>% 
  ungroup() %>%
  group_nest(treatment, group) %>% 
  # rowwise() %>% 
  # fit the models using all the data but group by strain, give some parameters
  mutate(params = list(c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA))) %>% 
  mutate(drc_norm_DMSO  = map2(data, params, ~ drm(.x$mean_norm_DMSO ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_norm_DMSO = map(drc_norm_DMSO, glance)) %>%
  mutate(tidy_norm_DMSO = map(drc_norm_DMSO, tidy)) 

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
  summarise(mean_raw_length = mean(raw_length),
            mean_norm_DMSO = mean(norm_DMSO)) %>% 
  ungroup() %>%
  group_nest(treatment) %>% 
  # rowwise() %>% 
  # fit the models using all the data but group by strain, give some parameters
  mutate(params = list(c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA))) %>% 
  mutate(drc_norm_DMSO  = map2(data, params, ~ drm(.x$mean_norm_DMSO ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_norm_DMSO = map(drc_norm_DMSO, glance)) %>%
  mutate(tidy_norm_DMSO = map(drc_norm_DMSO, tidy)) %>%
  mutate(group = '')


# get the ec50
ec50_norm_DMSO <- curves %>%
  unnest(tidy_norm_DMSO) %>%
  select(-drc_norm_DMSO, -data, -glance_norm_DMSO) %>%
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
  rename(curve_norm_DMSO = curve, estimate_norm_DMSO = estimate, std.error_norm_DMSO = std.error,
         statistic_norm_DMSO = statistic, p.value_norm_DMSO = p.value, xmax_norm_DMSO = xmax, xmin_norm_DMSO = xmin) %>%
  group_by(curve_norm_DMSO, treatment) %>%
  mutate(ec50_average = mean(estimate_norm_DMSO)) %>%
  ungroup() %>%
  mutate(route = case_when(
    curve_norm_DMSO == 'N2' ~ 'Wild type',
    curve_norm_DMSO %in% c('PR672', 'SP1234') ~ 'Amphid',
    curve_norm_DMSO %in% c('DC19', 'LC144') ~ 'Cuticle',
    curve_norm_DMSO %in% c('DA453', 'AE501') ~ 'Digestive'
  )) %>%
  mutate( genotype = case_when(
    curve_norm_DMSO == 'N2' ~ "Wild type",
    curve_norm_DMSO == 'SP1234' ~ "*dyf-2(m160)*",
    curve_norm_DMSO == 'DA453' ~ "*eat-2(ad453)*",
    curve_norm_DMSO == 'LC144' ~ "*agmo-1(e3016)*",
    curve_norm_DMSO == 'PR672' ~ "*che-1(p672)*",
    curve_norm_DMSO == 'DC19' ~ "*bus-5(br19)*",
    curve_norm_DMSO == 'AE501' ~ "*nhr-8(ok186)*"
  )) %>%
  mutate(route_group = paste(route, group)) %>%
  group_by(treatment, curve_norm_DMSO) %>%
  mutate(p_val_avg = mean(p.value_norm_DMSO), sd = sd(estimate_norm_DMSO))

ec50_norm_DMSO_all <- curves_all %>%
  unnest(tidy_norm_DMSO) %>%
  select(-drc_norm_DMSO, -data, -glance_norm_DMSO) %>%
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
  rename(curve_norm_DMSO = curve, estimate_norm_DMSO = estimate, std.error_norm_DMSO = std.error,
         statistic_norm_DMSO = statistic, p.value_norm_DMSO = p.value, xmax_norm_DMSO = xmax, xmin_norm_DMSO = xmin) %>%
  group_by(curve_norm_DMSO, treatment) %>%
  mutate(ec50_average = mean(estimate_norm_DMSO)) %>%
  ungroup() %>%
  mutate(route = case_when(
    curve_norm_DMSO == 'N2' ~ 'Wild type',
    curve_norm_DMSO %in% c('PR672', 'SP1234') ~ 'Amphid',
    curve_norm_DMSO %in% c('DC19', 'LC144') ~ 'Cuticle',
    curve_norm_DMSO %in% c('DA453', 'AE501') ~ 'Digestive'
  )) %>%
  mutate( genotype = case_when(
    curve_norm_DMSO == 'N2' ~ "Wild type",
    curve_norm_DMSO == 'SP1234' ~ "*dyf-2(m160)*",
    curve_norm_DMSO == 'DA453' ~ "*eat-2(ad453)*",
    curve_norm_DMSO == 'LC144' ~ "*agmo-1(e3016)*",
    curve_norm_DMSO == 'PR672' ~ "*che-1(p672)*",
    curve_norm_DMSO == 'DC19' ~ "*bus-5(br19)*",
    curve_norm_DMSO == 'AE501' ~ "*nhr-8(ok186)*"
  )) %>%
  mutate(route_group = paste(route, group))
#7
terms_norm_DMSO <- curves %>%
  unnest(tidy_norm_DMSO) %>%
  select(treatment, curve, term, estimate, group) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(strains = curve, b_norm_DMSO = b, c_norm_DMSO = c, d_norm_DMSO = d, e_norm_DMSO = e)

terms_norm_DMSO_all <- curves_all %>%
  unnest(tidy_norm_DMSO) %>%
  select(treatment, curve, term, estimate, group) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(strains = curve, b_norm_DMSO = b, c_norm_DMSO = c, d_norm_DMSO = d, e_norm_DMSO = e)

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
  left_join(terms_norm_DMSO) %>%
  # mutate(
  #   c = case_when(
  #     treatment == 'AlbendazoleSulfoxide' ~ sqrt(0.6),
  #     treatment == 'Levamisole' ~ sqrt(0.5),
  #     TRUE ~ c
  #   )) %>%
  # manually predict using the LL.4 model and the tidied terms
  mutate(pred_norm_DMSO = c_norm_DMSO + ((d_norm_DMSO - c_norm_DMSO) / (1 + exp(b_norm_DMSO*(log(conc) - log(e_norm_DMSO)))))) %>%
  select(-b_norm_DMSO, -c_norm_DMSO, -d_norm_DMSO) %>%
  select(treatment, strains, group, genotype, conc, pred_norm_DMSO) %>%
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  )) %>%
  mutate(route_group = paste(route, group))

predict_all <- newdata %>%
  left_join(terms_norm_DMSO_all) %>%
  # mutate(
  #   c = case_when(
  #     treatment == 'AlbendazoleSulfoxide' ~ sqrt(0.6),
  #     treatment == 'Levamisole' ~ sqrt(0.5),
  #     TRUE ~ c
  #   )) %>%
  # manually predict using the LL.4 model and the tidied terms
  mutate(pred_norm_DMSO = c_norm_DMSO + ((d_norm_DMSO - c_norm_DMSO) / (1 + exp(b_norm_DMSO*(log(conc) - log(e_norm_DMSO)))))) %>%
  select(-b_norm_DMSO, -c_norm_DMSO, -d_norm_DMSO) %>%
  select(treatment, strains, genotype, conc, pred_norm_DMSO) %>%
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  )) 
 

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
    geom_vline(data = ec50_norm_DMSO %>% filter(curve_norm_DMSO != 'N2', route == 'Amphid'),
               aes(xintercept = estimate_norm_DMSO, group = genotype, color = route_group), linetype = 'dashed', size = 0.3, alpha = 0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Amphid', group == 2),
              aes(x = conc, group = genotype, y = pred_norm_DMSO), color = '#4169E1', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Amphid', group == 3),
              aes(x = conc, group = genotype, y = pred_norm_DMSO), color = '#191970', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Amphid', group == 4),
              aes(x = conc, group = genotype, y = pred_norm_DMSO), color = '#0096FF', size = 0.3, alpha =0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Amphid', group == 2), 
                     aes(x = conc, y = mean_norm_DMSO), color = '#4169E1', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Amphid', group == 3), 
                     aes(x = conc, y = mean_norm_DMSO), color = '#191970', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Amphid', group == 4), 
                     aes(x = conc, y = mean_norm_DMSO), color = '#0096FF', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_hline(yintercept = 0, color = 'black', size = 0.25) +
     geom_vline(data = ec50_norm_DMSO_all %>% filter(curve_norm_DMSO != 'N2', route == 'Amphid'),
                aes(xintercept = estimate_norm_DMSO, group = genotype), linetype = 'dashed', color = 'black', size = 0.3, alpha = 1) +
     geom_line(data = predict_all %>% filter(route == 'Amphid'),
               aes(x = conc, group = genotype, y = pred_norm_DMSO), color = 'black', size = 0.75, alpha =1 ) +
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
    annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scale_y_continuous(limits = c(0, 1.45), expand = c(0, 0), breaks = seq(0, 1.2, 0.2)) +
    #scale_color_manual(values = c('#94C58C',  '#89CFF0', '#AA98A9', 'orange',  'yellow', 'red') )+
    scale_color_manual(values = c('#4169E1', '#191970', '#0096FF')) +
    scale_fill_manual(values = c('#4169E1', '#191970', '#0096FF')) +
    scale_shape(guide = 'none') +
    coord_cartesian(clip = "off") +
    
    labs(x = 'Drug concentration (µM)', y = "Normalized length",
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
    geom_vline(data = ec50_norm_DMSO %>% filter(curve_norm_DMSO != 'N2', route == 'Cuticle'),
               aes(xintercept = estimate_norm_DMSO, group = genotype, color = route_group), linetype = 'dashed', size = 0.3, alpha = 0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Cuticle', group == 1),
             aes(x = conc, group = genotype, y = pred_norm_DMSO), color = '#94C58C', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Cuticle', group == 2),
              aes(x = conc, group = genotype, y = pred_norm_DMSO), color = '#1A8828', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Cuticle', group == 3),
              aes(x = conc, group = genotype, y = pred_norm_DMSO), color = '#094F26', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Cuticle', group == 4),
              aes(x = conc, group = genotype, y = pred_norm_DMSO), color = '#2AAA8A', size = 0.3, alpha =0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Cuticle', group == 1), 
                     aes(x = conc, y = mean_norm_DMSO), color = '#94C58C', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Cuticle', group == 2), 
                     aes(x = conc, y = mean_norm_DMSO), color = '#1A8828', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Cuticle', group == 3), 
                     aes(x = conc, y = mean_norm_DMSO), color = '#094F26', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
                                                    plate != '20220210-p22-EJG_1187', route == 'Cuticle', group == 4), 
                     aes(x = conc, y = mean_norm_DMSO), color = '#2AAA8A', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_hline(yintercept = 0, color = 'black', size = 0.25) +
    geom_vline(data = ec50_norm_DMSO_all %>% filter(curve_norm_DMSO != 'N2', route == 'Cuticle'),
               aes(xintercept = estimate_norm_DMSO, group = genotype), linetype = 'dashed', color = 'black', size = 0.5) +
    geom_line(data = predict_all %>% filter(route == 'Cuticle'),
              aes(x = conc, group = genotype, y = pred_norm_DMSO), color = 'black', size = 0.75, alpha =1 ) +
    
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
    annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scale_y_continuous(limits = c(0, 1.45), expand = c(0, 0), breaks = seq(0, 1.2, 0.2)) +
    #scale_color_manual(values = c('#94C58C',  '#89CFF0', '#AA98A9', 'orange',  'yellow', 'red') )+
    scale_color_manual(values = c('#94C58C',  '#1A8828', '#094F26', '#2AAA8A')) +
    scale_fill_manual(values = c('#94C58C',  '#1A8828', '#094F26', '#2AAA8A')) +
    scale_shape(guide = 'none') +
    coord_cartesian(clip = "off") +
    
    labs(x = 'Drug concentration (µM)', y = "Normalized length",
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
    geom_vline(data = ec50_norm_DMSO %>% filter(curve_norm_DMSO != 'N2', route == 'Digestive'),
               aes(xintercept = estimate_norm_DMSO, group = genotype, color = route_group), linetype = 'dashed', size = 0.3, alpha = 0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Digestive', group == 1),
              aes(x = conc, group = genotype, y = pred_norm_DMSO), color = '#AA98A9', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Digestive', group == 2),
              aes(x = conc, group = genotype, y = pred_norm_DMSO), color = '#800080', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Digestive', group == 3),
              aes(x = conc, group = genotype, y = pred_norm_DMSO), color = '#702963', size = 0.3, alpha =0.75) +
    geom_line(data = predict %>% filter(strains != 'N2', route == 'Digestive', group == 4),
              aes(x = conc, group = genotype, y = pred_norm_DMSO), color = '#CCCCFF', size = 0.3, alpha =0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
        plate != '20220210-p22-EJG_1187', route == 'Digestive', group == 1), 
        aes(x = conc, y = mean_norm_DMSO), color = '#AA98A9', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
        plate != '20220210-p22-EJG_1187', route == 'Digestive', group == 2), 
        aes(x = conc, y = mean_norm_DMSO), color = '#800080', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
        plate != '20220210-p22-EJG_1187', route == 'Digestive', group == 3), 
        aes(x = conc, y = mean_norm_DMSO), color = '#702963', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
        plate != '20220210-p22-EJG_1187', route == 'Digestive', group == 4), 
        aes(x = conc, y = mean_norm_DMSO), color = '#CCCCFF', size = 0.5, width = 0.1, alpha = 0.75) +
    geom_hline(yintercept = 0, color = 'black', size = 0.25) +
    geom_line(data = predict_all %>% filter(route == 'Digestive'),
              aes(x = conc, group = genotype, y = pred_norm_DMSO), color = 'black', size = 0.75, alpha =1 ) +
    geom_vline(data = ec50_norm_DMSO_all %>% filter(curve_norm_DMSO != 'N2', route == 'Digestive'),
               aes(xintercept = estimate_norm_DMSO, group = genotype), linetype = 'dashed', color = 'black', size = 0.5, alpha = 1) +
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
    annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scale_y_continuous(limits = c(0, 1.45), expand = c(0, 0), breaks = seq(0, 1.2, 0.2)) +
    #scale_color_manual(values = c('#94C58C',  '#89CFF0', '#AA98A9', 'orange',  'yellow', 'red') )+
    scale_color_manual(values = c('#AA98A9', '#800080', '#702963', '#CCCCFF')) +
    scale_fill_manual(values = c('#AA98A9', '#800080', '#702963', '#CCCCFF')) +
    scale_shape(guide = 'none') +
    coord_cartesian(clip = "off") +
    
    labs(x = 'Drug concentration (µM)', y = "Normalized length",
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

(curve_plots_combined <- plot_grid(curve_plot_digestive, curve_plot_cuticle, curve_plot_amphid, align = "h", rel_widths = c(1, 1, 1), nrow = 1))

ggsave(here('plots_EGR/Figure_2/Figure_2.pdf'), curve_plots_combined, width = 300, height = 90, units = 'mm')
ggsave(here('plots_EGR/Figure_2/Figure_2.png'), curve_plots_combined, width = 300, height = 90, units = 'mm')

ggsave(here('/Users/elenagarncarz/Desktop/Figure_2_updated.tiff'), curve_plots_combined, width = 300, height = 90, units = 'mm', bg = "white")
ggsave(here('/Users/elenagarncarz/Desktop/Figure_2_updated.pdf'), curve_plots_combined, width = 225, height = 100, units = 'mm', bg = "white")
ggsave(here('/Users/elenagarncarz/Desktop/Figure_2_updated.png'), curve_plots_combined, width = 225, height = 100, units = 'mm', bg = "white")
# #get all EC values
# AZS_ECs <- ED(curves_all[[4]][[1]], c(10,50,90), interval="delta", type="relative")
# AZS_ECs <- as.data.frame(AZS_ECs) %>%
#   mutate(treatment = 'AZS') 
# 
# LEV_ECs <- ED(curves_all[[4]][[3]], c(10,50,90), interval="delta", type="relative")
# LEV_ECs <- as.data.frame(LEV_ECs) %>%
#   mutate(treatment = 'LEV') 
# 
# IVM_ECs <- ED(curves_all[[4]][[2]], c(10,50,90), interval="delta", type="relative")
# IVM_ECs <- as.data.frame(IVM_ECs) %>%
#   mutate(treatment = 'IVM') 
# 
# col_1 <-  c('AE501_10', 'AE501_50', 'AE501_90', 'DA453_10', 'DA453_50', 'DA453_90', 'DC19_10', 'DC19_50', 'DC19_90', 
#                            'LC144_10', 'LC144_50', 'LC144_90', 'N2_10', 'N2_50', 'N2_90')
# 
# col_2 <- c('AE501_10', 'AE501_50', 'AE501_90', 'DA453_10', 'DA453_50', 'DA453_90', 'LC144_10', 'LC144_50', 'LC144_90', 
#            'N2_10', 'N2_50', 'N2_90', 'PR672_10', 'PR672_50', 'PR672_90', 'SP1234_10', 'SP1234_50', 'SP1234_90')
#   
# col_3 <- c('AE501_10', 'AE501_50', 'AE501_90', 'LC144_10', 'LC144_50', 'LC144_90', 'N2_10', 'N2_50', 'N2_90', 
#            'PR672_10', 'PR672_50', 'PR672_90', 'SP1234_10', 'SP1234_50', 'SP1234_90')
#   
# col_4 <- c('DA453_10', 'DA453_50', 'DA453_90', 'DC19_10', 'DC19_50', 'DC19_90', 'PR672_10', 'PR672_50', 'PR672_90', 'SP1234_10', 'SP1234_50', 'SP1234_90')
# 
# AZS_ECs_1 <- ED(curves[[5]][[1]], c(10,50,90), interval="delta", type="relative")
# AZS_ECs_1 <- as.data.frame(AZS_ECs_1) %>%
#   mutate(treatment = 'AZS', group = 1)%>%
#   mutate(strains_ec = col_1)
# 
# AZS_ECs_2 <- ED(curves[[5]][[2]], c(10,50,90), interval="delta", type="relative")
# AZS_ECs_2 <- as.data.frame(AZS_ECs_2) %>%
#   mutate(treatment = 'AZS', group = 2) %>%
#   mutate(strains_ec = col_2)
# 
# AZS_ECs_3 <- ED(curves[[5]][[3]], c(10,50,90), interval="delta", type="relative")
# AZS_ECs_3 <- as.data.frame(AZS_ECs_3) %>%
#   mutate(treatment = 'AZS', group = 3) %>%
#   mutate(strains_ec = col_3)
# 
# AZS_ECs_4 <- ED(curves[[5]][[4]], c(10,50,90), interval="delta", type="relative")
# AZS_ECs_4 <- as.data.frame(AZS_ECs_4) %>%
#   mutate(treatment = 'AZS', group = 4) %>%
#   mutate(strains_ec = col_4)
# 
# col_4 <- c('DA453_10', 'DA453_50', 'DA453_90', 'DC19_10', 'DC19_50', 'DC19_90', 'LC144_10', 'LC144_50', 'LC144_90', 'N2_10', 'N2_50', 'N2_90', 'PR672_10', 'PR672_50', 'PR672_90', 'SP1234_10', 'SP1234_50', 'SP1234_90')
# 
# LEV_ECs_1 <- ED(curves[[5]][[5]], c(10,50,90), interval="delta", type="relative")
# LEV_ECs_1 <- as.data.frame(LEV_ECs_1) %>%
#   mutate(treatment = 'LEV', group = 1) %>%
#   mutate(strains_ec = col_1)
# 
# LEV_ECs_2 <- ED(curves[[5]][[6]], c(10,50,90), interval="delta", type="relative")
# LEV_ECs_2 <- as.data.frame(LEV_ECs_2) %>%
#   mutate(treatment = 'LEV', group = 2) %>%
#   mutate(strains_ec = col_2)
# 
# LEV_ECs_3 <- ED(curves[[5]][[7]], c(10,50,90), interval="delta", type="relative")
# LEV_ECs_3 <- as.data.frame(LEV_ECs_3) %>%
#   mutate(treatment = 'LEV', group = 3) %>%
#   mutate(strains_ec = col_3)

# LEV_ECs_4 <- ED(curves[[5]][[8]], c(10,50,90), interval="delta", type="relative")
# LEV_ECs_4 <- as.data.frame(LEV_ECs_4) %>%
#   mutate(treatment = 'LEV', group = 4) %>%
#   mutate(strains_ec = col_4)
# 
# col_4 <- c('DA453_10', 'DA453_50', 'DA453_90', 'DC19_10', 'DC19_50', 'DC19_90', 'N2_10', 'N2_50', 'N2_90', 'PR672_10', 'PR672_50', 'PR672_90', 'SP1234_10', 'SP1234_50', 'SP1234_90')
# 
# IVM_ECs_1 <- ED(curves[[5]][[9]], c(10,50,90), interval="delta", type="relative")
# IVM_ECs_1 <- as.data.frame(IVM_ECs_1) %>%
#   mutate(treatment = 'IVM', group = 1) %>%
#   mutate(strains_ec = col_1)
# 
# IVM_ECs_2 <- ED(curves[[5]][[10]], c(10,50,90), interval="delta", type="relative")
# IVM_ECs_2 <- as.data.frame(IVM_ECs_2) %>%
#   mutate(treatment = 'IVM', group = 2) %>%
#   mutate(strains_ec = col_2)
# 
# IVM_ECs_3 <- ED(curves[[5]][[11]], c(10,50,90), interval="delta", type="relative")
# IVM_ECs_3 <- as.data.frame(IVM_ECs_3) %>%
#   mutate(treatment = 'IVM', group = 3) %>%
#   mutate(strains_ec = col_3)
# 
# IVM_ECs_4 <- ED(curves[[5]][[12]], c(10,50,90), interval="delta", type="relative")
# IVM_ECs_4 <- as.data.frame(IVM_ECs_4) %>%
#   mutate(treatment = 'IVM', group = 4) %>%
#   mutate(strains_ec = col_4)
# 
# 
# ECs <- rbind(AZS_ECs_1, AZS_ECs_2, AZS_ECs_3, AZS_ECs_4, LEV_ECs_1,  LEV_ECs_2, LEV_ECs_3, LEV_ECs_4, IVM_ECs_1, IVM_ECs_2, IVM_ECs_3, IVM_ECs_4) %>%
#   group_by(strains_ec, treatment) %>%
#   mutate(mean_ec = mean(Estimate),
#          sd_ec = sd(Estimate))
# 
# ec50s_for_table <- ec50_norm_DMSO %>%
#   select(treatment, curve_norm_DMSO, ec50_average) %>%
#   unique()
# 
# N2_only_1 <- revised_data %>%
#   filter(strains == 'N2') %>%
#   select(-plate, -metadata_date, -metadata_plate) %>%
#   mutate(plate = '20230302-p01-EJG_3001', metadata_date = '20230301', metadata_plate = 'p01')
# 
# N2_only_2 <- revised_data %>%
#   filter(strains == 'N2') %>%
#   select(-plate, -metadata_date, -metadata_plate) %>%
#   mutate(plate = '20230302-p02-EJG_3002', metadata_date = '20230302', metadata_plate = 'p01')
# 
# N2_only_3 <- revised_data %>%
#   filter(strains == 'N2') %>%
#   select(-plate, -metadata_date, -metadata_plate) %>%
#   mutate(plate = '20230303-p03-EJG_3003', metadata_date = '20230303', metadata_plate = 'p01')
# 
# N2_only_4 <- revised_data %>%
#   filter(strains == 'N2') %>%
#   select(-plate, -metadata_date, -metadata_plate) %>%
#   mutate(plate = '20230304-p04-EJG_3004', metadata_date = '20230304', metadata_plate = 'p01')
# 
# revised_data_noN2 <- revised_data %>%
#   filter(strains != 'N2') 
# 
# revised_data_combinedN2 <- rbind(revised_data_noN2, N2_only_1, N2_only_2, N2_only_3, N2_only_4)
# 
# curves_pv_calc <- revised_data_combinedN2 %>% 
#   #left_join(reps) %>% 
#   # use only the most recent 3 replicates (not just those that used the same drug preparation)
#   filter(!treatment %in% c('DMSO', 'Untreated')) %>%
#   mutate(group = case_when(
#     metadata_date %in% c('20220311', '20220324', '20220131') ~ 1,
#     metadata_date %in% c('20220203', '20220408', '20220325') & strains != 'AE501' ~ 2,
#     metadata_date %in% c('20220203') & strains == 'AE501' ~ 4,
#     metadata_date %in% c('20220210') ~ 3,
#     metadata_date %in% c('20220331') ~ 4,
#     metadata_date %in% c('20230301') ~ 1,
#     metadata_date %in% c('20230302') ~ 2,
#     metadata_date %in% c('20230303') ~ 3,
#     metadata_date %in% c('20230304') ~ 4))  %>%
#   
#   # uncomment below to fit to the summarized well data
#   group_by(group, plate, well, treatment, strains, genotype, conc) %>%
#   summarise(mean_raw_length = mean(raw_length),
#             mean_norm_DMSO = mean(norm_DMSO)) %>% 
#   ungroup() %>%
#   group_nest(treatment, group) %>% 
#   # rowwise() %>% 
#   # fit the models using all the data but group by strain, give some parameters
#   mutate(params = list(c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA))) %>% 
#   mutate(drc_norm_DMSO  = map2(data, params, ~ drm(.x$mean_norm_DMSO ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
#   mutate(glance_norm_DMSO = map(drc_norm_DMSO, glance)) %>%
#   mutate(tidy_norm_DMSO = map(drc_norm_DMSO, tidy)) 
# 
# N2_only_5 <- revised_data %>%
#   filter(strains == 'N2', treatment == 'Ivermectin') %>%
#   select(norm_DMSO)
# 
# AE501_only <- revised_data %>%
#   filter(strains == 'PR672', treatment == 'Ivermectin') %>%
#   select(norm_DMSO) 
# 
# modelFit(curves_all[[4]][[1]])

# p_values_curves_AZS_1 <- curves_pv_calc[[7]][[1]] %>% filter(term == 'e') %>% mutate(treatment = 'AZS')
# p_values_curves_AZS_2 <- curves_pv_calc[[7]][[2]] %>% filter(term == 'e') %>% mutate(treatment = 'AZS')
# p_values_curves_AZS_3 <- curves_pv_calc[[7]][[3]] %>% filter(term == 'e') %>% mutate(treatment = 'AZS')
# p_values_curves_AZS_4 <- curves_pv_calc[[7]][[4]] %>% filter(term == 'e') %>% mutate(treatment = 'AZS')
# p_values_curves_IVM_1 <- curves_pv_calc[[7]][[5]] %>% filter(term == 'e') %>% mutate(treatment = 'IVM')
# p_values_curves_IVM_2 <- curves_pv_calc[[7]][[6]] %>% filter(term == 'e') %>% mutate(treatment = 'IVM')
# p_values_curves_IVM_3 <- curves_pv_calc[[7]][[7]] %>% filter(term == 'e') %>% mutate(treatment = 'IVM')
# p_values_curves_IVM_4 <- curves_pv_calc[[7]][[8]] %>% filter(term == 'e') %>% mutate(treatment = 'IVM')
# p_values_curves_LEV_1 <- curves_pv_calc[[7]][[9]] %>% filter(term == 'e') %>% mutate(treatment = 'LEV')
# p_values_curves_LEV_2 <- curves_pv_calc[[7]][[10]] %>% filter(term == 'e') %>% mutate(treatment = 'LEV')
# p_values_curves_LEV_3 <- curves_pv_calc[[7]][[11]] %>% filter(term == 'e') %>% mutate(treatment = 'LEV')
# p_values_curves_LEV_4 <- curves_pv_calc[[7]][[12]] %>% filter(term == 'e') %>% mutate(treatment = 'LEV')
# 
# p_values_curves <- rbind(p_values_curves_AZS_1, p_values_curves_AZS_2,p_values_curves_AZS_3, p_values_curves_AZS_4, p_values_curves_IVM_1, p_values_curves_IVM_2, 
#                          p_values_curves_IVM_3, p_values_curves_IVM_4, p_values_curves_LEV_1, p_values_curves_LEV_2, p_values_curves_LEV_3, p_values_curves_LEV_4) %>%
#   group_by(treatment, curve) %>%
#   mutate(pv_avg = mean(p.value)) %>%
#   select(treatment, curve, pv_avg) %>%
#   unique()

comps_norm_DMSO_all <- curves_all %>%                                      
  mutate(
    comp = map(drc_norm_DMSO, ~EDcomp(.x, percVec = c(50, 50))),
    comp_df = map(comp, ~as_tibble(.x, rownames = 'Comparison'))
  ) %>% 
  select(treatment, comp_df, group) %>% 
  unnest(cols = comp_df) %>% 
  filter(str_detect(Comparison, 'N2')) %>% 
  mutate(strains = str_remove_all(Comparison, str_c('N2', ':', '/', '50/50', sep = '|'))) %>% 
  select(strains, treatment, p.value =  `p-value`, group) %>% 
  group_by(treatment, strains) %>%
  mutate(avg_pv = mean(p.value))

comps_norm_DMSO <- curves_pv_calc %>%                                      
  mutate(
    comp = map(drc_norm_DMSO, ~EDcomp(.x, percVec = c(50, 50))),
    comp_df = map(comp, ~as_tibble(.x, rownames = 'Comparison'))
  ) %>% 
  select(treatment, comp_df) %>% 
  unnest(cols = comp_df) %>% 
  filter(str_detect(Comparison, 'N2')) %>% 
  mutate(strains = str_remove_all(Comparison, str_c('N2', ':', '/', '50/50', sep = '|'))) %>% 
  select(strains, treatment, p.value =  `p-value`) %>% 
  group_by(treatment, strains) %>%
  mutate(pv_avg = mean(p.value)) %>%
  select(-p.value) %>%
  unique()

ec50s <- ec50_norm_DMSO %>%
  select(treatment, curve_norm_DMSO, ec50_average, sd) %>%
  unique()

# revised_data_N2 <- revised_data %>%
#   filter(strains == 'N2', treatment== 'Levamisole')
# revised_data_AE501 <- revised_data %>%
#   filter(strains == 'AE501', treatment== 'Levamisole')
# revised_data_SP1234 <- revised_data %>%
#   filter(strains == 'SP1234', treatment== 'Levamisole')
# revised_data_DA453 <- revised_data %>%
#   filter(strains == 'DA453', treatment== 'Levamisole')
# revised_data_DC19 <- revised_data %>%
#   filter(strains == 'DC19', treatment== 'Levamisole')
# revised_data_PR672 <- revised_data %>%
#   filter(strains == 'PR672', treatment== 'Levamisole')
# revised_data_LC144<- revised_data %>%
#   filter(strains == 'LC144', treatment== 'Levamisole')
# 
# test <- drm(formula = revised_data_N2$norm_DMSO ~ revised_data_N2$conc, data = revised_data_N2, fct = LL.4())
# test1 <- drm(formula = revised_data_AE501$norm_DMSO ~ revised_data_AE501$conc, data = revised_data_AE501, fct = LL.4())
# test1 <- drm(formula = revised_data_SP1234$norm_DMSO ~ revised_data_SP1234$conc, data = revised_data_SP1234, fct = LL.4())
# test1 <- drm(formula = revised_data_DA453$norm_DMSO ~ revised_data_DA453$conc, data = revised_data_DA453, fct = LL.4())
# test1 <- drm(formula = revised_data_DC19$norm_DMSO ~ revised_data_DC19$conc, data = revised_data_DC19, fct = LL.4())
# test1 <- drm(formula = revised_data_PR672$norm_DMSO ~ revised_data_PR672$conc, data = revised_data_PR672, fct = LL.4())
# test1 <- drm(formula = revised_data_LC144$norm_DMSO ~ revised_data_LC144$conc, data = revised_data_LC144, fct = LL.4())
# 
# anova(test, test1)
# 
# (full_norm_DMSO <- plot_grid(curve_plot_norm_DMSO + remove_legend(), bstfun::as_ggplot(ec50_table_norm_DMSO),
#                              nrow = 1, rel_widths = c(1.5, 1),
#                              labels = 'auto'))
# 
# 
# ggsave(here('/Users/elenagr/Desktop/full_norm_DMSO.pdf'), full_norm_DMSO, width = 400, height = 450, units = 'mm')

#' (curve_plot_sqrt_norm_DMSO <- revised_data %>%
#'     filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin')) %>%
#'     mutate(route = case_when(
#'       strains == 'N2' ~ 'Wild type',
#'       strains %in% c('PR672', 'SP1234') ~ 'Amphid',
#'       strains %in% c('DC19', 'LC144') ~ 'Cuticle',
#'       strains %in% c('DA453', 'AE501') ~ 'Digestive'
#'     )) %>%
#'     # don't show N2 here
#'     filter(strains != 'N2') %>%
#'     ggplot() +
#'     geom_vline(data = ec50_pivoted_1 %>% filter(strains != 'N2', norm_method == 'estimate_sqrt_norm_DMSO'),
#'                aes(xintercept = estimate_value, group = genotype), linetype = 'dashed', color = '#212738', size = 0.3) +
#'     geom_vline(data = ec50_pivoted_2 %>% filter(strains != 'N2', norm_method == 'estimate_sqrt_norm_DMSO'),
#'                aes(xintercept = estimate_value, group = genotype), linetype = 'dashed', color = '#5C8492', size = 0.3) +
#'     geom_vline(data = ec50_pivoted_3 %>% filter(strains != 'N2', norm_method == 'estimate_sqrt_norm_DMSO'),
#'                aes(xintercept = estimate_value, group = genotype), linetype = 'dashed', color = '#498587', size = 0.3) +
#'     geom_vline(data = ec50_pivoted_4 %>% filter(strains != 'N2', norm_method == 'estimate_sqrt_norm_DMSO'),
#'                aes(xintercept = estimate_value, group = genotype), linetype = 'dashed', color = '#6D9DC5', size = 0.3) +
#'     geom_vline(data = ec50_average %>% filter(strains != 'N2', norm_method == 'estimate_sqrt_norm_DMSO'),
#'                aes(xintercept = estimate_avg, group = genotype), linetype = 'dashed', color = 'black', size = 0.6) +
#'     geom_line(data = predict_1 %>% filter(strains != 'N2'),
#'               aes(x = conc, group = genotype, y = pred_sqrt_norm_DMSO), #, color = group_route
#'               size = 0.3, alpha =0.5, color = '#212738') +
#'     geom_line(data = predict_2 %>% filter(strains != 'N2'),
#'               aes(x = conc, group = genotype, y = pred_sqrt_norm_DMSO), #, color = group_route
#'               size = 0.3, alpha =0.5, color = '#5C8492') +
#'     geom_line(data = predict_3 %>% filter(strains != 'N2', strains != 'DA453'),
#'               aes(x = conc, group = genotype, y = pred_sqrt_norm_DMSO), #, color = group_route
#'               size = 0.3, alpha =0.5, color = '#498587') +
#'     geom_line(data = predict_4 %>% filter(strains != 'N2'),
#'               aes(x = conc, group = genotype, y = pred_sqrt_norm_DMSO), #, color = group_route
#'               size = 0.3, alpha =0.5, color = '#6D9DC5') +
#'     geom_line(data = predict_all %>% filter(genotype != 'N2'),
#'               aes(x = conc, group = genotype, y = fit_line_sqrt_DMSO), #, color = group_route
#'               size = 0.75, color = 'black', alpha = 1) +
#'     
#'     geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
#'                                                     strains != 'N2',
#'                                                     plate != '20220210-p22-EJG_1187',
#'                                                     # rep <= 3
#'     ),
#'     aes(x = conc, y = mean_sqrt_norm_DMSO, color = plate),
#'     size = 0.5, width = 0.1, alpha = 0.6) +
#'     geom_hline(yintercept = 0.2, color = 'black', size = 0.25) +
#'     scale_x_log10(
#'       breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
#'       labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#'     annotation_logticks(sides = 'b', size = 0.25, outside = TRUE,
#'                         short = unit(0.1, "cm"),
#'                         mid = unit(0.15, "cm"),
#'                         long = unit(0.2, "cm")) +
#'     scale_y_continuous(limits = c(-0.2, 1.45), expand = c(0, 0), breaks = seq(-0.4, 1.2, 0.2)) +
#'     scale_color_manual(values = c('#212738', '#212738', '#212738', '#212738', '#212738', '#5C8492', '#5C8492', '#498587', '#5C8492','#498587', '#212738', '#5C8492', '#5C8492', '#498587', '#498587', '#5C8492', '#498587', 'blue')) +
#'     scale_fill_manual(values = c('#212738', '#212738', '#212738', '#212738', '#212738', '#5C8492', '#5C8492', '#498587', '#5C8492', '#498587', '#212738', '#5C8492', '#5C8492', '#498587', '#498587', '#5C8492', '#498587', 'blue')) +
#'     #scale_color_manual(values = c('#212738', '#622A7F', '#174636', '#6D9DC5', '#6F636B', '#2A7F62','#5C8492', '#82AC92', '#4985B7', '#C9BCF2', '#3DB88E', 'purple', 'orange', '#212738', '#622A7F', '#174636', '#6D9DC5', '#6F636B', '#2A7F62','#5C8492', '#82AC92', '#4985B7', '#C9BCF2', '#3DB88E', 'red', 'red', 'red', 'red', 'red', 'red', 'red','red', 'red')) +
#'     #scale_fill_manual(values = c('#212738', '#622A7F', '#174636', '#6D9DC5', '#6F636B', '#2A7F62','#5C8492', '#82AC92', '#4985B7', '#C9BCF2', '#3DB88E', 'purple', 'orange', '#212738', '#622A7F', '#174636', '#6D9DC5', '#6F636B', '#2A7F62','#5C8492', '#82AC92', '#4985B7', '#C9BCF2', '#3DB88E', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red')) +
#'     # scale_color_manual(values = c('black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'red', 'red', 'red')) +
#'     # scale_fill_manual(values = c('black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'red', 'red', 'red')) +
#'     #'#212738', '#212738', '#212738', '#212738', '#212738', '#5C8492', '#5C8492', '#498587', '#5C8492', '#498587', '#212738', '#5C8492', '#5C8492', '#498587', '#498587', '#5C8492', '#498587'
#'     scale_shape(guide = 'none') +
#'     coord_cartesian(clip = "off") +
#' 
#'     labs(x = '', y = "",
#'          color = 'Route', fill = 'Route', shape = 'Replicate') +
#'     facet_grid(rows = vars(genotype), cols = vars(treatment), scales = 'free_x',
#'                labeller = as_labeller(c(AlbendazoleSulfoxide = 'Albendazole sulfoxide',
#'                                         Ivermectin = 'Ivermectin',
#'                                         Levamisole = 'Levamisole',
#'                                         `*agmo-1(e3016)*` = '*agmo-1(e3016)*',
#'                                         `*bus-5(br19)*` = '*bus-5(br19)*',
#'                                         `*che-1(p672)*` = '*che-1(p672)*',
#'                                         `*dyf-2(m160)*` = '*dyf-2(m160)*',
#'                                         `*eat-2(ad453)*` = '*eat-2(ad453)*',
#'                                         `*nhr-8(ok186)*` = '*nhr-8(ok186)*'))) +
#'     theme_nw2() +
#'     theme(
#'       legend.position = 'none',
#'       axis.text.x = element_markdown(angle = 0, hjust = 0.5),
#'       axis.title.x = element_markdown(face = 'plain'),
#'       axis.title.y = element_markdown(face = 'plain')) +
#'     NULL)
#'  
#' 
#' comps_sqrt_norm_DMSO <- curves %>% 
#'   mutate(
#'   comp = map(drc_sqrt_norm_DMSO, ~EDcomp(.x, percVec = c(50, 50))),
#'   comp_df = map(comp, ~as_tibble(.x, rownames = 'Comparison'))
#'   ) %>% 
#'   select(treatment, comp_df) %>% 
#'   unnest(cols = comp_df) %>% 
#'   filter(str_detect(Comparison, 'N2')) %>% 
#'   mutate(strains = str_remove_all(Comparison, str_c('N2', ':', '/', '50/50', sep = '|'))) %>% 
#'   select(strains, treatment, p.value =  `p-value`) 
#' 
#' AE501_AZS_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "AE501" & treatment == 'AlbendazoleSulfoxide')
#' DA453_AZS_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "DA453" & treatment == 'AlbendazoleSulfoxide')
#' DC19_AZS_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "DC19" & treatment == 'AlbendazoleSulfoxide')
#' LC144_AZS_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "LC144" & treatment == 'AlbendazoleSulfoxide')
#' PR672_AZS_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "PR672" & treatment == 'AlbendazoleSulfoxide')
#' SP1234_AZS_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "SP1234" & treatment == 'AlbendazoleSulfoxide')
#' 
#' AE501_IVM_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "AE501" & treatment == 'Ivermectin')
#' DA453_IVM_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "DA453" & treatment == 'Ivermectin')
#' DC19_IVM_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "DC19" & treatment == 'Ivermectin')
#' LC144_IVM_pv_sqrt_norm <- subset(comps_norm_DMSO, strains == "LC144" & treatment == 'Ivermectin')
#' PR672_IVM_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "PR672" & treatment == 'Ivermectin')
#' SP1234_IVM_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "SP1234" & treatment == 'Ivermectin')
#' 
#' AE501_LEV_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "AE501" & treatment == 'Levamisole')
#' DA453_LEV_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "DA453" & treatment == 'Levamisole')
#' DC19_LEV_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "DC19" & treatment == 'Levamisole')
#' LC144_LEV_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "LC144" & treatment == 'Levamisole')
#' PR672_LEV_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "PR672" & treatment == 'Levamisole')
#' SP1234_LEV_pv_sqrt_norm <- subset(comps_sqrt_norm_DMSO, strains == "SP1234" & treatment == 'Levamisole')
#' 
#' ec50_sqrt_norm_DMSO_for_table <- curves %>%
#'   unnest(tidy_sqrt_norm_DMSO) %>%
#'   select(-drc_sqrt_norm_DMSO, -data, -glance_sqrt_norm_DMSO, -drc_norm_min, -tidy_norm_min, -glance_norm_min, -drc_norm_DMSO, -glance_norm_DMSO, -tidy_norm_DMSO, -drc_sqrt_norm_min,
#'          -glance_sqrt_norm_min, -tidy_sqrt_norm_min) %>%
#'   filter(term == 'e') %>%
#'   select(-term, -params) %>%
#'   mutate(xmax = estimate + (0.434 * std.error / estimate),
#'          xmin = estimate - (0.434 * std.error / estimate)) %>%
#'   mutate(
#'     xmin = case_when(
#'       is.na(xmin) & treatment == 'Ivermectin' ~ 0.00025,
#'       is.na(xmin) & treatment == 'AlbendazoleSulfoxide' ~ 0.15,
#'       TRUE ~ xmin),
#'     xmax = case_when(
#'       is.na(xmax) & treatment == 'Ivermectin' ~ 1,
#'       is.na(xmax) & treatment == 'AlbendazoleSulfoxide' ~ 450,
#'       TRUE ~ xmax
#'     )
#'   ) %>%
#'   select( -group) %>%
#'   rename(estimate_sqrt_norm_DMSO = estimate, std.error_sqrt_norm_DMSO = std.error,
#'          statistic_sqrt_norm_DMSO = statistic, p.value_sqrt_norm_DMSO = p.value, xmax_sqrt_norm_DMSO = xmax, xmin_sqrt_norm_DMSO = xmin)
#' 
#' tab_sqrt_norm_DMSO  <- ec50_sqrt_norm_DMSO_for_table %>% 
#'   rename(strains = curve) %>%
#'   select(-xmax_sqrt_norm_DMSO , -xmin_sqrt_norm_DMSO , -p.value_sqrt_norm_DMSO) %>% 
#'   mutate(p.value = case_when(
#'     strains == 'AE501' & treatment == 'AlbendazoleSulfoxide' ~ AE501_AZS_pv_sqrt_norm$p.value,
#'     strains == 'DA453' & treatment == 'AlbendazoleSulfoxide' ~ DA453_AZS_pv_sqrt_norm$p.value,
#'     strains == 'DC19' & treatment == 'AlbendazoleSulfoxide' ~ DC19_AZS_pv_sqrt_norm$p.value,
#'     strains == 'LC144' & treatment == 'AlbendazoleSulfoxide' ~ LC144_AZS_pv_sqrt_norm$p.value,
#'     strains == 'PR672' & treatment == 'AlbendazoleSulfoxide' ~ PR672_AZS_pv_sqrt_norm$p.value,
#'     strains == 'SP1234' & treatment == 'AlbendazoleSulfoxide' ~ SP1234_AZS_pv_sqrt_norm$p.value,
#'     
#'     strains == 'AE501' & treatment == 'Ivermectin' ~ AE501_IVM_pv_sqrt_norm$p.value,
#'     strains == 'DA453' & treatment == 'Ivermectin' ~ DA453_IVM_pv_sqrt_norm$p.value,
#'     strains == 'DC19' & treatment == 'Ivermectin' ~ DC19_IVM_pv_sqrt_norm$p.value,
#'     strains == 'LC144' & treatment == 'Ivermectin' ~ LC144_IVM_pv_sqrt_norm$p.value,
#'     strains == 'PR672' & treatment == 'Ivermectin' ~ PR672_IVM_pv_sqrt_norm$p.value,
#'     strains == 'SP1234' & treatment == 'Ivermectin' ~ SP1234_IVM_pv_sqrt_norm$p.value,
#'     
#'     strains == 'AE501' & treatment == 'Levamisole' ~ AE501_LEV_pv_sqrt_norm$p.value,
#'     strains == 'DA453' & treatment == 'Levamisole' ~ DA453_LEV_pv_sqrt_norm$p.value,
#'     strains == 'DC19' & treatment == 'Levamisole' ~ DC19_LEV_pv_sqrt_norm$p.value,
#'     strains == 'LC144' & treatment == 'Levamisole' ~ LC144_LEV_pv_sqrt_norm$p.value,
#'     strains == 'PR672' & treatment == 'Levamisole' ~ PR672_LEV_pv_sqrt_norm$p.value,
#'     strains == 'SP1234' & treatment == 'Levamisole' ~ SP1234_LEV_pv_sqrt_norm$p.value,
#'   )) %>%
#'   arrange(-estimate_sqrt_norm_DMSO)  %>% 
#'   arrange(fct_relevel(strains, 'N2')) %>% 
#'   select(-statistic_sqrt_norm_DMSO) %>% 
#'   mutate(p.value = case_when(
#'     p.value > 0.05 ~ as.character(round(p.value, 3)),
#'     p.value < 0.0001 ~ paste(as.character(round(p.value, 3)), '****', sep=""), 
#'     p.value < 0.001 ~ paste(as.character(round(p.value, 3)), '***', sep=""), 
#'     p.value < 0.01 ~ paste(as.character(round(p.value, 3)), '**', sep=""), 
#'     p.value < 0.05 ~ paste(as.character(round(p.value, 3)), '*', sep=""), 
#'     TRUE ~ as.character(p.value)
#'   )) %>%
#'   mutate(genotype = case_when(
#'     strains == 'N2' ~ 'Wild type',
#'     strains == 'AE501' ~ '*nhr-8(ok186)*',
#'     strains == 'DA453' ~ '*eat-2(ad453)*',
#'     strains == 'DC19' ~ '*bus-5(br19)*',
#'     strains == 'LC144' ~ '*agmo-1(e3016)*',
#'     strains == 'PR672' ~ '*che-1(p672)*',
#'     strains == 'SP1234' ~ '*dyf-2(m160)*'
#'   ))
#' 
#' ec50_table_sqrt_norm_DMSO  <- tab_sqrt_norm_DMSO  %>% 
#'   mutate(
#'     treatment = case_when(
#'       treatment == 'AlbendazoleSulfoxide' ~ 'Albendazole sulfoxide',
#'       TRUE ~ treatment),
#'     std_error = case_when(
#'       treatment == 'Ivermectin' ~ as.character(round(std.error_sqrt_norm_DMSO, 4)),
#'       TRUE ~ as.character(round(std.error_sqrt_norm_DMSO, 2))),
#'     # p.value = round(p.value, 4),
#'     estimate = case_when(
#'       treatment == 'Ivermectin' ~ as.character(round(estimate_sqrt_norm_DMSO, 4)),
#'       TRUE ~ as.character(round(estimate_sqrt_norm_DMSO, 2))
#'     )) %>% 
#'   select(-strains, -estimate_sqrt_norm_DMSO, -std.error_sqrt_norm_DMSO) %>% 
#'   gt(groupname_col = "treatment",
#'      # rowname_col = "genotype"
#'   ) %>%
#'   tab_options(
#'     table.align = "left",
#'     table.font.size = 10) %>%
#'   tab_header(title = "Dose-Response Estimates") %>%
#'   # set the font
#'   tab_style(
#'     style = cell_text(font = google_font("Helvetica")),
#'     locations = list(
#'       cells_column_labels(),
#'       cells_stub(),
#'       cells_title(groups = "title"),
#'       cells_body(columns = everything())
#'     )) %>%
#'   # rownames
#'   tab_style(
#'     style = cell_text(style = 'italic'),
#'     locations = cells_row_groups()
#'   ) %>%
#'   # columns
#'   cols_label(
#'     treatment = md("**Treatment**"),
#'     genotype = md("**Genotype**"),
#'     # strains = md("**Strain name**"),
#'     estimate = md("**EC<sub>50</sub>**"),
#'     std_error = md("**SE**"), 
#'     p.value = md("**_p_<br>(N2 == strain)**")) %>%
#'   # numeric columns
#'   # fmt_symbol_first(column = estimate, symbol = " µM", last_row_n = 8, gfont = "Helvetica", decimals = 2) %>%
#'   fmt_markdown(column = genotype, rows = everything()) %>% 
#'   # fmt_scientific(columns = p.value, rows = p.value < 0.001) %>% 
#'   cols_align(align = 'left', columns = everything()) %>%
#'   # color sig p-value
#'   # tab_style(
#'   #   style = cell_fill(color = 'grey', alpha = 0.5),
#'   #   locations = cells_body(
#'   #     rows = p.value < 0.05,
#'   #     columns = p.value)
#'   #   ) %>% 
#'   # color strains
#'   tab_style(
#'     style = cell_text(color = '#6D9DC5'),
#'     locations = cells_body(
#'       rows = genotype == '*agmo-1(e3016)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#6D9DC5'),
#'     locations = cells_body(
#'       rows = genotype == '*bus-5(br19)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#212738'),
#'     locations = cells_body(
#'       rows = genotype == '*che-1(p672)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#212738'),
#'     locations = cells_body(
#'       rows = genotype == '*dyf-2(m160)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#2A7F62'),
#'     locations = cells_body(
#'       rows = genotype == '*eat-2(ad453)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#2A7F62'),
#'     locations = cells_body(
#'       rows = genotype == '*nhr-8(ok186)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#A5243D'),
#'     locations = cells_body(
#'       rows = genotype == 'N2',
#'       columns = genotype)
#'   )
#' 
#' (full_sqrt_norm_DMSO <- plot_grid(curve_plot_sqrt_norm_DMSO + remove_legend(), bstfun::as_ggplot(ec50_table_sqrt_norm_DMSO),
#'                              nrow = 1, rel_widths = c(1.5, 1),
#'                              labels = 'auto'))
#' 
#' (curve_plot_norm_min <- revised_data %>%
#'     filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin')) %>%
#'     mutate(route = case_when(
#'       strains == 'N2' ~ 'Wild type',
#'       strains %in% c('PR672', 'SP1234') ~ 'Amphid',
#'       strains %in% c('DC19', 'LC144') ~ 'Cuticle',
#'       strains %in% c('DA453', 'AE501') ~ 'Digestive'
#'     )) %>%
#'     # don't show N2 here
#'     filter(strains != 'N2') %>%
#'     ggplot() +
#'     geom_vline(data = ec50_pivoted_1 %>% filter(strains != 'N2', norm_method == 'estimate_norm_min'),
#'                aes(xintercept = estimate_value, group = genotype), linetype = 'dashed', color = '#212738', size = 0.3) +
#'     geom_vline(data = ec50_pivoted_2 %>% filter(strains != 'N2', norm_method == 'estimate_norm_min'),
#'                aes(xintercept = estimate_value, group = genotype), linetype = 'dashed', color = '#5C8492', size = 0.3) +
#'     geom_vline(data = ec50_pivoted_3 %>% filter(strains != 'N2', norm_method == 'estimate_norm_min'),
#'                aes(xintercept = estimate_value, group = genotype), linetype = 'dashed', color = '#498587', size = 0.3) +
#'     geom_vline(data = ec50_pivoted_4 %>% filter(strains != 'N2', norm_method == 'estimate_norm_min'),
#'                aes(xintercept = estimate_value, group = genotype), linetype = 'dashed', color = '#6D9DC5', size = 0.3) +
#'     geom_vline(data = ec50_average %>% filter(strains != 'N2', norm_method == 'estimate_norm_min'),
#'                aes(xintercept = estimate_avg, group = genotype), linetype = 'dashed', color = 'black', size = 0.6) +
#'     geom_line(data = predict_1 %>% filter(strains != 'N2'),
#'               aes(x = conc, group = genotype, y = pred_norm_min), #, color = group_route
#'               size = 0.3, alpha =0.5, color = '#212738') +
#'     geom_line(data = predict_2 %>% filter(strains != 'N2'),
#'               aes(x = conc, group = genotype, y = pred_norm_min), #, color = group_route
#'               size = 0.3, alpha =0.5, color = '#5C8492') +
#'     geom_line(data = predict_3 %>% filter(strains != 'N2', strains != 'DA453'),
#'               aes(x = conc, group = genotype, y = pred_norm_min), #, color = group_route
#'               size = 0.3, alpha =0.5, color = '#498587') +
#'     geom_line(data = predict_4 %>% filter(strains != 'N2'),
#'               aes(x = conc, group = genotype, y = pred_norm_min), #, color = group_route
#'               size = 0.3, alpha =0.5, color = '#6D9DC5') +
#'     geom_line(data = predict_all %>% filter(genotype != 'N2'),
#'               aes(x = conc, group = genotype, y = fit_line_min), #, color = group_route
#'               size = 0.75, color = 'black', alpha = 1) +
#'     geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
#'                                                     strains != 'N2',
#'                                                     plate != '20220210-p22-EJG_1187',
#'                                                     # rep <= 3
#'     ),
#'     aes(x = conc, y = mean_norm_min, color = plate),
#'     size = 0.5, width = 0.1, alpha = 0.6) +
#'     geom_hline(yintercept = 0.2, color = 'black', size = 0.25) +
#'     scale_x_log10(
#'       breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
#'       labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#'     annotation_logticks(sides = 'b', size = 0.25, outside = TRUE,
#'                         short = unit(0.1, "cm"),
#'                         mid = unit(0.15, "cm"),
#'                         long = unit(0.2, "cm")) +
#'     scale_y_continuous(limits = c(-0.2, 1.45), expand = c(0, 0), breaks = seq(-0.5, 1.2, 0.2)) +
#'     scale_color_manual(values = c('#212738', '#212738', '#212738', '#212738', '#212738', '#5C8492', '#5C8492', '#498587', '#5C8492','#498587', '#212738', '#5C8492', '#5C8492', '#498587', '#498587', '#5C8492', '#498587', 'blue')) +
#'     scale_fill_manual(values = c('#212738', '#212738', '#212738', '#212738', '#212738', '#5C8492', '#5C8492', '#498587', '#5C8492', '#498587', '#212738', '#5C8492', '#5C8492', '#498587', '#498587', '#5C8492', '#498587', 'blue')) +
#'     #scale_color_manual(values = c('#212738', '#622A7F', '#174636', '#6D9DC5', '#6F636B', '#2A7F62','#5C8492', '#82AC92', '#4985B7', '#C9BCF2', '#3DB88E', 'purple', 'orange', '#212738', '#622A7F', '#174636', '#6D9DC5', '#6F636B', '#2A7F62','#5C8492', '#82AC92', '#4985B7', '#C9BCF2', '#3DB88E', 'red', 'red', 'red', 'red', 'red', 'red', 'red','red', 'red')) +
#'     #scale_fill_manual(values = c('#212738', '#622A7F', '#174636', '#6D9DC5', '#6F636B', '#2A7F62','#5C8492', '#82AC92', '#4985B7', '#C9BCF2', '#3DB88E', 'purple', 'orange', '#212738', '#622A7F', '#174636', '#6D9DC5', '#6F636B', '#2A7F62','#5C8492', '#82AC92', '#4985B7', '#C9BCF2', '#3DB88E', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red')) +
#'     # scale_color_manual(values = c('black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'red', 'red', 'red')) +
#'     # scale_fill_manual(values = c('black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'red', 'red', 'red')) +
#'     #'#212738', '#212738', '#212738', '#212738', '#212738', '#5C8492', '#5C8492', '#498587', '#5C8492', '#498587', '#212738', '#5C8492', '#5C8492', '#498587', '#498587', '#5C8492', '#498587'
#'     scale_shape(guide = 'none') +
#'     coord_cartesian(clip = "off") +
#' 
#'     labs(x = 'Concentration (µM)', y = "Normalized length",
#'          color = 'Route', fill = 'Route', shape = 'Replicate') +
#'     facet_grid(rows = vars(genotype), cols = vars(treatment), scales = 'free_x',
#'                labeller = as_labeller(c(AlbendazoleSulfoxide = 'Albendazole sulfoxide',
#'                                         Ivermectin = 'Ivermectin',
#'                                         Levamisole = 'Levamisole',
#'                                         `*agmo-1(e3016)*` = '*agmo-1(e3016)*',
#'                                         `*bus-5(br19)*` = '*bus-5(br19)*',
#'                                         `*che-1(p672)*` = '*che-1(p672)*',
#'                                         `*dyf-2(m160)*` = '*dyf-2(m160)*',
#'                                         `*eat-2(ad453)*` = '*eat-2(ad453)*',
#'                                         `*nhr-8(ok186)*` = '*nhr-8(ok186)*'))) +
#'     theme_nw2() +
#'     theme(
#'       legend.position = 'none',
#'       axis.text.x = element_markdown(angle = 0, hjust = 0.5),
#'       axis.title.x = element_markdown(face = 'plain'),
#'       axis.title.y = element_markdown(face = 'plain')) +
#'     NULL)
#' 
#' 
#' 
#' comps_norm_min <- curves %>% 
#'   mutate(
#'     comp = map(drc_norm_DMSO, ~EDcomp(.x, percVec = c(50, 50))),
#'     comp_df = map(comp, ~as_tibble(.x, rownames = 'Comparison'))
#'   ) %>% 
#'   select(treatment, comp_df) %>% 
#'   unnest(cols = comp_df) %>% 
#'   filter(str_detect(Comparison, 'N2')) %>% 
#'   mutate(strains = str_remove_all(Comparison, str_c('N2', ':', '/', '50/50', sep = '|'))) %>% 
#'   select(strains, treatment, p.value =  `p-value`) 
#' 
#' AE501_AZS_pv_norm_min <- subset(comps_norm_min, strains == "AE501" & treatment == 'AlbendazoleSulfoxide')
#' DA453_AZS_pv_norm_min <- subset(comps_norm_min, strains == "DA453" & treatment == 'AlbendazoleSulfoxide')
#' DC19_AZS_pv_norm_min <- subset(comps_norm_min, strains == "DC19" & treatment == 'AlbendazoleSulfoxide')
#' LC144_AZS_pv_norm_min <- subset(comps_norm_min, strains == "LC144" & treatment == 'AlbendazoleSulfoxide')
#' PR672_AZS_pv_norm_min <- subset(comps_norm_min, strains == "PR672" & treatment == 'AlbendazoleSulfoxide')
#' SP1234_AZS_pv_norm_min <- subset(comps_norm_min, strains == "SP1234" & treatment == 'AlbendazoleSulfoxide')
#' 
#' AE501_IVM_pv_norm_min <- subset(comps_norm_min, strains == "AE501" & treatment == 'Ivermectin')
#' DA453_IVM_pv_norm_min <- subset(comps_norm_min, strains == "DA453" & treatment == 'Ivermectin')
#' DC19_IVM_pv_norm_min <- subset(comps_norm_min, strains == "DC19" & treatment == 'Ivermectin')
#' LC144_IVM_pv_norm_min <- subset(comps_norm_min, strains == "LC144" & treatment == 'Ivermectin')
#' PR672_IVM_pv_norm_min <- subset(comps_norm_min, strains == "PR672" & treatment == 'Ivermectin')
#' SP1234_IVM_pv_norm_min <- subset(comps_norm_min, strains == "SP1234" & treatment == 'Ivermectin')
#' 
#' AE501_LEV_pv_norm_min <- subset(comps_norm_min, strains == "AE501" & treatment == 'Levamisole')
#' DA453_LEV_pv_norm_min <- subset(comps_norm_min, strains == "DA453" & treatment == 'Levamisole')
#' DC19_LEV_pv_norm_min <- subset(comps_norm_min, strains == "DC19" & treatment == 'Levamisole')
#' LC144_LEV_pv_norm_min <- subset(comps_norm_min, strains == "LC144" & treatment == 'Levamisole')
#' PR672_LEV_pv_norm_min <- subset(comps_norm_min, strains == "PR672" & treatment == 'Levamisole')
#' SP1234_LEV_pv_norm_min <- subset(comps_norm_min, strains == "SP1234" & treatment == 'Levamisole')
#' 
#' ec50_norm_min_for_table <- curves %>%
#'   unnest(tidy_norm_min) %>%
#'   select(-drc_norm_min, -data, -glance_norm_min, -drc_norm_DMSO, -glance_norm_DMSO, -tidy_norm_DMSO, -drc_sqrt_norm_DMSO,
#'          -glance_sqrt_norm_DMSO, -tidy_sqrt_norm_DMSO, -drc_sqrt_norm_min, -glance_sqrt_norm_min, -tidy_sqrt_norm_min) %>%
#'   filter(term == 'e') %>%
#'   select(-term, -params) %>%
#'   mutate(xmax = estimate + (0.434 * std.error / estimate),
#'          xmin = estimate - (0.434 * std.error / estimate)) %>%
#'   mutate(
#'     xmin = case_when(
#'       is.na(xmin) & treatment == 'Ivermectin' ~ 0.00025,
#'       is.na(xmin) & treatment == 'AlbendazoleSulfoxide' ~ 0.15,
#'       TRUE ~ xmin),
#'     xmax = case_when(
#'       is.na(xmax) & treatment == 'Ivermectin' ~ 1,
#'       is.na(xmax) & treatment == 'AlbendazoleSulfoxide' ~ 450,
#'       TRUE ~ xmax
#'     )
#'   ) %>%
#'   select( -group) %>%
#'   rename(estimate_norm_min = estimate, std.error_norm_min = std.error,
#'          statistic_norm_min = statistic, p.value_norm_min = p.value, xmax_norm_min = xmax, xmin_norm_min = xmin)
#' 
#' 
#' tab_norm_min  <- ec50_norm_min_for_table %>% 
#'   rename(strains = curve) %>%
#'   select(-xmax_norm_min , -xmin_norm_min , -p.value_norm_min) %>% 
#'   mutate(p.value = case_when(
#'     strains == 'AE501' & treatment == 'AlbendazoleSulfoxide' ~ AE501_AZS_pv_norm_min$p.value,
#'     strains == 'DA453' & treatment == 'AlbendazoleSulfoxide' ~ DA453_AZS_pv_norm_min$p.value,
#'     strains == 'DC19' & treatment == 'AlbendazoleSulfoxide' ~ DC19_AZS_pv_norm_min$p.value,
#'     strains == 'LC144' & treatment == 'AlbendazoleSulfoxide' ~ LC144_AZS_pv_norm_min$p.value,
#'     strains == 'PR672' & treatment == 'AlbendazoleSulfoxide' ~ PR672_AZS_pv_norm_min$p.value,
#'     strains == 'SP1234' & treatment == 'AlbendazoleSulfoxide' ~ SP1234_AZS_pv_norm_min$p.value,
#'     
#'     strains == 'AE501' & treatment == 'Ivermectin' ~ AE501_IVM_pv_norm_min$p.value,
#'     strains == 'DA453' & treatment == 'Ivermectin' ~ DA453_IVM_pv_norm_min$p.value,
#'     strains == 'DC19' & treatment == 'Ivermectin' ~ DC19_IVM_pv_norm_min$p.value,
#'     strains == 'LC144' & treatment == 'Ivermectin' ~ LC144_IVM_pv_norm_min$p.value,
#'     strains == 'PR672' & treatment == 'Ivermectin' ~ PR672_IVM_pv_norm_min$p.value,
#'     strains == 'SP1234' & treatment == 'Ivermectin' ~ SP1234_IVM_pv_norm_min$p.value,
#'     
#'     strains == 'AE501' & treatment == 'Levamisole' ~ AE501_LEV_pv_norm_min$p.value,
#'     strains == 'DA453' & treatment == 'Levamisole' ~ DA453_LEV_pv_norm_min$p.value,
#'     strains == 'DC19' & treatment == 'Levamisole' ~ DC19_LEV_pv_norm_min$p.value,
#'     strains == 'LC144' & treatment == 'Levamisole' ~ LC144_LEV_pv_norm_min$p.value,
#'     strains == 'PR672' & treatment == 'Levamisole' ~ PR672_LEV_pv_norm_min$p.value,
#'     strains == 'SP1234' & treatment == 'Levamisole' ~ SP1234_LEV_pv_norm_min$p.value,
#'   )) %>%
#'   arrange(-estimate_norm_min)  %>% 
#'   arrange(fct_relevel(strains, 'N2')) %>% 
#'   select(-statistic_norm_min) %>% 
#'   mutate(p.value = case_when(
#'     p.value > 0.05 ~ as.character(round(p.value, 3)),
#'     p.value < 0.0001 ~ paste(as.character(round(p.value, 3)), '****', sep=""), 
#'     p.value < 0.001 ~ paste(as.character(round(p.value, 3)), '***', sep=""), 
#'     p.value < 0.01 ~ paste(as.character(round(p.value, 3)), '**', sep=""), 
#'     p.value < 0.05 ~ paste(as.character(round(p.value, 3)), '*', sep=""), 
#'     TRUE ~ as.character(p.value)
#'   )) %>%
#'   mutate(genotype = case_when(
#'     strains == 'N2' ~ 'Wild type',
#'     strains == 'AE501' ~ '*nhr-8(ok186)*',
#'     strains == 'DA453' ~ '*eat-2(ad453)*',
#'     strains == 'DC19' ~ '*bus-5(br19)*',
#'     strains == 'LC144' ~ '*agmo-1(e3016)*',
#'     strains == 'PR672' ~ '*che-1(p672)*',
#'     strains == 'SP1234' ~ '*dyf-2(m160)*'
#'   ))
#' 
#' ec50_table_norm_min  <- tab_norm_min  %>% 
#'   mutate(
#'     treatment = case_when(
#'       treatment == 'AlbendazoleSulfoxide' ~ 'Albendazole sulfoxide',
#'       TRUE ~ treatment),
#'     std_error = case_when(
#'       treatment == 'Ivermectin' ~ as.character(round(std.error_norm_min, 4)),
#'       TRUE ~ as.character(round(std.error_norm_min, 2))),
#'     # p.value = round(p.value, 4),
#'     estimate = case_when(
#'       treatment == 'Ivermectin' ~ as.character(round(estimate_norm_min, 4)),
#'       TRUE ~ as.character(round(estimate_norm_min, 2))
#'     )) %>% 
#'   select(-strains, -estimate_norm_min, -std.error_norm_min) %>% 
#'   gt(groupname_col = "treatment",
#'      # rowname_col = "genotype"
#'   ) %>%
#'   tab_options(
#'     table.align = "left",
#'     table.font.size = 10) %>%
#'   tab_header(title = "Dose-Response Estimates") %>%
#'   # set the font
#'   tab_style(
#'     style = cell_text(font = google_font("Helvetica")),
#'     locations = list(
#'       cells_column_labels(),
#'       cells_stub(),
#'       cells_title(groups = "title"),
#'       cells_body(columns = everything())
#'     )) %>%
#'   # rownames
#'   tab_style(
#'     style = cell_text(style = 'italic'),
#'     locations = cells_row_groups()
#'   ) %>%
#'   # columns
#'   cols_label(
#'     treatment = md("**Treatment**"),
#'     genotype = md("**Genotype**"),
#'     # strains = md("**Strain name**"),
#'     estimate = md("**EC<sub>50</sub>**"),
#'     std_error = md("**SE**"), 
#'     p.value = md("**_p_<br>(N2 == strain)**")) %>%
#'   # numeric columns
#'   # fmt_symbol_first(column = estimate, symbol = " µM", last_row_n = 8, gfont = "Helvetica", decimals = 2) %>%
#'   fmt_markdown(column = genotype, rows = everything()) %>% 
#'   # fmt_scientific(columns = p.value, rows = p.value < 0.001) %>% 
#'   cols_align(align = 'left', columns = everything()) %>%
#'   # color sig p-value
#'   # tab_style(
#'   #   style = cell_fill(color = 'grey', alpha = 0.5),
#'   #   locations = cells_body(
#'   #     rows = p.value < 0.05,
#'   #     columns = p.value)
#'   #   ) %>% 
#'   # color strains
#'   tab_style(
#'     style = cell_text(color = '#6D9DC5'),
#'     locations = cells_body(
#'       rows = genotype == '*agmo-1(e3016)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#6D9DC5'),
#'     locations = cells_body(
#'       rows = genotype == '*bus-5(br19)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#212738'),
#'     locations = cells_body(
#'       rows = genotype == '*che-1(p672)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#212738'),
#'     locations = cells_body(
#'       rows = genotype == '*dyf-2(m160)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#2A7F62'),
#'     locations = cells_body(
#'       rows = genotype == '*eat-2(ad453)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#2A7F62'),
#'     locations = cells_body(
#'       rows = genotype == '*nhr-8(ok186)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#A5243D'),
#'     locations = cells_body(
#'       rows = genotype == 'N2',
#'       columns = genotype)
#'   )
#' 
#' (full_norm_min <- plot_grid(curve_plot_norm_min + remove_legend(), bstfun::as_ggplot(ec50_table_norm_min),
#'                              nrow = 1, rel_widths = c(1.5, 1),
#'                              labels = 'auto'))
#' 
#' 
#' 
#' (curve_plot_sqrt_norm_min <- revised_data %>%
#'     filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin')) %>%
#'     mutate(route = case_when(
#'       strains == 'N2' ~ 'Wild type',
#'       strains %in% c('PR672', 'SP1234') ~ 'Amphid',
#'       strains %in% c('DC19', 'LC144') ~ 'Cuticle',
#'       strains %in% c('DA453', 'AE501') ~ 'Digestive'
#'     )) %>%
#'     # don't show N2 here
#'     filter(strains != 'N2') %>%
#'     ggplot() +
#'     geom_vline(data = ec50_pivoted_1 %>% filter(strains != 'N2', norm_method == 'estimate_sqrt_norm_min'),
#'                aes(xintercept = estimate_value, group = genotype), linetype = 'dashed', color = '#212738', size = 0.3) +
#'     geom_vline(data = ec50_pivoted_2 %>% filter(strains != 'N2', norm_method == 'estimate_sqrt_norm_min'),
#'                aes(xintercept = estimate_value, group = genotype), linetype = 'dashed', color = '#5C8492', size = 0.3) +
#'     geom_vline(data = ec50_pivoted_3 %>% filter(strains != 'N2', norm_method == 'estimate_sqrt_norm_min'),
#'                aes(xintercept = estimate_value, group = genotype), linetype = 'dashed', color = '#498587', size = 0.3) +
#'     geom_vline(data = ec50_pivoted_4 %>% filter(strains != 'N2', norm_method == 'estimate_sqrt_norm_min'),
#'                aes(xintercept = estimate_value, group = genotype), linetype = 'dashed', color = '#6D9DC5', size = 0.3) +
#'     geom_vline(data = ec50_average %>% filter(strains != 'N2', norm_method == 'estimate_sqrt_norm_min'),
#'                aes(xintercept = estimate_avg, group = genotype), linetype = 'dashed', color = 'black', size = 0.6) +
#'     geom_line(data = predict_1 %>% filter(strains != 'N2'),
#'               aes(x = conc, group = genotype, y = pred_norm_min), #, color = group_route
#'               size = 0.3, alpha =0.5, color = '#212738') +
#'     geom_line(data = predict_2 %>% filter(strains != 'N2'),
#'               aes(x = conc, group = genotype, y = pred_sqrt_norm_min), #, color = group_route
#'               size = 0.3, alpha =0.5, color = '#5C8492') +
#'     geom_line(data = predict_3 %>% filter(strains != 'N2', strains != 'DA453'),
#'               aes(x = conc, group = genotype, y = pred_sqrt_norm_min), #, color = group_route
#'               size = 0.3, alpha =0.5, color = '#498587') +
#'     geom_line(data = predict_4 %>% filter(strains != 'N2'),
#'               aes(x = conc, group = genotype, y = pred_sqrt_norm_min), #, color = group_route
#'               size = 0.3, alpha =0.5, color = '#6D9DC5') +
#'     geom_line(data = predict_all %>% filter(genotype != 'N2'),
#'               aes(x = conc, group = genotype, y = fit_line_sqrt_min), #, color = group_route
#'               size = 0.75, color = 'black', alpha = 1) +
#'     geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
#'                                                     strains != 'N2',
#'                                                     plate != '20220210-p22-EJG_1187',
#'                                                     # rep <= 3
#'     ),
#'     aes(x = conc, y = mean_sqrt_norm_min, color = plate),
#'     size = 0.5, width = 0.1, alpha = 0.6) +
#'     geom_hline(yintercept = 0.2, color = 'black', size = 0.25) +
#'     scale_x_log10(
#'       breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
#'       labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#'     annotation_logticks(sides = 'b', size = 0.25, outside = TRUE,
#'                         short = unit(0.1, "cm"),
#'                         mid = unit(0.15, "cm"),
#'                         long = unit(0.2, "cm")) +
#'     scale_y_continuous(limits = c(-0.2, 1.45), expand = c(0, 0), breaks = seq(-0.5, 1.2, 0.2)) +
#'     scale_color_manual(values = c('#212738', '#212738', '#212738', '#212738', '#212738', '#5C8492', '#5C8492', '#498587', '#5C8492','#498587', '#212738', '#5C8492', '#5C8492', '#498587', '#498587', '#5C8492', '#498587', 'blue')) +
#'     scale_fill_manual(values = c('#212738', '#212738', '#212738', '#212738', '#212738', '#5C8492', '#5C8492', '#498587', '#5C8492', '#498587', '#212738', '#5C8492', '#5C8492', '#498587', '#498587', '#5C8492', '#498587', 'blue')) +
#'     #scale_color_manual(values = c('#212738', '#622A7F', '#174636', '#6D9DC5', '#6F636B', '#2A7F62','#5C8492', '#82AC92', '#4985B7', '#C9BCF2', '#3DB88E', 'purple', 'orange', '#212738', '#622A7F', '#174636', '#6D9DC5', '#6F636B', '#2A7F62','#5C8492', '#82AC92', '#4985B7', '#C9BCF2', '#3DB88E', 'red', 'red', 'red', 'red', 'red', 'red', 'red','red', 'red')) +
#'     #scale_fill_manual(values = c('#212738', '#622A7F', '#174636', '#6D9DC5', '#6F636B', '#2A7F62','#5C8492', '#82AC92', '#4985B7', '#C9BCF2', '#3DB88E', 'purple', 'orange', '#212738', '#622A7F', '#174636', '#6D9DC5', '#6F636B', '#2A7F62','#5C8492', '#82AC92', '#4985B7', '#C9BCF2', '#3DB88E', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red')) +
#'     # scale_color_manual(values = c('black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'red', 'red', 'red')) +
#'     # scale_fill_manual(values = c('black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'red', 'red', 'red')) +
#'     #'#212738', '#212738', '#212738', '#212738', '#212738', '#5C8492', '#5C8492', '#498587', '#5C8492', '#498587', '#212738', '#5C8492', '#5C8492', '#498587', '#498587', '#5C8492', '#498587'
#'     scale_shape(guide = 'none') +
#'     coord_cartesian(clip = "off") +
#' 
#'     labs(x = 'Concentration (µM)', y = "",
#'          color = 'Route', fill = 'Route', shape = 'Replicate') +
#'     facet_grid(rows = vars(genotype), cols = vars(treatment), scales = 'free_x',
#'                labeller = as_labeller(c(AlbendazoleSulfoxide = 'Albendazole sulfoxide',
#'                                         Ivermectin = 'Ivermectin',
#'                                         Levamisole = 'Levamisole',
#'                                         `*agmo-1(e3016)*` = '*agmo-1(e3016)*',
#'                                         `*bus-5(br19)*` = '*bus-5(br19)*',
#'                                         `*che-1(p672)*` = '*che-1(p672)*',
#'                                         `*dyf-2(m160)*` = '*dyf-2(m160)*',
#'                                         `*eat-2(ad453)*` = '*eat-2(ad453)*',
#'                                         `*nhr-8(ok186)*` = '*nhr-8(ok186)*'))) +
#'     theme_nw2() +
#'     theme(
#'       legend.position = 'none',
#'       axis.text.x = element_markdown(angle = 0, hjust = 0.5),
#'       axis.title.x = element_markdown(face = 'plain'),
#'       axis.title.y = element_markdown(face = 'plain')) +
#'     NULL)
#' 
#' comps_sqrt_norm_min <- curves %>% 
#'   mutate(
#'     comp = map(drc_sqrt_norm_DMSO, ~EDcomp(.x, percVec = c(50, 50))),
#'     comp_df = map(comp, ~as_tibble(.x, rownames = 'Comparison'))
#'   ) %>% 
#'   select(treatment, comp_df) %>% 
#'   unnest(cols = comp_df) %>% 
#'   filter(str_detect(Comparison, 'N2')) %>% 
#'   mutate(strains = str_remove_all(Comparison, str_c('N2', ':', '/', '50/50', sep = '|'))) %>% 
#'   select(strains, treatment, p.value =  `p-value`) 
#' 
#' AE501_AZS_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "AE501" & treatment == 'AlbendazoleSulfoxide')
#' DA453_AZS_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "DA453" & treatment == 'AlbendazoleSulfoxide')
#' DC19_AZS_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "DC19" & treatment == 'AlbendazoleSulfoxide')
#' LC144_AZS_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "LC144" & treatment == 'AlbendazoleSulfoxide')
#' PR672_AZS_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "PR672" & treatment == 'AlbendazoleSulfoxide')
#' SP1234_AZS_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "SP1234" & treatment == 'AlbendazoleSulfoxide')
#' 
#' AE501_IVM_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "AE501" & treatment == 'Ivermectin')
#' DA453_IVM_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "DA453" & treatment == 'Ivermectin')
#' DC19_IVM_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "DC19" & treatment == 'Ivermectin')
#' LC144_IVM_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "LC144" & treatment == 'Ivermectin')
#' PR672_IVM_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "PR672" & treatment == 'Ivermectin')
#' SP1234_IVM_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "SP1234" & treatment == 'Ivermectin')
#' 
#' AE501_LEV_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "AE501" & treatment == 'Levamisole')
#' DA453_LEV_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "DA453" & treatment == 'Levamisole')
#' DC19_LEV_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "DC19" & treatment == 'Levamisole')
#' LC144_LEV_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "LC144" & treatment == 'Levamisole')
#' PR672_LEV_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "PR672" & treatment == 'Levamisole')
#' SP1234_LEV_pv_sqrt_norm_min <- subset(comps_sqrt_norm_min, strains == "SP1234" & treatment == 'Levamisole')
#' 
#' ec50_sqrt_norm_min_for_table <- curves %>%
#'   unnest(tidy_sqrt_norm_min) %>%
#'   select(-drc_sqrt_norm_min, -data, -glance_sqrt_norm_min, -drc_norm_min, -tidy_norm_min, -glance_norm_min, -drc_norm_DMSO, -glance_norm_DMSO, -tidy_norm_DMSO, -drc_sqrt_norm_DMSO,
#'          -glance_sqrt_norm_DMSO, -tidy_sqrt_norm_DMSO) %>%
#'   filter(term == 'e') %>%
#'   select(-term, -params) %>%
#'   mutate(xmax = estimate + (0.434 * std.error / estimate),
#'          xmin = estimate - (0.434 * std.error / estimate)) %>%
#'   mutate(
#'     xmin = case_when(
#'       is.na(xmin) & treatment == 'Ivermectin' ~ 0.00025,
#'       is.na(xmin) & treatment == 'AlbendazoleSulfoxide' ~ 0.15,
#'       TRUE ~ xmin),
#'     xmax = case_when(
#'       is.na(xmax) & treatment == 'Ivermectin' ~ 1,
#'       is.na(xmax) & treatment == 'AlbendazoleSulfoxide' ~ 450,
#'       TRUE ~ xmax
#'     )
#'   ) %>%
#'   select(-group) %>%
#'   rename(estimate_sqrt_norm_min = estimate, std.error_sqrt_norm_min = std.error,
#'          statistic_sqrt_norm_min = statistic, p.value_sqrt_norm_min = p.value, xmax_sqrt_norm_min = xmax, xmin_sqrt_norm_min = xmin)
#' 
#' 
#' tab_sqrt_norm_min  <- ec50_sqrt_norm_min_for_table %>% 
#'   rename(strains = curve) %>%
#'   select(-xmax_sqrt_norm_min , -xmin_sqrt_norm_min , -p.value_sqrt_norm_min) %>% 
#'   mutate(p.value = case_when(
#'     strains == 'AE501' & treatment == 'AlbendazoleSulfoxide' ~ AE501_AZS_pv_sqrt_norm_min$p.value,
#'     strains == 'DA453' & treatment == 'AlbendazoleSulfoxide' ~ DA453_AZS_pv_sqrt_norm_min$p.value,
#'     strains == 'DC19' & treatment == 'AlbendazoleSulfoxide' ~ DC19_AZS_pv_sqrt_norm_min$p.value,
#'     strains == 'LC144' & treatment == 'AlbendazoleSulfoxide' ~ LC144_AZS_pv_sqrt_norm_min$p.value,
#'     strains == 'PR672' & treatment == 'AlbendazoleSulfoxide' ~ PR672_AZS_pv_sqrt_norm_min$p.value,
#'     strains == 'SP1234' & treatment == 'AlbendazoleSulfoxide' ~ SP1234_AZS_pv_sqrt_norm_min$p.value,
#'     
#'     strains == 'AE501' & treatment == 'Ivermectin' ~ AE501_IVM_pv_sqrt_norm_min$p.value,
#'     strains == 'DA453' & treatment == 'Ivermectin' ~ DA453_IVM_pv_sqrt_norm_min$p.value,
#'     strains == 'DC19' & treatment == 'Ivermectin' ~ DC19_IVM_pv_sqrt_norm_min$p.value,
#'     strains == 'LC144' & treatment == 'Ivermectin' ~ LC144_IVM_pv_sqrt_norm_min$p.value,
#'     strains == 'PR672' & treatment == 'Ivermectin' ~ PR672_IVM_pv_sqrt_norm_min$p.value,
#'     strains == 'SP1234' & treatment == 'Ivermectin' ~ SP1234_IVM_pv_sqrt_norm_min$p.value,
#'     
#'     strains == 'AE501' & treatment == 'Levamisole' ~ AE501_LEV_pv_sqrt_norm_min$p.value,
#'     strains == 'DA453' & treatment == 'Levamisole' ~ DA453_LEV_pv_sqrt_norm_min$p.value,
#'     strains == 'DC19' & treatment == 'Levamisole' ~ DC19_LEV_pv_sqrt_norm_min$p.value,
#'     strains == 'LC144' & treatment == 'Levamisole' ~ LC144_LEV_pv_sqrt_norm_min$p.value,
#'     strains == 'PR672' & treatment == 'Levamisole' ~ PR672_LEV_pv_sqrt_norm_min$p.value,
#'     strains == 'SP1234' & treatment == 'Levamisole' ~ SP1234_LEV_pv_sqrt_norm_min$p.value,
#'   )) %>%
#'   arrange(-estimate_sqrt_norm_min)  %>% 
#'   arrange(fct_relevel(strains, 'N2')) %>% 
#'   select(-statistic_sqrt_norm_min) %>% 
#'   mutate(p.value = case_when(
#'     p.value > 0.05 ~ as.character(round(p.value, 3)),
#'     p.value < 0.0001 ~ paste(as.character(round(p.value, 3)), '****', sep=""), 
#'     p.value < 0.001 ~ paste(as.character(round(p.value, 3)), '***', sep=""), 
#'     p.value < 0.01 ~ paste(as.character(round(p.value, 3)), '**', sep=""), 
#'     p.value < 0.05 ~ paste(as.character(round(p.value, 3)), '*', sep=""), 
#'     TRUE ~ as.character(p.value)
#'   )) %>%
#'   # mutate(p.value = case_when(
#'   #   p.value > 0.05 ~ as.character(round(p.value, 3)),
#'   #   p.value < 0.0001 ~ '****',
#'   #   p.value < 0.001 ~ '***',
#'   #   p.value < 0.01 ~ '**',
#'   #   p.value < 0.05 ~ '*',
#'   #   TRUE ~ as.character(p.value)
#'   # )) %>%
#'   mutate(genotype = case_when(
#'     strains == 'N2' ~ 'Wild type',
#'     strains == 'AE501' ~ '*nhr-8(ok186)*',
#'     strains == 'DA453' ~ '*eat-2(ad453)*',
#'     strains == 'DC19' ~ '*bus-5(br19)*',
#'     strains == 'LC144' ~ '*agmo-1(e3016)*',
#'     strains == 'PR672' ~ '*che-1(p672)*',
#'     strains == 'SP1234' ~ '*dyf-2(m160)*'
#'   ))
#' 
#' (ec50_table_sqrt_norm_min  <- tab_sqrt_norm_min  %>% 
#'   mutate(
#'     treatment = case_when(
#'       treatment == 'AlbendazoleSulfoxide' ~ 'Albendazole sulfoxide',
#'       TRUE ~ treatment),
#'     std_error = case_when(
#'       treatment == 'Ivermectin' ~ as.character(round(std.error_sqrt_norm_min, 4)),
#'       TRUE ~ as.character(round(std.error_sqrt_norm_min, 2))),
#'     # p.value = round(p.value, 4),
#'     estimate = case_when(
#'       treatment == 'Ivermectin' ~ as.character(round(estimate_sqrt_norm_min, 4)),
#'       TRUE ~ as.character(round(estimate_sqrt_norm_min, 2))
#'     )) %>% 
#'   select(-strains, -estimate_sqrt_norm_min, -std.error_sqrt_norm_min) %>% 
#'   gt(groupname_col = "treatment",
#'      # rowname_col = "genotype"
#'   ) %>%
#'   tab_options(
#'     table.align = "left",
#'     table.font.size = 10) %>%
#'   tab_header(title = "Dose-Response Estimates") %>%
#'   # set the font
#'   tab_style(
#'     style = cell_text(font = google_font("Helvetica")),
#'     locations = list(
#'       cells_column_labels(),
#'       cells_stub(),
#'       cells_title(groups = "title"),
#'       cells_body(columns = everything())
#'     )) %>%
#'   # rownames
#'   tab_style(
#'     style = cell_text(style = 'italic'),
#'     locations = cells_row_groups()
#'   ) %>%
#'   # columns
#'   cols_label(
#'     treatment = md("**Treatment**"),
#'     genotype = md("**Genotype**"),
#'     # strains = md("**Strain name**"),
#'     estimate = md("**EC<sub>50</sub>**"),
#'     std_error = md("**SE**"), 
#'     p.value = md("**_p_<br>(N2 == strain)**")) %>%
#'   # numeric columns
#'   # fmt_symbol_first(column = estimate, symbol = " µM", last_row_n = 8, gfont = "Helvetica", decimals = 2) %>%
#'   fmt_markdown(column = genotype, rows = everything()) %>% 
#'   # fmt_scientific(columns = p.value, rows = p.value < 0.001) %>% 
#'   cols_align(align = 'left', columns = everything()) %>%
#'   # color sig p-value
#'   # tab_style(
#'   #   style = cell_fill(color = 'grey', alpha = 0.5),
#'   #   locations = cells_body(
#'   #     rows = p.value < 0.05,
#'   #     columns = p.value)
#'   #   ) %>% 
#'   # color strains
#'   tab_style(
#'     style = cell_text(color = '#6D9DC5'),
#'     locations = cells_body(
#'       rows = genotype == '*agmo-1(e3016)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#6D9DC5'),
#'     locations = cells_body(
#'       rows = genotype == '*bus-5(br19)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#212738'),
#'     locations = cells_body(
#'       rows = genotype == '*che-1(p672)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#212738'),
#'     locations = cells_body(
#'       rows = genotype == '*dyf-2(m160)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#2A7F62'),
#'     locations = cells_body(
#'       rows = genotype == '*eat-2(ad453)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#2A7F62'),
#'     locations = cells_body(
#'       rows = genotype == '*nhr-8(ok186)*',
#'       columns = genotype)
#'   ) %>%
#'   tab_style(
#'     style = cell_text(color = '#A5243D'),
#'     locations = cells_body(
#'       rows = genotype == 'N2',
#'       columns = genotype)
#'   ))
#' 
#' joined_test <- tab_sqrt_norm_min %>%
#'   select()
#' 
#' (full_sqrt_norm_min <- plot_grid(curve_plot_sqrt_norm_min + remove_legend(), bstfun::as_ggplot(ec50_table_sqrt_norm_min),
#'                              nrow = 1, rel_widths = c(1.5, 1),
#'                              labels = 'auto'))
#' 
#' ggsave(here('/Users/elenagr/Desktop/full_sqrt_norm_min.pdf'), full_sqrt_norm_min, width = 400, height = 450, units = 'mm')
#' 
#' (full_norm_min <- plot_grid(curve_plot_norm_min + remove_legend(), bstfun::as_ggplot(ec50_table_norm_min),
#'                                  nrow = 1, rel_widths = c(1.5, 1),
#'                                  labels = 'auto'))
#' 
#' ggsave(here('/Users/elenagr/Desktop/full_norm_min.pdf'), full_norm_min, width = 400, height = 450, units = 'mm')
#' 
#' (full_sqrt_norm_DMSO <- plot_grid(curve_plot_sqrt_norm_DMSO + remove_legend(), bstfun::as_ggplot(ec50_table_sqrt_norm_DMSO),
#'                                  nrow = 1, rel_widths = c(1.5, 1),
#'                                  labels = 'auto'))
#' 
#' ggsave(here('/Users/elenagr/Desktop/full_sqrt_norm_DMSO.pdf'), full_sqrt_norm_DMSO, width = 400, height = 450, units = 'mm')
#' 
#' (full_norm_DMSO <- plot_grid(curve_plot_norm_DMSO + remove_legend(), bstfun::as_ggplot(ec50_table_norm_DMSO),
#'                                  nrow = 1, rel_widths = c(1.5, 1),
#'                                  labels = 'auto'))
#' 
#' ggsave(here('/Users/elenagr/Desktop/full_norm_DMSO.pdf'), full_norm_DMSO, width = 400, height = 450, units = 'mm')
#' 
#' (Fig_2_comp <- ggarrange(full_norm_DMSO, full_sqrt_norm_DMSO, full_norm_min, full_sqrt_norm_min))
#' 
#' 
#' ggsave(here('/Users/elenagr/Desktop/Fig_2_comp.pdf'), Fig_2_comp, width = 400, height = 450, units = 'mm')
#' ggsave(here('/Users/elenagr/Desktop/Fig_2_comp.png'), Fig_2_comp, width = 400, height = 450, units = 'mm')
#' 
#' 
