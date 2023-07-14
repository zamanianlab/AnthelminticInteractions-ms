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


# tables
library(gt)
library(gtExtras)

# stats
library(drc)
library(broom)

# misc
library(conflicted)
library(here)
library(harrypotter)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

dr_data <- read_rds(here('data_EGR/dr_data_EGR.rds')) %>% # when box : 'data_EGR/dr_data_EGR.rds'
  filter(case_when(
    metadata_date == '20211115' & metadata_plate %in% c('p01', 'p02', 'p04', 'p05') & treatment %in% c('DMSO', 'Untreated', 'AlbendazoleSulfoxide') ~ TRUE,
    metadata_date == '20211118' & metadata_plate %in% c('p01', 'p02', 'p03', 'p04') & treatment %in% c('DMSO', 'Untreated', 'AlbendazoleSulfoxide', 'Levamisole') ~ TRUE,
    metadata_date == '20211208' & metadata_plate %in% c('p02') & treatment %in% c('DMSO', 'Untreated', 'AlbendazoleSulfoxide', 'Levamisole') ~ TRUE,
    metadata_date == '20211208' & metadata_plate %in% c('p01') & treatment %in% c('DMSO', 'Untreated', 'AlbendazoleSulfoxide') ~ TRUE,
    metadata_date == '20211216' & metadata_plate %in% c('p14', 'p27', 'p28') ~ TRUE,
    metadata_date == '20220131' & metadata_plate %in% c('p01', 'p02', 'p03', 'p04') ~ TRUE,
    metadata_date == '20220203' & metadata_plate %in% c('p01', 'p02', 'p03', 'p04') ~ TRUE,
    metadata_date == '20220210' & metadata_plate %in% c('p01', 'p08', 'p17', 'p22', 'p23') ~ TRUE,
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

# calculate all the averages for the highest concentrations - to be the average minimum length in normalizations
IVM_high_mean <- trimmed_data %>%
  filter(treatment == 'Ivermectin', conc == '0.5uM') %>%
  group_by(plate) %>% 
  summarise(IVM_high_mean = mean(raw_length, na.rm = FALSE),
            sqrt_IVM_high_mean = mean(sqrt_length, na.rm = FALSE))  %>%
  ungroup()

LEV_high_mean <- trimmed_data %>%
  filter(treatment == 'Levamisole', conc == '50uM') %>%
  group_by(plate) %>% 
  summarise(LEV_high_mean = mean(raw_length, na.rm = FALSE),
            sqrt_LEV_high_mean = mean(sqrt_length, na.rm = FALSE)) %>% 
  ungroup()

AZS_high_mean <- trimmed_data %>%
  filter(treatment == 'AlbendazoleSulfoxide', conc == '300uM') %>%
  group_by(plate) %>% 
  summarise(AZS_high_mean = mean(raw_length, na.rm = FALSE),
            sqrt_AZS_high_mean = mean(sqrt_length, na.rm = FALSE)) %>%
  ungroup()

DMSO_mean <- trimmed_data %>%
  filter(treatment == 'DMSO') %>%
  group_by(plate) %>% 
  summarise(DMSO_mean = mean(raw_length, na.rm = FALSE),
            sqrt_DMSO_mean = mean(sqrt_length, na.rm = FALSE)) %>%
  ungroup()

revised_data <- trimmed_data %>% 
  left_join(DMSO_mean) %>%
  mutate(conc = case_when(
    treatment == 'DMSO' ~ 0.01,
    treatment == 'Untreated' ~ NA_real_,
    treatment != 'DMSO' ~ as.numeric(str_remove(conc, 'uM'))
  )) %>%
  mutate(norm_DMSO = raw_length / DMSO_mean) 

# to add three supplementary normalizations to the revised dataframe
revised_data <- revised_data %>% 
  left_join(LEV_high_mean) %>%
  left_join(IVM_high_mean) %>%
  left_join(AZS_high_mean) %>%
  left_join(DMSO_mean) %>%
  mutate(min_value = case_when(
    treatment == 'Ivermectin' ~ IVM_high_mean,
    treatment == 'Levamisole' ~ LEV_high_mean,
    treatment == 'AlbendazoleSulfoxide' ~ AZS_high_mean,
    treatment == 'DMSO' ~ DMSO_mean
  )) %>%
  mutate(sqrt_min_value = case_when(
    treatment == 'Ivermectin' ~ sqrt_IVM_high_mean,
    treatment == 'Levamisole' ~ sqrt_LEV_high_mean,
    treatment == 'AlbendazoleSulfoxide' ~ sqrt_AZS_high_mean,
    treatment == 'DMSO' ~ sqrt_DMSO_mean
  )) %>%
  mutate(sqrt_norm_DMSO = sqrt_length / sqrt_DMSO_mean) %>%
  mutate(x_minus_highconc = case_when(
    treatment == 'Ivermectin' ~ raw_length - IVM_high_mean,
    treatment == 'Levamisole' ~ raw_length - LEV_high_mean,
    treatment == 'AlbendazoleSulfoxide' ~ raw_length - AZS_high_mean,
    treatment == 'DMSO' ~ raw_length - DMSO_mean
  )) %>%
  mutate(sqrt_x_minus_highconc = case_when(
    treatment == 'Ivermectin' ~ sqrt_length - sqrt_IVM_high_mean,
    treatment == 'Levamisole' ~ sqrt_length - sqrt_LEV_high_mean,
    treatment == 'AlbendazoleSulfoxide' ~ sqrt_length - sqrt_AZS_high_mean,
    treatment == 'DMSO' ~ sqrt_length - sqrt_DMSO_mean
  )) %>%
  mutate(DMSO_minus_highconc = case_when(
    treatment == 'Ivermectin' ~ DMSO_mean - IVM_high_mean,
    treatment == 'Levamisole' ~ DMSO_mean - LEV_high_mean,
    treatment == 'AlbendazoleSulfoxide' ~ DMSO_mean - AZS_high_mean,
    treatment == 'DMSO' ~ DMSO_mean - DMSO_mean
  )) %>%
  mutate(sqrt_DMSO_minus_highconc = case_when(
    treatment == 'Ivermectin' ~ sqrt_DMSO_mean - sqrt_IVM_high_mean,
    treatment == 'Levamisole' ~ sqrt_DMSO_mean - sqrt_LEV_high_mean,
    treatment == 'AlbendazoleSulfoxide' ~ sqrt_DMSO_mean - sqrt_AZS_high_mean,
    treatment == 'DMSO' ~ sqrt_DMSO_mean - sqrt_DMSO_mean
  )) %>%
  mutate(norm_min = x_minus_highconc / DMSO_minus_highconc) %>%
  mutate(sqrt_norm_min = sqrt_x_minus_highconc / sqrt_DMSO_minus_highconc)

well_summary <- revised_data %>% 
  select(plate, well, conc, metadata_date, strains, treatment, genotype, raw_length, norm_DMSO, sqrt_length, min_value, sqrt_norm_DMSO, norm_min, sqrt_norm_min) %>% 
  group_by(metadata_date, plate, genotype, strains, well, treatment, conc) %>% 
  summarise(mean_raw_length = mean(raw_length),
            mean_sqrt_length = mean(sqrt_length),
            mean_norm_DMSO = mean(norm_DMSO),
            mean_sqrt_norm_DMSO = mean(sqrt_norm_DMSO),
            mean_norm_min = mean(norm_min),
            mean_sqrt_norm_min = mean(sqrt_norm_min)
  ) %>% 
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  ))


# plot all DMSO and Untreated for every strain and color by date -----------------
well_summary <- well_summary %>%
  mutate(drug_stock = case_when(
    metadata_date %in% c('20211115', '20211118', '20211208', '20211216') ~ 'before',
    metadata_date %in% c('20220131', '20220203', '20220210', '20220311', '20220324', '20220325', '20220331', '20220408') ~ 'after',
  )) %>%
  mutate(group = case_when(
    plate %in% c('20211115-p01-ECG_959', '20211216-p27-EJG_1081', '20220203-p01-EJG_1142') ~ '1',
    plate %in% c('20211118-p01-EJG_975', '20211216-p28-EJG_1082', '20220131-p03-EJG_1135') ~ '2',
    plate %in% c('20211208-p01-EJG_1022', '20220203-p02-EJG_1143', '20220210-p22-EJG_1187') ~ '3',
    plate %in% c('20211216-p14-EJG_1074', '20220210-p15-EJG_1180', '20220210-p08-EJG_1172') ~ '4',
    plate %in% c('20220131-p04-EJG_1136', '20220324-p05-EJG_1292', '20220325-p01-EJG_1299') ~ '5',
    plate %in% c('20220210-p01-EJG_1178', '20220131-p02-EJG_1134', '20211115-p05-ECG_963') ~ '6',
    plate %in% c('20220331-p08-EJG_1318', '20220331-p02-EJG_1312', '20220210-p17-EJG_1182', '20211208-p02-EJG_1023') ~ '7',
    plate %in% c('20220408-p01-MGC_1351', '20220210-p23-EJG_1188', '20211118-p04-EJG_978') ~ '8',
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

curves <- revised_data %>%
  left_join(reps) %>%
  # use only the most recent 3 replicates (not just those that used the same drug preparation)
  filter(!treatment %in% c('DMSO', 'Untreated')) %>%
  mutate(group = case_when(
    plate %in% c('20211115-p01-ECG_959', '20211216-p27-EJG_1081', '20220203-p01-EJG_1142') ~ '1',
    plate %in% c('20211118-p01-EJG_975', '20211216-p28-EJG_1082', '20220131-p03-EJG_1135') ~ '2',
    plate %in% c('20211208-p01-EJG_1022', '20220203-p02-EJG_1143', '20220210-p22-EJG_1187') ~ '3',
    plate %in% c('20211216-p14-EJG_1074', '20220210-p15-EJG_1180', '20220210-p08-EJG_1172') ~ '4',
    plate %in% c('20220131-p04-EJG_1136', '20220324-p05-EJG_1292', '20220325-p01-EJG_1299') ~ '5',
    plate %in% c('20220210-p01-EJG_1178', '20220131-p02-EJG_1134', '20211115-p05-ECG_963') ~ '6',
    plate %in% c('20220331-p08-EJG_1318', '20220331-p02-EJG_1312', '20220210-p17-EJG_1182', '20211208-p02-EJG_1023') ~ '7',
    plate %in% c('20220408-p01-MGC_1351', '20220210-p23-EJG_1188', '20211118-p04-EJG_978') ~ '8',
  )) %>%
  # uncomment below to fit to the summarized well data
  group_by(plate, well, treatment, strains, genotype, conc, group) %>%
  summarise(mean_raw_length = mean(raw_length),
            mean_sqrt_length = mean(sqrt_length),
            mean_norm_DMSO = mean(norm_DMSO),
            mean_sqrt_norm_DMSO = mean(sqrt_norm_DMSO),
            mean_norm_min = mean(norm_min),
            mean_sqrt_norm_min = mean(sqrt_norm_min)
  ) %>%
  ungroup() %>%
  group_nest(treatment, group) %>%
  # rowwise() %>%
  # fit the models using all the data but group by strain, give some parameters
  mutate(params = list(c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA))) %>%
  mutate(drc_norm_DMSO  = map2(data, params, ~ drm(.x$mean_norm_DMSO ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_norm_DMSO = map(drc_norm_DMSO, glance)) %>%
  mutate(tidy_norm_DMSO = map(drc_norm_DMSO, tidy)) %>%
  mutate(drc_sqrt_norm_DMSO = map2(data, params, ~ drm(.x$mean_sqrt_norm_DMSO ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_sqrt_norm_DMSO = map(drc_sqrt_norm_DMSO, glance)) %>%
  mutate(tidy_sqrt_norm_DMSO = map(drc_sqrt_norm_DMSO, tidy)) %>%
  mutate(drc_norm_min = map2(data, params, ~ drm(.x$mean_norm_min ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_norm_min = map(drc_norm_min, glance)) %>%
  mutate(tidy_norm_min = map(drc_norm_min, tidy)) %>%
  mutate(drc_sqrt_norm_min = map2(data, params, ~ drm(.x$mean_sqrt_norm_min ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_sqrt_norm_min = map(drc_sqrt_norm_min, glance)) %>%
  mutate(tidy_sqrt_norm_min = map(drc_sqrt_norm_min, tidy))

curves_all <- revised_data %>%
  left_join(reps) %>%
  # use only the most recent 3 replicates (not just those that used the same drug preparation)
  filter(!treatment %in% c('DMSO', 'Untreated')) %>%
  # uncomment below to fit to the summarized well data
  group_by(plate, treatment, strains, genotype, conc) %>%
  summarise(mean_raw_length = mean(raw_length),
            mean_sqrt_length = mean(sqrt_length),
            mean_norm_DMSO = mean(norm_DMSO),
            mean_sqrt_norm_DMSO = mean(sqrt_norm_DMSO),
            mean_norm_min = mean(norm_min),
            mean_sqrt_norm_min = mean(sqrt_norm_min)
  ) %>%
  ungroup() %>%
  group_nest(treatment) %>%
  # rowwise() %>%
  # fit the models using all the data but group by strain, give some parameters
  mutate(params = list(c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA))) %>%
  mutate(drc_norm_DMSO  = map2(data, params, ~ drm(.x$mean_norm_DMSO ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_norm_DMSO = map(drc_norm_DMSO, glance)) %>%
  mutate(tidy_norm_DMSO = map(drc_norm_DMSO, tidy)) %>%
  mutate(drc_sqrt_norm_DMSO = map2(data, params, ~ drm(.x$mean_sqrt_norm_DMSO ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_sqrt_norm_DMSO = map(drc_sqrt_norm_DMSO, glance)) %>%
  mutate(tidy_sqrt_norm_DMSO = map(drc_sqrt_norm_DMSO, tidy)) %>%
  mutate(drc_norm_min = map2(data, params, ~ drm(.x$mean_norm_min ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_norm_min = map(drc_norm_min, glance)) %>%
  mutate(tidy_norm_min = map(drc_norm_min, tidy)) %>%
  mutate(drc_sqrt_norm_min = map2(data, params, ~ drm(.x$mean_sqrt_norm_min ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_sqrt_norm_min = map(drc_sqrt_norm_min, glance)) %>%
  mutate(tidy_sqrt_norm_min = map(drc_sqrt_norm_min, tidy))


ec50_norm_DMSO <- curves %>%
  unnest(tidy_norm_DMSO) %>%
  select(-drc_norm_DMSO, -data, -glance_norm_DMSO, -drc_sqrt_norm_DMSO, -glance_sqrt_norm_DMSO, -tidy_sqrt_norm_DMSO, 
         -drc_norm_min, -glance_norm_min, -tidy_norm_min, -drc_sqrt_norm_min, -glance_sqrt_norm_min, -tidy_sqrt_norm_min) %>% 
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
  filter(curve_norm_DMSO == 'N2')

ec50_sqrt_norm_DMSO <- curves %>%
  unnest(tidy_sqrt_norm_DMSO) %>%
  select(-drc_sqrt_norm_DMSO, -data, -glance_sqrt_norm_DMSO, -drc_norm_DMSO, -glance_norm_DMSO, -tidy_norm_DMSO,
         -drc_norm_min, -glance_norm_min, -tidy_norm_min, -drc_sqrt_norm_min, -glance_sqrt_norm_min, -tidy_sqrt_norm_min) %>%
  filter(term == 'e') %>%
  select(-term, -params) %>%
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
  select(-treatment, -group) %>%
  rename(curve_sqrt_norm_DMSO = curve, estimate_sqrt_norm_DMSO = estimate, std.error_sqrt_norm_DMSO = std.error,
          statistic_sqrt_norm_DMSO = statistic, p.value_sqrt_norm_DMSO = p.value, xmax_sqrt_norm_DMSO = xmax, xmin_sqrt_norm_DMSO = xmin) %>%
  filter(curve_sqrt_norm_DMSO == 'N2')
   
ec50_norm_min <- curves %>%
  unnest(tidy_norm_min) %>%
  select(-drc_norm_min, -data, -glance_norm_min, -drc_norm_DMSO, -glance_norm_DMSO, -tidy_norm_DMSO, -drc_sqrt_norm_DMSO,
         -glance_sqrt_norm_DMSO, -tidy_sqrt_norm_DMSO, -drc_sqrt_norm_min, -glance_sqrt_norm_min, -tidy_sqrt_norm_min) %>%
  filter(term == 'e') %>%
  select(-term, -params) %>%
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
  select(-treatment, -group) %>%
  rename(curve_norm_min = curve, estimate_norm_min = estimate, std.error_norm_min = std.error,
         statistic_norm_min = statistic, p.value_norm_min = p.value, xmax_norm_min = xmax, xmin_norm_min = xmin) %>%
  filter(curve_norm_min == 'N2')

ec50_sqrt_norm_min <- curves %>%
  unnest(tidy_sqrt_norm_min) %>%
  select(-drc_sqrt_norm_min, -data, -glance_sqrt_norm_min, -drc_norm_min, -tidy_norm_min, -glance_norm_min, -drc_norm_DMSO, -glance_norm_DMSO, -tidy_norm_DMSO, -drc_sqrt_norm_DMSO,
         -glance_sqrt_norm_DMSO, -tidy_sqrt_norm_DMSO) %>%
  filter(term == 'e') %>%
  select(-term, -params) %>%
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
  select(-treatment, -group) %>%
  rename(curve_sqrt_norm_min = curve, estimate_sqrt_norm_min = estimate, std.error_sqrt_norm_min = std.error,
         statistic_sqrt_norm_min = statistic, p.value_sqrt_norm_min = p.value, xmax_sqrt_norm_min = xmax, xmin_sqrt_norm_min = xmin) %>%
  filter(curve_sqrt_norm_min == 'N2')

EC50s_combined <- cbind(ec50_norm_DMSO, ec50_sqrt_norm_DMSO, ec50_norm_min, ec50_sqrt_norm_min) 
  

ec50_norm_DMSO_all <- curves_all %>%
  unnest(tidy_norm_DMSO) %>%
  filter(curve == 'N2', curve == 'N2') %>%
  select(-drc_norm_DMSO, -data, -glance_norm_DMSO, -drc_sqrt_norm_DMSO, -glance_sqrt_norm_DMSO, -tidy_sqrt_norm_DMSO, 
         -drc_norm_min, -glance_norm_min, -tidy_norm_min, -drc_sqrt_norm_min, -glance_sqrt_norm_min, -tidy_sqrt_norm_min) %>% 
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
         statistic_norm_DMSO = statistic, p.value_norm_DMSO = p.value, xmax_norm_DMSO = xmax, xmin_norm_DMSO = xmin) 


ec50_sqrt_norm_DMSO_all <- curves_all %>%
  unnest(tidy_sqrt_norm_DMSO) %>%
  select(-drc_sqrt_norm_DMSO, -data, -glance_sqrt_norm_DMSO, -drc_norm_DMSO, -glance_norm_DMSO, -tidy_norm_DMSO,
         -drc_norm_min, -glance_norm_min, -tidy_norm_min, -drc_sqrt_norm_min, -glance_sqrt_norm_min, -tidy_sqrt_norm_min) %>%
  filter(term == 'e', curve == 'N2') %>%
  select(-term, -params, -curve) %>%
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
  select(-treatment) %>%
  rename( estimate_sqrt_norm_DMSO = estimate, std.error_sqrt_norm_DMSO = std.error,
          statistic_sqrt_norm_DMSO = statistic, p.value_sqrt_norm_DMSO = p.value, xmax_sqrt_norm_DMSO = xmax, xmin_sqrt_norm_DMSO = xmin)

ec50_norm_min_all <- curves_all %>%
  unnest(tidy_norm_min) %>%
  select(-drc_norm_min, -data, -glance_norm_min, -drc_norm_DMSO, -glance_norm_DMSO, -tidy_norm_DMSO, -drc_sqrt_norm_DMSO,
         -glance_sqrt_norm_DMSO, -tidy_sqrt_norm_DMSO, -drc_sqrt_norm_min, -glance_sqrt_norm_min, -tidy_sqrt_norm_min) %>%
  filter(term == 'e', curve == 'N2') %>%
  select(-term, -params, -curve) %>%
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
  select(-treatment) %>%
  rename(estimate_norm_min = estimate, std.error_norm_min = std.error,
         statistic_norm_min = statistic, p.value_norm_min = p.value, xmax_norm_min = xmax, xmin_norm_min = xmin)

ec50_sqrt_norm_min_all <- curves_all %>%
  unnest(tidy_sqrt_norm_min) %>%
  select(-drc_sqrt_norm_min, -data, -glance_sqrt_norm_min, -drc_norm_min, -tidy_norm_min, -glance_norm_min, -drc_norm_DMSO, -glance_norm_DMSO, -tidy_norm_DMSO, -drc_sqrt_norm_DMSO,
         -glance_sqrt_norm_DMSO, -tidy_sqrt_norm_DMSO) %>%
  filter(term == 'e', curve == 'N2') %>%
  select(-term, -params, -curve) %>%
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
  select(-treatment) %>%
  rename(estimate_sqrt_norm_min = estimate, std.error_sqrt_norm_min = std.error,
         statistic_sqrt_norm_min = statistic, p.value_sqrt_norm_min = p.value, xmax_sqrt_norm_min = xmax, xmin_sqrt_norm_min = xmin)

EC50s_combined_all <- cbind(ec50_norm_DMSO_all, ec50_sqrt_norm_DMSO_all, ec50_norm_min_all, ec50_sqrt_norm_min_all) %>% select(treatment, estimate_norm_DMSO, estimate_sqrt_norm_DMSO, estimate_norm_min, estimate_sqrt_norm_min)

terms_norm_DMSO <- curves %>% 
  unnest(tidy_norm_DMSO) %>% 
  select(treatment, curve, term, estimate, group) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  rename(strains = curve, b_norm_DMSO = b, c_norm_DMSO = c, d_norm_DMSO = d, e_norm_DMSO = e)

terms_sqrt_norm_DMSO <- curves %>%
  unnest(tidy_sqrt_norm_DMSO) %>%
  select(treatment, curve, term, estimate, group) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(b_sqrt_norm_DMSO = b, c_sqrt_norm_DMSO = c, d_sqrt_norm_DMSO = d, e_sqrt_norm_DMSO = e)
 
terms_norm_min <- curves %>%
  unnest(tidy_norm_min) %>%
  select(treatment, curve, term, estimate, group) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename( b_norm_min = b, c_norm_min = c, d_norm_min = d, e_norm_min = e)

terms_sqrt_norm_min <- curves %>%
  unnest(tidy_sqrt_norm_min) %>%
  select(treatment, curve, term, estimate, group) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(b_sqrt_norm_min = b, c_sqrt_norm_min = c, d_sqrt_norm_min = d, e_sqrt_norm_min = e)

terms_norm_DMSO_all <- curves_all %>% 
  unnest(tidy_norm_DMSO) %>% 
  select(treatment, curve, term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  rename(strains = curve, b_norm_DMSO = b, c_norm_DMSO = c, d_norm_DMSO = d, e_norm_DMSO = e)

terms_sqrt_norm_DMSO_all <- curves_all %>%
  unnest(tidy_sqrt_norm_DMSO) %>%
  select(treatment, curve, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  select(-treatment, -curve) %>%
  rename(b_sqrt_norm_DMSO = b, c_sqrt_norm_DMSO = c, d_sqrt_norm_DMSO = d, e_sqrt_norm_DMSO = e)

terms_norm_min_all <- curves_all %>%
  unnest(tidy_norm_min) %>%
  select(treatment, curve, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  select(-treatment, -curve) %>%
  rename( b_norm_min = b, c_norm_min = c, d_norm_min = d, e_norm_min = e)

terms_sqrt_norm_min_all <- curves_all %>%
  unnest(tidy_sqrt_norm_min) %>%
  select(treatment, curve, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  select(-treatment, -curve) %>%
  rename(b_sqrt_norm_min = b, c_sqrt_norm_min = c, d_sqrt_norm_min = d, e_sqrt_norm_min = e)

terms <- cbind(terms_norm_DMSO_all, terms_sqrt_norm_DMSO_all, terms_norm_min_all, terms_sqrt_norm_min_all)

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


predict_norm_DMSO <- newdata %>%
  left_join(terms_norm_DMSO) %>% 
  filter(strains == 'N2') %>%
  # mutate(
  #   c = case_when(
  #     treatment == 'AlbendazoleSulfoxide' ~ sqrt(0.6),
  #     treatment == 'Levamisole' ~ sqrt(0.5),
  #     TRUE ~ c
  #   )) %>% 
  # manually predict using the LL.4 model and the tidied terms
  mutate(pred_norm_DMSO = c_norm_DMSO + ((d_norm_DMSO - c_norm_DMSO) / (1 + exp(b_norm_DMSO*(log(conc) - log(e_norm_DMSO)))))) %>%
  select(-b_norm_DMSO, -c_norm_DMSO, -d_norm_DMSO) %>%
  select(treatment, strains, genotype, conc, pred_norm_DMSO, group) %>%
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  )) 

predict_sqrt_norm_DMSO <- newdata %>%
  left_join(terms_sqrt_norm_DMSO) %>% 
  filter(curve == 'N2') %>%
  # mutate(
  #   c = case_when(
  #     treatment == 'AlbendazoleSulfoxide' ~ sqrt(0.6),
  #     treatment == 'Levamisole' ~ sqrt(0.5),
  #     TRUE ~ c
  #   )) %>% 
  # manually predict using the LL.4 model and the tidied terms
  mutate(pred_sqrt_norm_DMSO = c_sqrt_norm_DMSO + ((d_sqrt_norm_DMSO - c_sqrt_norm_DMSO) / (1 + exp(b_sqrt_norm_DMSO*(log(conc) - log(e_sqrt_norm_DMSO)))))) %>%
  select(-b_sqrt_norm_DMSO, -c_sqrt_norm_DMSO, -d_sqrt_norm_DMSO) %>%
  select(treatment, strains, genotype, conc, pred_sqrt_norm_DMSO, group) %>%
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  )) 

predict_sqrt_norm_min <- newdata %>%
  left_join(terms_sqrt_norm_min) %>% 
  filter(curve == 'N2') %>%
  # mutate(
  #   c = case_when(
  #     treatment == 'AlbendazoleSulfoxide' ~ sqrt(0.6),
  #     treatment == 'Levamisole' ~ sqrt(0.5),
  #     TRUE ~ c
  #   )) %>% 
  # manually predict using the LL.4 model and the tidied terms
  mutate(pred_sqrt_norm_min = c_sqrt_norm_min + ((d_sqrt_norm_min - c_sqrt_norm_min) / (1 + exp(b_sqrt_norm_min*(log(conc) - log(e_sqrt_norm_min)))))) %>%
  select(-b_sqrt_norm_min, -c_sqrt_norm_min, -d_sqrt_norm_min) %>%
  select(treatment, strains, genotype, conc, pred_sqrt_norm_min, group) %>%
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  )) 

predict_norm_min <- newdata %>%
  left_join(terms_norm_min) %>% 
  filter(curve == 'N2') %>%
  # mutate(
  #   c = case_when(
  #     treatment == 'AlbendazoleSulfoxide' ~ sqrt(0.6),
  #     treatment == 'Levamisole' ~ sqrt(0.5),
  #     TRUE ~ c
  #   )) %>% 
  # manually predict using the LL.4 model and the tidied terms
  mutate(pred_norm_min = c_norm_min + ((d_norm_min - c_norm_min) / (1 + exp(b_norm_min*(log(conc) - log(e_norm_min)))))) %>%
  select(-b_norm_min, -c_norm_min, -d_norm_min) %>%
  select(treatment, strains, genotype, conc, pred_norm_min, group) %>%
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  )) 

predict_all <- newdata %>%
  left_join(terms) %>%
  # mutate(
  #   c = case_when(
  #     treatment == 'AlbendazoleSulfoxide' ~ sqrt(0.6),
  #     treatment == 'Levamisole' ~ sqrt(0.5),
  #     TRUE ~ c
  #   )) %>%
  # manually predict using the LL.4 model and the tidied terms
  mutate(pred_norm_DMSO = c_norm_DMSO + ((d_norm_DMSO - c_norm_DMSO) / (1 + exp(b_norm_DMSO*(log(conc) - log(e_norm_DMSO)))))) %>%
  select(-b_norm_DMSO, -c_norm_DMSO, -d_norm_DMSO) %>%
  mutate(pred_sqrt_norm_DMSO = c_sqrt_norm_DMSO + ((d_sqrt_norm_DMSO - c_sqrt_norm_DMSO) / (1 + exp(b_sqrt_norm_DMSO*(log(conc) - log(e_sqrt_norm_DMSO)))))) %>%
  select(-b_sqrt_norm_DMSO, -c_sqrt_norm_DMSO, -d_sqrt_norm_DMSO) %>%
  mutate(pred_norm_min = c_norm_min + ((d_norm_min - c_norm_min) / (1 + exp(b_norm_min*(log(conc) - log(e_norm_min)))))) %>%
  select(-b_norm_min, -c_norm_min, -d_norm_min) %>%
  mutate(pred_sqrt_norm_min = c_sqrt_norm_min + ((d_sqrt_norm_min - c_sqrt_norm_min) / (1 + exp(b_sqrt_norm_min*(log(conc) - log(e_sqrt_norm_min)))))) %>%
  select(-b_sqrt_norm_min, -c_sqrt_norm_min, -d_sqrt_norm_min) %>%
  select(treatment, strains, conc, genotype, pred_norm_DMSO, pred_sqrt_norm_DMSO, pred_norm_min, pred_sqrt_norm_min) %>%
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  ))



#make dataframes containing ec50 values for labelling on plots

ann_text_A_norm_DMSO <- data.frame(lab = "EC[50] *=* 15.7", treatment = factor('AlbendazoleSulfoxide',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))
ann_text_I_norm_DMSO <- data.frame(lab = "EC[50] *=* 8.55", treatment = factor('Ivermectin',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))
ann_text_L_norm_DMSO <- data.frame(lab = "EC[50] *=* 6.83", treatment = factor('Levamisole',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))

ann_text_A_sqrt_norm_DMSO <- data.frame(lab = "EC[50] *=* 18.6", treatment = factor('AlbendazoleSulfoxide',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))
ann_text_I_sqrt_norm_DMSO <- data.frame(lab = "EC[50] *=* 8.73", treatment = factor('Ivermectin',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))
ann_text_L_sqrt_norm_DMSO <- data.frame(lab = "EC[50] *=* 7.73", treatment = factor('Levamisole',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))

ann_text_A_norm_min <- data.frame(lab = "EC[50] *=* 15.3", treatment = factor('AlbendazoleSulfoxide',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))
ann_text_I_norm_min <- data.frame(lab = "EC[50] *=* 8.75", treatment = factor('Ivermectin',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))
ann_text_L_norm_min <- data.frame(lab = "EC[50] *=* 6.51", treatment = factor('Levamisole',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))

ann_text_A_sqrt_norm_min <- data.frame(lab = "EC[50] *=* 19.5", treatment = factor('AlbendazoleSulfoxide',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))
ann_text_I_sqrt_norm_min <- data.frame(lab = "EC[50] *=* 9.03", treatment = factor('Ivermectin',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))
ann_text_L_sqrt_norm_min <- data.frame(lab = "EC[50] *=* 7.20", treatment = factor('Levamisole',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))

#plot all 4 normalization schemes
#figure 1B
(n2_plot_norm_DMSO <- revised_data %>%
    filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
           strains == 'N2') %>% 
    ggplot() +
    geom_line(data = predict_norm_DMSO %>% filter(strains == 'N2', group == '1'), aes(x = conc, group = genotype, y = pred_norm_DMSO),linewidth = 0.5, color = 'paleturquoise4', alpha =0.7) +
    geom_line(data = predict_norm_DMSO %>% filter(strains == 'N2', group == '2'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'slateblue4', alpha =0.7) +
    geom_line(data = predict_norm_DMSO %>% filter(strains == 'N2', group == '3'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'darkorange4', alpha =0.7) +
    geom_line(data = predict_norm_DMSO %>% filter(strains == 'N2', group == '4'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'chocolate', alpha =0.7) +
    geom_line(data = predict_norm_DMSO %>% filter(strains == 'N2', group == '5'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'chartreuse4', alpha =0.7) +
    geom_line(data = predict_norm_DMSO %>% filter(strains == 'N2', group == '6'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'cornflowerblue', alpha =0.7) +
    geom_line(data = predict_norm_DMSO %>% filter(strains == 'N2', group == '7'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'darkorchid', alpha =0.7) +
    geom_line(data = predict_norm_DMSO %>% filter(strains == 'N2', group == '8'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'rosybrown4', alpha =0.7) +
    
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 1),
                     aes(x = conc, y = mean_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'paleturquoise4', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 2),
                     aes(x = conc, y = mean_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'slateblue4', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 3),
                     aes(x = conc, y = mean_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'darkorange4', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 4),
                     aes(x = conc, y = mean_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'chocolate', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 5),
                     aes(x = conc, y = mean_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'chartreuse4', shape = 16) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 6),
                     aes(x = conc, y = mean_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'cornflowerblue', shape = 16) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 7),
                     aes(x = conc, y = mean_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'darkorchid', shape = 16) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 8),
                     aes(x = conc, y = mean_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'rosybrown4', shape = 16) +
    geom_vline(data = EC50s_combined %>% filter(group == 1), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'paleturquoise4', alpha = 0.5) + 
    geom_vline(data = EC50s_combined %>% filter(group == 2), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'slateblue4', alpha = 0.5) + 
    geom_vline(data = EC50s_combined %>% filter(group == 3), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'darkorange4', alpha = 0.5) + #AZS only
    geom_vline(data = EC50s_combined %>% filter(group == 4), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'chocolate', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 5), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'chartreuse4', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 6), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'cornflowerblue', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 7), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'darkorchid', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 8), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'rosybrown4', alpha = 0.5) +
    geom_vline(data = EC50s_combined_all, aes(xintercept = estimate_norm_DMSO, group = treatment), size = 1, linetype = 'dashed', color = 'black', alpha = 1) +
    
    geom_line(data = predict_all %>% filter(strains == 'N2'), aes(x = conc, y = pred_norm_DMSO), color = 'black', size = 1, alpha =1) + 
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
    annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scale_y_continuous(limits = c(0, 1.25), expand = c(0, 0)) +
    # scale_color_manual(values = c('paleturquoise4', 'slateblue4', 'paleturquoise4', 'slateblue4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4')) + 
    # scale_fill_manual(values = c('paleturquoise4', 'slateblue4', 'paleturquoise4', 'slateblue4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4')) +
    scale_shape(guide = 'none') +
    coord_cartesian(clip = "off") +
    
    labs(x = 'Drug concentration (µM)', y = "Normalized length",
         color = 'Strain', fill = 'Strain', tag = "A") +
    facet_grid(cols = vars(treatment), scales = 'free_x',
               labeller = as_labeller(c(AlbendazoleSulfoxide = 'AZS',
                                        Ivermectin = 'IVM',
                                        Levamisole = 'LEV'))) +
    theme_nw2() +
    geom_richtext(data = ann_text_A_norm_DMSO, label = "EC<sub>50</sub> = 15.7 µM", aes(x = 0.8, y = 0.05), size = 3) + #value from ec50_norm_DMSO_all, albendazolesulfoxide
    geom_richtext(data = ann_text_I_norm_DMSO, label = "EC<sub>50</sub> = 8.55 nM", aes(x = 0.001, y = 0.05), size = 3) + #value from ec50_norm_DMSO_all, ivermectin
    geom_richtext(data = ann_text_L_norm_DMSO, label = "EC<sub>50</sub> = 6.83 µM", aes(x = 0.1, y = 0.05), size = 3) + #value from ec50_norm_DMSO_all, levamisole
    theme(
      legend.position = 'empty',
      axis.text.x = element_markdown(angle = 0, hjust = 0.5),
      axis.title.x = element_markdown(face = 'plain'),
      axis.title.y = element_markdown(face = 'plain'),
      plot.title = element_text(hjust = 0.5),
      plot.tag = element_text(size = 10),
      strip.text = element_text(color = ('black'), size = 15),
      plot.margin = unit(c(0.5,0.5,0.5,2), 'lines') 
    ) +
    NULL)


(n2_plot_sqrt_norm_DMSO <- revised_data %>%
    filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
           strains == 'N2') %>% 
    ggplot() +
    geom_line(data = predict_sqrt_norm_DMSO %>% filter(strains == 'N2', group == '1'), aes(x = conc, group = treatment, y = pred_sqrt_norm_DMSO),linewidth = 0.5, color = 'paleturquoise4', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_DMSO %>% filter(strains == 'N2', group == '2'), aes(x = conc, group = genotype, y = pred_sqrt_norm_DMSO), linewidth = 0.5, color = 'slateblue4', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_DMSO %>% filter(strains == 'N2', group == '3'), aes(x = conc, group = genotype, y = pred_sqrt_norm_DMSO), linewidth = 0.5, color = 'darkorange4', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_DMSO %>% filter(strains == 'N2', group == '4'), aes(x = conc, group = genotype, y = pred_sqrt_norm_DMSO), linewidth = 0.5, color = 'chocolate', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_DMSO %>% filter(strains == 'N2', group == '5'), aes(x = conc, group = genotype, y = pred_sqrt_norm_DMSO), linewidth = 0.5, color = 'chartreuse4', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_DMSO %>% filter(strains == 'N2', group == '6'), aes(x = conc, group = genotype, y = pred_sqrt_norm_DMSO), linewidth = 0.5, color = 'cornflowerblue', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_DMSO %>% filter(strains == 'N2', group == '7'), aes(x = conc, group = genotype, y = pred_sqrt_norm_DMSO), linewidth = 0.5, color = 'darkorchid', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_DMSO %>% filter(strains == 'N2', group == '8'), aes(x = conc, group = genotype, y = pred_sqrt_norm_DMSO), linewidth = 0.5, color = 'rosybrown4', alpha =0.7) +
    
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 1),
                     aes(x = conc, y = mean_sqrt_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'paleturquoise4', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 2),
                     aes(x = conc, y = mean_sqrt_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'slateblue4', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 3),
                     aes(x = conc, y = mean_sqrt_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'darkorange4', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 4),
                     aes(x = conc, y = mean_sqrt_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'chocolate', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 5),
                     aes(x = conc, y = mean_sqrt_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'chartreuse4', shape = 16) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 6),
                     aes(x = conc, y = mean_sqrt_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'cornflowerblue', shape = 16) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 7),
                     aes(x = conc, y = mean_sqrt_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'darkorchid', shape = 16) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 8),
                     aes(x = conc, y = mean_sqrt_norm_DMSO, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'rosybrown4', shape = 16) +
    geom_vline(data = EC50s_combined %>% filter(group == 1), aes(xintercept = estimate_sqrt_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'paleturquoise4', alpha = 0.5) + 
    geom_vline(data = EC50s_combined %>% filter(group == 2), aes(xintercept = estimate_sqrt_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'slateblue4', alpha = 0.5) + 
    geom_vline(data = EC50s_combined %>% filter(group == 3), aes(xintercept = estimate_sqrt_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'darkorange4', alpha = 0.5) + #AZS only
    geom_vline(data = EC50s_combined %>% filter(group == 4), aes(xintercept = estimate_sqrt_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'chocolate', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 5), aes(xintercept = estimate_sqrt_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'chartreuse4', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 6), aes(xintercept = estimate_sqrt_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'cornflowerblue', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 7), aes(xintercept = estimate_sqrt_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'darkorchid', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 8), aes(xintercept = estimate_sqrt_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'rosybrown4', alpha = 0.5) +
    geom_vline(data = EC50s_combined_all, aes(xintercept = estimate_sqrt_norm_DMSO, group = treatment), linetype = 'dashed', size = 0.5, color = 'black', alpha = 1) +
    #geom_vline(data = EC50s_combined_all %>% filter(s), aes(xintercept = estimate_sqrt_norm_DMSO, group = treatment), size = 1, linetype = 'dashed', color = 'black', alpha = 0.75) +
    
    geom_line(data = predict_all %>% filter(strains == 'N2'), aes(x = conc, y = pred_sqrt_norm_DMSO), color = 'black', size = 1, alpha =1) + 
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
    annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scale_y_continuous(limits = c(0, 1.25), expand = c(0, 0)) +
    # scale_color_manual(values = c('paleturquoise4', 'slateblue4', 'paleturquoise4', 'slateblue4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4')) + 
    # scale_fill_manual(values = c('paleturquoise4', 'slateblue4', 'paleturquoise4', 'slateblue4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4')) +
    scale_shape(guide = 'none') +
    coord_cartesian(clip = "off") +
    
    labs(x = 'Drug concentration (µM)', y = "Normalized length",
         color = 'Strain', fill = 'Strain', tag = "C") +
    facet_grid(cols = vars(treatment), scales = 'free_x',
               labeller = as_labeller(c(AlbendazoleSulfoxide = 'AZS',
                                        Ivermectin = 'IVM',
                                        Levamisole = 'LEV'))) +
    theme_nw2() +
    theme_nw2() +
    geom_richtext(data = ann_text_A_sqrt_norm_DMSO, label = "EC<sub>50</sub> = 18.6 µM", aes(x = 0.8, y = 0.05), size = 3) + #value from ec50_sqrt_norm_DMSO_all, albendazolesulfoxide
    geom_richtext(data = ann_text_I_sqrt_norm_DMSO, label = "EC<sub>50</sub> = 8.73 nM", aes(x = 0.001, y = 0.05), size = 3) + #value from ec50_sqrt_norm_DMSO_all, ivermectin
    geom_richtext(data = ann_text_L_sqrt_norm_DMSO, label = "EC<sub>50</sub> = 7.73 µM", aes(x = 0.1, y = 0.05), size = 3) + #value from ec50_sqrt_norm_DMSO_all, levamisole
    theme(
      legend.position = 'empty',
      axis.text.x = element_markdown(angle = 0, hjust = 0.5),
      axis.title.x = element_markdown(face = 'plain'),
      axis.title.y = element_markdown(face = 'plain'),
      plot.title = element_text(hjust = 0.5),
      plot.tag = element_text(size = 10),
      strip.text = element_text(color = ('black'), size = 15),
      plot.margin = unit(c(0.5,0.5,0.5,2), 'lines') 
    ) +
    NULL)


 
#ggsave(here('/Users/elenagr/Desktop/Fig_1_supp_SQ_DMSO.pdf'), n2_plot_sqrt_norm_DMSO, width = 240, height = 80, units = 'mm')
#ggsave(here('/Users/elenagr/Desktop/Fig_1_supp_SQ_DMSO.png'), n2_plot_sqrt_norm_DMSO, width = 240, height = 80, units = 'mm')

(n2_plot_norm_min <- revised_data %>%
    filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
           strains == 'N2') %>% 
    ggplot() +
    geom_line(data = predict_norm_min %>% filter(strains == 'N2', group == '1'), aes(x = conc, group = treatment, y = pred_norm_min),linewidth = 0.5, color = 'paleturquoise4', alpha =0.7) +
    geom_line(data = predict_norm_min %>% filter(strains == 'N2', group == '2'), aes(x = conc, group = genotype, y = pred_norm_min), linewidth = 0.5, color = 'slateblue4', alpha =0.7) +
    geom_line(data = predict_norm_min %>% filter(strains == 'N2', group == '3'), aes(x = conc, group = genotype, y = pred_norm_min), linewidth = 0.5, color = 'darkorange4', alpha =0.7) +
    geom_line(data = predict_norm_min %>% filter(strains == 'N2', group == '4'), aes(x = conc, group = genotype, y = pred_norm_min), linewidth = 0.5, color = 'chocolate', alpha =0.7) +
    geom_line(data = predict_norm_min %>% filter(strains == 'N2', group == '5'), aes(x = conc, group = genotype, y = pred_norm_min), linewidth = 0.5, color = 'chartreuse4', alpha =0.7) +
    geom_line(data = predict_norm_min %>% filter(strains == 'N2', group == '6'), aes(x = conc, group = genotype, y = pred_norm_min), linewidth = 0.5, color = 'cornflowerblue', alpha =0.7) +
    geom_line(data = predict_norm_min %>% filter(strains == 'N2', group == '7'), aes(x = conc, group = genotype, y = pred_norm_min), linewidth = 0.5, color = 'darkorchid', alpha =0.7) +
    geom_line(data = predict_norm_min %>% filter(strains == 'N2', group == '8'), aes(x = conc, group = genotype, y = pred_norm_min), linewidth = 0.5, color = 'rosybrown4', alpha =0.7) +
    
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 1),
                     aes(x = conc, y = mean_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'paleturquoise4', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 2),
                     aes(x = conc, y = mean_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'slateblue4', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 3),
                     aes(x = conc, y = mean_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'darkorange4', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 4),
                     aes(x = conc, y = mean_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'chocolate', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 5),
                     aes(x = conc, y = mean_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'chartreuse4', shape = 16) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 6),
                     aes(x = conc, y = mean_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'cornflowerblue', shape = 16) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 7),
                     aes(x = conc, y = mean_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'darkorchid', shape = 16) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 8),
                     aes(x = conc, y = mean_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'rosybrown4', shape = 16) +
    geom_vline(data = EC50s_combined %>% filter(group == 1), aes(xintercept = estimate_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'paleturquoise4', alpha = 0.5) + 
    geom_vline(data = EC50s_combined %>% filter(group == 2), aes(xintercept = estimate_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'slateblue4', alpha = 0.5) + 
    geom_vline(data = EC50s_combined %>% filter(group == 3), aes(xintercept = estimate_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'darkorange4', alpha = 0.5) + #AZS only
    geom_vline(data = EC50s_combined %>% filter(group == 4), aes(xintercept = estimate_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'chocolate', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 5), aes(xintercept = estimate_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'chartreuse4', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 6), aes(xintercept = estimate_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'cornflowerblue', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 7), aes(xintercept = estimate_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'darkorchid', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 8), aes(xintercept = estimate_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'rosybrown4', alpha = 0.5) +
    geom_vline(data = EC50s_combined_all, aes(xintercept = estimate_norm_min, group = treatment), linetype = 'dashed', size = 0.5, color = 'black', alpha = 1) +
    #geom_vline(data = EC50s_combined_all %>% filter(s), aes(xintercept = estimate_sqrt_norm_DMSO, group = treatment), size = 1, linetype = 'dashed', color = 'black', alpha = 0.75) +
    
    geom_line(data = predict_all %>% filter(strains == 'N2'), aes(x = conc, y = pred_norm_min), color = 'black', size = 1, alpha =1) + 
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
    annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scale_y_continuous(limits = c(0, 1.25), expand = c(0, 0)) +
    # scale_color_manual(values = c('paleturquoise4', 'slateblue4', 'paleturquoise4', 'slateblue4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4')) + 
    # scale_fill_manual(values = c('paleturquoise4', 'slateblue4', 'paleturquoise4', 'slateblue4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4')) +
    scale_shape(guide = 'none') +
    coord_cartesian(clip = "off") +
    
    labs(x = 'Drug concentration (µM)', y = "Normalized length",
         color = 'Strain', fill = 'Strain', tag = "B") +
    facet_grid(cols = vars(treatment), scales = 'free_x',
               labeller = as_labeller(c(AlbendazoleSulfoxide = 'AZS',
                                        Ivermectin = 'IVM',
                                        Levamisole = 'LEV'))) +
    theme_nw2() +
    geom_richtext(data = ann_text_A_norm_min, label = "EC<sub>50</sub> = 15.3 µM", aes(x = 0.8, y = 0.05), size = 3) + #value from ec50_norm_min_all, albendazolesulfoxide
    geom_richtext(data = ann_text_I_norm_min, label = "EC<sub>50</sub> = 8.75 nM", aes(x = 0.001, y = 0.05), size = 3) + #value from ec50_norm_min_all, ivermectin
    geom_richtext(data = ann_text_L_norm_min, label = "EC<sub>50</sub> = 6.51 µM", aes(x = 0.1, y = 0.05), size = 3) + #value from ec50_norm_minall, levamisole
    theme(
      legend.position = 'empty',
      axis.text.x = element_markdown(angle = 0, hjust = 0.5),
      axis.title.x = element_markdown(face = 'plain'),
      axis.title.y = element_markdown(face = 'plain'),
      plot.title = element_text(hjust = 0.5),
      plot.tag = element_text(size = 10),
      strip.text = element_text(color = ('black'), size = 15),
      plot.margin = unit(c(0.5,0.5,0.5,2), 'lines') 
    ) +
    NULL)

#ggsave(here('/Users/elenagr/Desktop/Fig_1_supp_min.pdf'), n2_plot_norm_min, width = 240, height = 80, units = 'mm')
#ggsave(here('/Users/elenagr/Desktop/Fig_1_supp_min.png'), n2_plot_norm_min, width = 240, height = 80, units = 'mm')

(n2_plot_sqrt_norm_min <- revised_data %>%
    filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
           strains == 'N2') %>% 
    ggplot() +
    geom_line(data = predict_sqrt_norm_min %>% filter(strains == 'N2', group == '1'), aes(x = conc, group = treatment, y = pred_sqrt_norm_min),linewidth = 0.5, color = 'paleturquoise4', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_min %>% filter(strains == 'N2', group == '2'), aes(x = conc, group = genotype, y = pred_sqrt_norm_min), linewidth = 0.5, color = 'slateblue4', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_min %>% filter(strains == 'N2', group == '3'), aes(x = conc, group = genotype, y = pred_sqrt_norm_min), linewidth = 0.5, color = 'darkorange4', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_min %>% filter(strains == 'N2', group == '4'), aes(x = conc, group = genotype, y = pred_sqrt_norm_min), linewidth = 0.5, color = 'chocolate', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_min %>% filter(strains == 'N2', group == '5'), aes(x = conc, group = genotype, y = pred_sqrt_norm_min), linewidth = 0.5, color = 'chartreuse4', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_min %>% filter(strains == 'N2', group == '6'), aes(x = conc, group = genotype, y = pred_sqrt_norm_min), linewidth = 0.5, color = 'cornflowerblue', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_min %>% filter(strains == 'N2', group == '7'), aes(x = conc, group = genotype, y = pred_sqrt_norm_min), linewidth = 0.5, color = 'darkorchid', alpha =0.7) +
    geom_line(data = predict_sqrt_norm_min %>% filter(strains == 'N2', group == '8'), aes(x = conc, group = genotype, y = pred_sqrt_norm_min), linewidth = 0.5, color = 'rosybrown4', alpha =0.7) +
    
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 1),
                     aes(x = conc, y = mean_sqrt_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'paleturquoise4', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 2),
                     aes(x = conc, y = mean_sqrt_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'slateblue4', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 3),
                     aes(x = conc, y = mean_sqrt_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'darkorange4', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 4),
                     aes(x = conc, y = mean_sqrt_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'chocolate', shape = 17) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 5),
                     aes(x = conc, y = mean_sqrt_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'chartreuse4', shape = 16) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 6),
                     aes(x = conc, y = mean_sqrt_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'cornflowerblue', shape = 16) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 7),
                     aes(x = conc, y = mean_sqrt_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'darkorchid', shape = 16) +
    geom_quasirandom(data = well_summary %>% filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'), strains == 'N2', group == 8),
                     aes(x = conc, y = mean_sqrt_norm_min, color = plate), size = 1, width = 0.1, alpha = 0.5, color = 'rosybrown4', shape = 16) +
    geom_vline(data = EC50s_combined %>% filter(group == 1), aes(xintercept = estimate_sqrt_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'paleturquoise4', alpha = 0.5) + 
    geom_vline(data = EC50s_combined %>% filter(group == 2), aes(xintercept = estimate_sqrt_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'slateblue4', alpha = 0.5) + 
    geom_vline(data = EC50s_combined %>% filter(group == 3), aes(xintercept = estimate_sqrt_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'darkorange4', alpha = 0.5) + #AZS only
    geom_vline(data = EC50s_combined %>% filter(group == 4), aes(xintercept = estimate_sqrt_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'chocolate', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 5), aes(xintercept = estimate_sqrt_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'chartreuse4', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 6), aes(xintercept = estimate_sqrt_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'cornflowerblue', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 7), aes(xintercept = estimate_sqrt_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'darkorchid', alpha = 0.5) +
    geom_vline(data = EC50s_combined %>% filter(group == 8), aes(xintercept = estimate_sqrt_norm_min, group = treatment), size = 0.5, linetype = 'dashed', color = 'rosybrown4', alpha = 0.5) +
    geom_vline(data = EC50s_combined_all, aes(xintercept = estimate_sqrt_norm_min, group = treatment), linetype = 'dashed', size = 0.5, color = 'black', alpha = 1) +
    #geom_vline(data = EC50s_combined_all %>% filter(s), aes(xintercept = estimate_sqrt_norm_DMSO, group = treatment), size = 1, linetype = 'dashed', color = 'black', alpha = 0.75) +
    
    geom_line(data = predict_all %>% filter(strains == 'N2'), aes(x = conc, y = pred_sqrt_norm_min), color = 'black', size = 1, alpha = 1) + 
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
    annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scale_y_continuous(limits = c(0, 1.25), expand = c(0, 0)) +
    # scale_color_manual(values = c('paleturquoise4', 'slateblue4', 'paleturquoise4', 'slateblue4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4')) + 
    # scale_fill_manual(values = c('paleturquoise4', 'slateblue4', 'paleturquoise4', 'slateblue4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4')) +
    scale_shape(guide = 'none') +
    coord_cartesian(clip = "off") +
    
    labs(x = 'Drug concentration (µM)', y = "Normalized length",
         color = 'Strain', fill = 'Strain', tag = "D") +
    facet_grid(cols = vars(treatment), scales = 'free_x',
               labeller = as_labeller(c(AlbendazoleSulfoxide = 'AZS',
                                        Ivermectin = 'IVM',
                                        Levamisole = 'LEV'))) +
    theme_nw2() +
    theme_nw2() +
    geom_richtext(data = ann_text_A_sqrt_norm_min, label = "EC<sub>50</sub> = 19.5 µM", aes(x = 0.8, y = 0.05), size = 3) + #value from ec50_sqrt_norm_min_all, albendazolesulfoxide
    geom_richtext(data = ann_text_I_sqrt_norm_min, label = "EC<sub>50</sub> = 9.03 nM", aes(x = 0.001, y = 0.05), size = 3) + #value from ec50_sqrt_norm_min_all, ivermectin
    geom_richtext(data = ann_text_L_sqrt_norm_min, label = "EC<sub>50</sub> = 7.20 µM", aes(x = 0.1, y = 0.05), size = 3) + #value from ec50_nsqrt_orm_min_all, levamisole
    theme(
      legend.position = 'empty',
      axis.text.x = element_markdown(angle = 0, hjust = 0.5),
      axis.title.x = element_markdown(face = 'plain'),
      axis.title.y = element_markdown(face = 'plain'),
      plot.title = element_text(hjust = 0.5),
      plot.tag = element_text(size = 10),
      strip.text = element_text(color = ('black'), size = 15),
      plot.margin = unit(c(0.5,0.5,0.5,2), 'lines') 
    ) +
    NULL)
                                       
#ggsave(here('/Users/elenagr/Desktop/Fig_1_supp_SQ_min.pdf'), n2_plot_sqrt_norm_min, width = 240, height = 80, units = 'mm')
#ggsave(here('/Users/elenagr/Desktop/Fig_1_supp_SQ_min.png'), n2_plot_sqrt_norm_min, width = 240, height = 80, units = 'mm')


(combined_plot <- plot_grid(n2_plot_norm_DMSO,
          n2_plot_norm_min,
          n2_plot_sqrt_norm_DMSO,
          n2_plot_sqrt_norm_min,
          nrow = 2, ncol = 2, align = 'h', axis = 'l', rel_heights = c(1, 1, 1.2), base_height = 10, base_width = 14))

ggsave(here('/Users/elenagarncarz/Desktop/Supp_Fig_1_updated.pdf'), combined_plot, width = 500, height = 300, units = 'mm')
ggsave(here('/Users/elenagarncarz/Desktop/Supp_Fig_1_updated.png'), combined_plot, width = 500, height = 300, units = 'mm')
