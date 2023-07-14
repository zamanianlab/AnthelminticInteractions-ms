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

well_summary <- revised_data %>% 
  select(plate, well, conc, metadata_date, strains, treatment, genotype, raw_length, norm_DMSO) %>% # sqrt_length, min_value, sqrt_norm_DMSO, norm_min, sqrt_norm_min
  group_by(metadata_date, plate, genotype, strains, well, treatment, conc) %>% 
  summarise(mean_raw_length = mean(raw_length),
            # mean_sqrt_length = mean(sqrt_length),
            mean_norm_DMSO = mean(norm_DMSO)
            # mean_sqrt_norm_DMSO = mean(sqrt_norm_DMSO),
            # mean_norm_min = mean(norm_min),
            # mean_sqrt_norm_min = mean(sqrt_norm_min)
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


### dose-responses
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
            # mean_sqrt_length = mean(sqrt_length),
            mean_norm_DMSO = mean(norm_DMSO),
            # mean_sqrt_norm_DMSO = mean(sqrt_norm_DMSO),
            # mean_norm_min = mean(norm_min),
            # mean_sqrt_norm_min = mean(sqrt_norm_min)
            ) %>%
  ungroup() %>%
  group_nest(treatment, group) %>%
  # rowwise() %>%
  # fit the models using all the data but group by strain, give some parameters
  mutate(params = list(c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), 
                       c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), 
                       c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), 
                       c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), 
                       c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA), 
                       c(NA, NA, NA, NA), c(NA, NA, NA, NA))) %>%
  mutate(drc_norm_DMSO  = map2(data, params, ~ drm(.x$mean_norm_DMSO ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_norm_DMSO = map(drc_norm_DMSO, glance)) %>%
  mutate(tidy_norm_DMSO = map(drc_norm_DMSO, tidy))
  # mutate(group = '1')

curve_all <- revised_data %>%
  left_join(reps) %>%
  # use only the most recent 3 replicates (not just those that used the same drug preparation)
  filter(!treatment %in% c('DMSO', 'Untreated')) %>%
  # uncomment below to fit to the summarized well data
  group_by(plate, treatment, strains, genotype, conc) %>%
  summarise(mean_raw_length = mean(raw_length),
            # mean_sqrt_length = mean(sqrt_length),
            mean_norm_DMSO = mean(norm_DMSO),
            # mean_sqrt_norm_DMSO = mean(sqrt_norm_DMSO),
            # mean_norm_min = mean(norm_min),
            # mean_sqrt_norm_min = mean(sqrt_norm_min)
  ) %>%
  ungroup() %>%
  group_nest(treatment) %>%
  # rowwise() %>%
  # fit the models using all the data but group by strain, give some parameters
  mutate(params = list(c(NA, NA, NA, NA), c(NA, NA, NA, NA), c(NA, NA, NA, NA))) %>%
  mutate(drc_norm_DMSO  = map2(data, params, ~ drm(.x$mean_norm_DMSO ~ .x$conc, .x$strains, fct = LL.4(fixed = .y)))) %>%
  mutate(glance_norm_DMSO = map(drc_norm_DMSO, glance)) %>%
  mutate(tidy_norm_DMSO = map(drc_norm_DMSO, tidy))


# get the ec50
ec50_norm_DMSO <- curves %>%
  unnest(tidy_norm_DMSO) %>%
  filter(curve == 'N2') %>%
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
         statistic_norm_DMSO = statistic, p.value_norm_DMSO = p.value, xmax_norm_DMSO = xmax, xmin_norm_DMSO = xmin) 

ec50_norm_DMSO_all <- curve_all %>%
  unnest(tidy_norm_DMSO) %>%
  filter(curve == 'N2') %>%
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
         statistic_norm_DMSO = statistic, p.value_norm_DMSO = p.value, xmax_norm_DMSO = xmax, xmin_norm_DMSO = xmin) 

# terms_1 <- cbind(terms_norm_DMSO_1, terms_sqrt_norm_DMSO_1, terms_norm_min_1, terms_sqrt_norm_min_1)
# terms_2 <- cbind(terms_norm_DMSO_2, terms_sqrt_norm_DMSO_2, terms_norm_min_2, terms_sqrt_norm_min_2)
# terms_3 <- cbind(terms_norm_DMSO_3, terms_sqrt_norm_DMSO_3, terms_norm_min_3, terms_sqrt_norm_min_3)
# terms_4 <- cbind(terms_norm_DMSO_4, terms_sqrt_norm_DMSO_4, terms_norm_min_4, terms_sqrt_norm_min_4)
# terms_5 <- cbind(terms_norm_DMSO_5, terms_sqrt_norm_DMSO_5, terms_norm_min_5, terms_sqrt_norm_min_5)
# terms_6 <- cbind(terms_norm_DMSO_6, terms_sqrt_norm_DMSO_6, terms_norm_min_6, terms_sqrt_norm_min_6)
# terms_7 <- cbind(terms_norm_DMSO_7, terms_sqrt_norm_DMSO_7, terms_norm_min_7, terms_sqrt_norm_min_7)
# terms_8 <- cbind(terms_norm_DMSO_8, terms_sqrt_norm_DMSO_8, terms_norm_min_8, terms_sqrt_norm_min_8)
# terms <- cbind(terms_norm_DMSO, terms_sqrt_norm_DMSO, terms_norm_min, terms_sqrt_norm_min)

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

terms_norm_DMSO <- curves %>% 
  unnest(tidy_norm_DMSO) %>% 
  select(treatment, curve, term, estimate, group) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  rename(strains = curve, b_norm_DMSO = b, c_norm_DMSO = c, d_norm_DMSO = d, e_norm_DMSO = e)

terms_norm_DMSO_all <- curve_all %>% 
  unnest(tidy_norm_DMSO) %>% 
  select(treatment, curve, term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  rename(strains = curve, b_norm_DMSO = b, c_norm_DMSO = c, d_norm_DMSO = d, e_norm_DMSO = e)

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
  # mutate(pred_sqrt_norm_DMSO = c_sqrt_norm_DMSO + ((d_sqrt_norm_DMSO - c_sqrt_norm_DMSO) / (1 + exp(b_sqrt_norm_DMSO*(log(conc) - log(e_sqrt_norm_DMSO)))))) %>%
  # select(-b_sqrt_norm_DMSO, -c_sqrt_norm_DMSO, -d_sqrt_norm_DMSO) %>%
  # mutate(pred_norm_min = c_norm_min + ((d_norm_min - c_norm_min) / (1 + exp(b_norm_min*(log(conc) - log(e_norm_min)))))) %>%
  # select(-b_norm_min, -c_norm_min, -d_norm_min) %>%
  # mutate(pred_sqrt_norm_min = c_sqrt_norm_min + ((d_sqrt_norm_min - c_sqrt_norm_min) / (1 + exp(b_sqrt_norm_min*(log(conc) - log(e_sqrt_norm_min)))))) %>%
  # select(-b_sqrt_norm_min, -c_sqrt_norm_min, -d_sqrt_norm_min) %>%
  # select(treatment, strains, conc, genotype, pred_norm_DMSO, pred_sqrt_norm_DMSO, pred_norm_min, pred_sqrt_norm_min) %>%
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  ))

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

  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  ))


ec50_average <- ec50_norm_DMSO %>%
  filter(curve_norm_DMSO == 'N2') %>%
  group_by(curve_norm_DMSO, treatment) %>%
  mutate(sd = sd(estimate_norm_DMSO)) %>%
  ungroup() %>%
  group_by(curve_norm_DMSO, treatment, sd) %>%
  summarise(estimate_avg = mean(estimate_norm_DMSO)) %>%
  mutate(genotype = 'Wild type') %>%
  mutate(route = 'Wild type') 

predict_average <- predict %>%
  filter(strains == 'N2') %>%
  group_by(treatment, strains, conc) %>%
  mutate(average = mean(pred_norm_DMSO, na.rm = TRUE))


(control_well_plot_n2 <- well_summary %>%
    filter(treatment %in% c('DMSO', 'Untreated')) %>%
    filter(strains == "N2") %>%
    replace(is.na(.), 0) %>%
    ggplot(aes(x = genotype, y = mean_raw_length)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 320), breaks = seq(0, 320, 60)) +
    geom_quasirandom(aes(color = plate, shape = drug_stock), alpha = 1, size = 0.5) +
    stat_summary(aes(fill = plate, shape = drug_stock, color = plate), geom = 'point', fun = median, size = 1.5) +
    #stat_summary(aes(fill = plate), geom = 'point', fun = median, size = 2, shape = 24, color = 'black') +
    geom_text(aes(label = ifelse(metadata_date %in% c('0131'),as.character(metadata_date),'')), hjust=10, vjust=0, size = 1) +
    # stat_summary(aes(fill = assay_date),
    #              geom = 'point', fun = mean, size = 3, shape = 22, color = 'black') +
    #scale_x_log10() +
    #facet_grid(rows = vars(strains), cols = vars(treatment), scales = 'free') +
    #rainbow
    scale_color_manual(values = c('paleturquoise4', 'slateblue4', 'darkorange4', 'chocolate', 'darkorchid', 'cornflowerblue', 'chartreuse4', 'rosybrown4', 'black')) + #243da5
    scale_fill_manual(values = c('paleturquoise4', 'slateblue4', 'darkorange4', 'chocolate', 'darkorchid', 'cornflowerblue', 'chartreuse4', 'rosybrown4', 'black')) +
    theme_nw2() +
    theme(legend.position = 'none', axis.ticks = element_blank(), axis.text.x = element_blank(), plot.tag = element_text(size = 10)) +
    guides(fill='none', color='none') +
    labs(tag = "A") +
    ylab('Raw length') +
    xlab("1% DMSO") +
    NULL)

#ggsave(here('/Users/elenagarncarz/Desktop/Fig_1_DMSO_plot.pdf'), control_well_plot_n2, width = 150, height = 200, units = 'mm')
#ggsave(here('/Users/elenagarncarz/Desktop/Fig_1_DMSO_plot.png'), control_well_plot_n2, width = 150, height = 200, units = 'mm')

ann_text_A <- data.frame(lab = "EC[50] *=* 15.7", treatment = factor('AlbendazoleSulfoxide',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))
ann_text_I <- data.frame(lab = "EC[50] *=* 8.55", treatment = factor('Ivermectin',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))
ann_text_L <- data.frame(lab = "EC[50] *=* 6.83", treatment = factor('Levamisole',levels = c("AlbendazoleSulfoxide","Ivermectin","Levamisole")))

(n2_plot <- revised_data %>%
    filter(treatment %in% c('AlbendazoleSulfoxide', 'Levamisole', 'Ivermectin'),
           strains == 'N2') %>% 
    ggplot() +
    geom_line(data = predict %>% filter(strains == 'N2', group == '1'), aes(x = conc, group = genotype, y = pred_norm_DMSO),linewidth = 0.5, color = 'paleturquoise4', alpha =0.7) +
    geom_line(data = predict %>% filter(strains == 'N2', group == '2'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'slateblue4', alpha =0.7) +
    geom_line(data = predict %>% filter(strains == 'N2', group == '3'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'darkorange4', alpha =0.7) +
    geom_line(data = predict %>% filter(strains == 'N2', group == '4'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'chocolate', alpha =0.7) +
    geom_line(data = predict %>% filter(strains == 'N2', group == '5'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'chartreuse4', alpha =0.7) +
    geom_line(data = predict %>% filter(strains == 'N2', group == '6'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'cornflowerblue', alpha =0.7) +
    geom_line(data = predict %>% filter(strains == 'N2', group == '7'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'darkorchid', alpha =0.7) +
    geom_line(data = predict %>% filter(strains == 'N2', group == '8'), aes(x = conc, group = genotype, y = pred_norm_DMSO), linewidth = 0.5, color = 'rosybrown4', alpha =0.7) +
    
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
    geom_vline(data = ec50_norm_DMSO %>% filter(curve_norm_DMSO == 'N2', group == 1), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'paleturquoise4', alpha = 0.5) + 
    geom_vline(data = ec50_norm_DMSO %>% filter(curve_norm_DMSO == 'N2', group == 2), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'slateblue4', alpha = 0.5) + 
    geom_vline(data = ec50_norm_DMSO %>% filter(curve_norm_DMSO == 'N2', group == 3), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'darkorange4', alpha = 0.5) + #AZS only
    geom_vline(data = ec50_norm_DMSO %>% filter(curve_norm_DMSO == 'N2', group == 4), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'chocolate', alpha = 0.5) +
    geom_vline(data = ec50_norm_DMSO %>% filter(curve_norm_DMSO == 'N2', group == 5), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'chartreuse4', alpha = 0.5) +
    geom_vline(data = ec50_norm_DMSO %>% filter(curve_norm_DMSO == 'N2', group == 6), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'cornflowerblue', alpha = 0.5) +
    geom_vline(data = ec50_norm_DMSO %>% filter(curve_norm_DMSO == 'N2', group == 7), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'darkorchid', alpha = 0.5) +
    geom_vline(data = ec50_norm_DMSO %>% filter(curve_norm_DMSO == 'N2', group == 8), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 0.5, linetype = 'dashed', color = 'rosybrown4', alpha = 0.5) +
    geom_vline(data = ec50_norm_DMSO_all %>% filter(curve_norm_DMSO == 'N2'), aes(xintercept = estimate_norm_DMSO, group = treatment), size = 1, linetype = 'dashed', color = 'black', alpha = 1) +
    
    geom_line(data = predict_all %>% filter(strains == 'N2'), aes(x = conc, y = pred_norm_DMSO), color = 'black', size = 1, alpha =1) + 
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
      labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
    annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scale_y_continuous(limits = c(0.2, 1.25), expand = c(0, 0)) +
    # scale_color_manual(values = c('paleturquoise4', 'slateblue4', 'paleturquoise4', 'slateblue4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4')) + 
    # scale_fill_manual(values = c('paleturquoise4', 'slateblue4', 'paleturquoise4', 'slateblue4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4', 'darkorange4', 'chocolate', 'chartreuse4', 'cornflowerblue', 'darkorchid', 'rosybrown4')) +
    scale_shape(guide = 'none') +
    coord_cartesian(clip = "off") +
    
    labs(x = 'Drug concentration (µM)', y = "Normalized length",
         color = 'Strain', fill = 'Strain', tag = "B") + 
    facet_grid(cols = vars(treatment), scales = 'free_x',
               labeller = as_labeller(c(AlbendazoleSulfoxide = 'Albendazole sulfoxide (AZS)',
                                        Ivermectin = 'Ivermectin (IVM)',
                                        Levamisole = 'Levamisole (LEV)'))) +
    theme_nw2() +
    geom_richtext(data = ann_text_A,label = "EC<sub>50</sub> = 15.7 µM", aes(x = 0.6, y = 0.305), size = 1.8) + #value from ec50_norm_DMSO_all, albendazolesulfoxide
    geom_richtext(data = ann_text_I,label = "EC<sub>50</sub> = 8.55 nM", aes(x = 0.001, y = 0.305), size = 1.8) + #value from ec50_norm_DMSO_all, ivermectin
    geom_richtext(data = ann_text_L,label = "EC<sub>50</sub> = 6.83 µM", aes(x = 0.1, y = 0.305), size = 1.8) + #value from ec50_norm_DMSO_all, levamisole
    theme(
      legend.position = 'empty',
      axis.text.x = element_markdown(angle = 0, hjust = 0.5, size = 10),
      axis.text.y = element_markdown(angle = 0, hjust = 0.5, size = 10),
      axis.title.x = element_markdown(face = 'plain', size = 10),
      axis.title.y = element_markdown(face = 'plain', size = 10),
      plot.title = element_text(hjust = 0.5, size = 20),
      plot.tag = element_text(size = 10),
      strip.text = element_text(color = ('black'), size = 30),
      plot.margin = unit(c(0.5,0.5,0.5,2), 'lines') 
      ) +
    NULL)


#ggsave(here('/Users/elenagarncarz/Desktop/Fig_1_curves.pdf'), n2_plot, width = 450, height = 200, units = 'mm')
#ggsave(here('/Users/elenagarncarz/Desktop/Fig_1_curves.png'), n2_plot, width = 450, height = 200, units = 'mm')

(Fig_1_combined <- plot_grid(control_well_plot_n2, n2_plot, align = "h", axis = 'b', rel_widths = c(1, 4.5)))
#print(here())
ggsave(here('plots_EGR/Figure_1/Figure_1_updated.pdf'), Fig_1_combined, width = 240, height = 80, units = 'mm')
ggsave(here('plots_EGR/Figure_1/Figure_1_updated.png'), Fig_1_combined, width = 240, height = 80, units = 'mm')

ec50s_averages <- ec50_norm_DMSO %>%
  group_by(treatment) %>%
  mutate(ec50_avg = mean(estimate_norm_DMSO),
         stdev = sd(estimate_norm_DMSO))
