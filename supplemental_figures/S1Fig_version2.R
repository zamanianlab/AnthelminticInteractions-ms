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


dc19_data <- revised_data %>% filter(strains == 'DC19', !treatment %in% c('DMSO', 'Untreated'))
# dc19_predict <- predict %>% filter(strains == 'DC19')
dc19_predict_1 <- predict %>% filter(strains == 'DC19', group == '1', !treatment %in% c('DMSO', 'Untreated'))
dc19_predict_4 <- predict %>% filter(strains == 'DC19', group == '4', !treatment %in% c('DMSO', 'Untreated'))
dc19_ec50 <- ec50_norm_DMSO %>% filter(curve_norm_DMSO == 'DC19', !treatment %in% c('DMSO', 'Untreated'))
dc19_well_summary <- well_summary %>% filter(strains == 'DC19', !treatment %in% c('DMSO', 'Untreated'))

#0331 is group 4 and glass tube and blue
#0131 is group 1 and plate and black
(dc19_plot <- dc19_data  %>%
  
  ggplot() +
  
  geom_line(data = dc19_predict_1, aes(x = conc, group = genotype, y = pred_norm_DMSO), size = 1, color = 'black') +
  geom_line(data = dc19_predict_4, aes(x = conc, group = genotype, y = pred_norm_DMSO), size = 1, color = 'cornflowerblue') +
  
  geom_quasirandom(data = dc19_well_summary %>% filter(plate == '20220131-p01-EJG_1133'), aes(x = conc, y = mean_norm_DMSO, color = plate), size = 0.5, width = 0.1, alpha = 0.75, color = 'black') +
  geom_quasirandom(data = dc19_well_summary %>% filter(plate == '20220331-p05-EJG_1315'), aes(x = conc, y = mean_norm_DMSO, color = plate), size = 0.5, width = 0.1, alpha = 0.75, color = 'cornflowerblue') +
  
  geom_vline(data = dc19_ec50 %>% filter(group == '1'), aes(xintercept = estimate_norm_DMSO), size = 0.75, linetype = 'dashed', color = 'black') +
  geom_vline(data = dc19_ec50 %>% filter(group == '4'), aes(xintercept = estimate_norm_DMSO), size = 0.75, linetype = 'dashed', color = 'cornflowerblue') +
  
  scale_x_log10(
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
    labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
  
  annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                      short = unit(0.1, "cm"),
                      mid = unit(0.15, "cm"),
                      long = unit(0.2, "cm")) +
  
  scale_y_continuous(limits = c(0.2, 1.37), expand = c(0, 0), breaks = seq(0.4, 1.2, 0.2)) +


  scale_shape(guide = 'none') +
  
  coord_cartesian(clip = "off") +
  
  labs(x = 'Drug concentration (µM)', y = "Normalized length",
       color = 'Strain', fill = 'Strain', shape = 'Replicate') +

  facet_grid(cols = vars(treatment), scales = 'free_x',
             labeller = as_labeller(c(AlbendazoleSulfoxide = 'AZS',
                                      Ivermectin = 'IVM',
                                      Levamisole = 'LEV'))) +
  theme_nw2() +
  
  theme(
    legend.position = 'empty',
    axis.text.x = element_markdown(angle = 0, hjust = 0.5),
    axis.title.x = element_markdown(face = 'plain'),
    axis.title.y = element_markdown(face = 'plain'),
    plot.background = element_rect(fill = 'white')) +
  # remove_legend() +
  NULL)

ggsave(here('/Users/elenagarncarz/Desktop/Supp_2_bus5_updated.tiff'), dc19_plot, width = 200, height = 100, units = 'mm')
ggsave(here('/Users/elenagarncarz/Desktop/Supp_2_bus5_updated.png'), dc19_plot, width = 200, height = 80, units = 'mm')

