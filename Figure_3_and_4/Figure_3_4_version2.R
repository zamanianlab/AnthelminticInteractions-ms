# data wrangling/plotting
library(tidyverse)
library(metR)

# other plotting
library(ggbeeswarm)
library(ggtext)
library(cowplot)
library(ZamanianLabThemes)
library(patchwork)

# tables
# library(gt)
# library(gtExtras)

# stats
library(drc)
library(broom)
library(synergyfinder)

# misc
library(conflicted)
library(here)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

# analysis ----------------------------------------------------------------
print(here())
# sqrt transform of poisson distribution
trans_data <- read_rds(here('data_EGR/iso_data_EGR.rds'))%>% 
  #trans_data <- read_rds('/Users/mzamanian/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Nic/project-drug_access/data_EGR/iso_data_EGR.rds') %>% 
  mutate(sqrt_length = sqrt(area_shape_major_axis_length)) %>%
  rename(raw_length = area_shape_major_axis_length) %>%
  select(plate, well, conc1, conc2, strains, treatment1, treatment2, metadata_date, raw_length, sqrt_length)

# DMSO length for normalization
DMSO_mean <- trans_data %>%
  filter(treatment1 %in% c('DMSO', 'Untreated')) %>%
  group_by(metadata_date, strains) %>% 
  summarise(DMSO_mean = mean(raw_length, na.rm = FALSE)) %>%
  ungroup()

# normalize
normalized_data <- trans_data %>% 
  left_join(DMSO_mean) %>%
  mutate(norm_DMSO = raw_length / DMSO_mean)

# sanity check
# normalized_data %>%
#   filter(treatment1 == 'DMSO') %>%
#   ggplot(aes(x = strains, y = raw_length)) +
#   geom_quasirandom(aes(color = metadata_date), alpha = 0.5) +
#   theme_nw2() +
#   NULL
# 
# # should be centered at 1
# normalized_data %>%
#   filter(treatment1 == 'DMSO') %>%
#   ggplot(aes(x = strains, y = norm_DMSO)) +
#   geom_quasirandom(aes(color = metadata_date), alpha = 0.5) +
#   theme_nw2() +
#   NULL

well_summary <- normalized_data %>% 
  mutate(across(contains('conc'), ~ case_when(
    treatment1 == 'DMSO' ~ 0,
    TRUE ~ .
  ))) %>% 
  group_by(metadata_date, plate, strains, well, treatment1, treatment2, conc1, conc2) %>% 
  summarise(mean_raw_length = mean(raw_length),
            mean_norm_DMSO = mean(norm_DMSO)) %>%
  mutate(route = case_when(
    strains == 'N2' ~ 'Wild type',
    strains %in% c('PR672', 'SP1234') ~ 'Amphid',
    strains %in% c('DC19', 'LC144') ~ 'Cuticle',
    strains %in% c('DA453', 'AE501') ~ 'Digestive'
  )) %>%
  mutate(genotype = case_when(
    strains == 'N2' ~ 'Wild type',
    strains == 'AE501' ~ '*nhr-8(ok186)*',
    strains == 'DA453' ~ '*eat-2(ad453)*',
    strains == 'DC19' ~ '*bus-5(br19)*',
    strains == 'LC144' ~ '*agmo-1(e3016)*',
    strains == 'PR672' ~ '*che-1(p672)*',
    strains == 'SP1234' ~ '*dyf-2(m160)*'
  )) 

reps <- well_summary %>% 
  ungroup() %>% 
  filter(treatment1 != 'DMSO') %>% 
  distinct(metadata_date, genotype, strains, treatment1, treatment2) %>% 
  arrange(desc(metadata_date)) %>% 
  group_by(genotype, strains, treatment1, treatment2) %>% 
  mutate(rep = row_number())

# summarize every combination for every strain
all_strains <- well_summary %>% 
  group_by(strains, genotype, treatment1, treatment2, conc1, conc2) %>% 
  summarize(mean_well_raw_length = mean(mean_raw_length),
            mean_well_norm_DMSO_length = mean(mean_norm_DMSO)) %>%
  mutate(
    textcolor = case_when(
      mean_well_raw_length > 0.7 ~ 'black',
      mean_well_norm_DMSO_length > 0.7 ~ 'black',
      TRUE ~ 'white'
    ),
    genotype = factor(genotype, levels = c('Wild type', '*eat-2(ad453)*', '*nhr-8(ok186)*', '*agmo-1(e3016)*', '*bus-5(br19)*', '*che-1(p672)*', '*dyf-2(m160)*'))) 

synergy <- normalized_data %>% 
  mutate(
    inhibition_norm_DMSO = ((1 - norm_DMSO) * 100)) %>% 
  select(plate, strains, well, contains(c('treatment', 'conc')), norm_DMSO, contains(c('inhibition')), raw_length, sqrt_length) %>% 
  mutate(
    conc1 = case_when(
      treatment1 == 'DMSO' ~ 0,
      TRUE ~ conc1),
    conc2 = case_when(
      treatment1 == 'DMSO' ~ 0,
      TRUE ~ conc2)) %>%
  group_by(plate, well, strains, conc1, conc2, treatment1, treatment2) %>%
  summarise(mean_inhibition_norm_DMSO = mean(inhibition_norm_DMSO))

(azs_lev_norm_DMSO <- all_strains %>% 
    filter(
      treatment1 == 'AlbendazoleSulfoxide',
      treatment2 == 'Levamisole',
      conc1 > 0.78124,
      conc1 <  101,
      conc2 > 0.0390624, 
      conc2 < 5.1
    ) %>% 
    ggplot(aes(x = conc1, y = conc2)) +
    metR::geom_contour_fill(aes(z = mean_well_norm_DMSO_length), breaks = MakeBreaks(binwidth = 0.03)) +
    scale_x_log10(
      breaks = c(1, 10, 100),
      labels = c(1, 10, 100),
      expand = c(0, 0)
    ) +
    scale_y_log10(
      breaks = c(0.1, 1),
      labels = c(0.1, 1),
      expand = c(0, 0)
    ) +
    annotation_logticks(sides = 'lb', size = 0.25, outside = TRUE,
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scale_fill_viridis_c(option = 'magma', limits = c(0, 1.5), breaks = c(0, 0.75, 1.5)) +
    coord_fixed(clip = "off") +
    facet_wrap(facets = vars(genotype), nrow = 1) +
    labs(x = 'Albendazole sulfoxide (µM)', y = 'Levamisole (µM)', fill = 'Worm length\n(% control)', tag = "A") +
    #theme_dark() +
    theme(
      legend.position = 'right',
      legend.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x = element_markdown(size = 10),
      strip.background =element_rect(fill="white"),
      axis.text.y = element_text(margin = margin(r = 5))
    ) +
    NULL)

# save_plot(here('Desktop/A_L_raw_norm_DMSO.pdf'),
#           azs_lev_norm_DMSO,
#           base_height = 4, base_width = 14)

# synergy
azs_lev_s_norm_DMSO <- synergy %>%
  ungroup() %>% 
  filter(treatment1 %in% c('DMSO', 'AlbendazoleSulfoxide'),
         treatment2 %in% c(NA, 'Levamisole'),
         conc1 <  200,
         conc2 < 10) %>% 
  # keep the DMSO data but recode the treatments
  mutate(treatment1 = 'AlbendazoleSulfoxide',
         treatment2 = 'Levamisole') %>% 
  group_by(strains) %>% 
  mutate(block_id = cur_group_id()) %>% 
  ungroup() %>% 
  select(block_id, drug1 = treatment1, drug2 = treatment2,
         conc1, conc2, response = mean_inhibition_norm_DMSO) %>%
  mutate(conc_unit1 = 'µM', 
         conc_unit2 = 'µM') 


azs_lev_s_res_norm_DMSO <- ReshapeData(
  data = azs_lev_s_norm_DMSO,
  data_type = "inhibition",
  impute = TRUE,
  impute_method = NULL,
  noise = TRUE,
  iteration = 10,
  seed = 1)


azs_lev_s_res_norm_DMSO <- CalculateSynergy(
  data = azs_lev_s_res_norm_DMSO,
  method = 'ZIP')

azs_lev_s_res_norm_DMSO <- CalculateSensitivity(
  data = azs_lev_s_res_norm_DMSO,
  correct_baseline = "non",
  iteration = 10)

azs_lev_fit_norm_DMSO <- azs_lev_s_res_norm_DMSO$synergy_scores %>% 
  mutate(genotype = case_when(
    block_id == 1 ~ '*nhr-8(ok186)*',
    block_id == 2 ~ '*eat-2(ad453)*',
    block_id == 3 ~ '*bus-5(br19)*',
    block_id == 4 ~ '*agmo-1(e3016)*',
    block_id == 5 ~ 'Wild type',
    block_id == 6 ~ '*che-1(p672)*',
    block_id == 7 ~ '*dyf-2(m160)*',
  )) %>% 
  mutate(genotype = factor(genotype, levels = c('Wild type', '*eat-2(ad453)*', '*nhr-8(ok186)*', '*agmo-1(e3016)*', '*bus-5(br19)*', '*che-1(p672)*', '*dyf-2(m160)*'))) %>%
  filter(
    conc1 > 0.78124,
    conc1 <  100.1,
    conc2 > 0.0390624, 
    conc2 < 5.1
  ) %>% 
  pivot_longer(cols = c(ZIP_fit, ZIP_synergy), names_to = 'ZIP_metric', values_to = 'value')

azs_lev_output <- data.frame(azs_lev_s_res_norm_DMSO$drug_pairs) %>%
  mutate(genotype = case_when(
    block_id == 1 ~ '*nhr-8(ok186)*',
    block_id == 2 ~ '*eat-2(ad453)*',
    block_id == 3 ~ '*bus-5(br19)*',
    block_id == 4 ~ '*agmo-1(e3016)*',
    block_id == 5 ~ 'Wild type',
    block_id == 6 ~ '*che-1(p672)*',
    block_id == 7 ~ '*dyf-2(m160)*',
  )) 

azs_lev_responses <- data.frame(azs_lev_s_res_norm_DMSO$response) %>%
  mutate(genotype = case_when(
    block_id == 1 ~ '*nhr-8(ok186)*',
    block_id == 2 ~ '*eat-2(ad453)*',
    block_id == 3 ~ '*bus-5(br19)*',
    block_id == 4 ~ '*agmo-1(e3016)*',
    block_id == 5 ~ 'Wild type',
    block_id == 6 ~ '*che-1(p672)*',
    block_id == 7 ~ '*dyf-2(m160)*',
  )) %>%
  mutate(treatment = 'azs_lev')

# (azs_lev_fit_plot_norm_DMSO <- azs_lev_fit_norm_DMSO  %>% 
#     filter(ZIP_metric == 'ZIP_fit') %>% 
#     ggplot() +
#     metR::geom_contour_fill(aes(x = conc1, y = conc2, z = value),
#     ) +
#     scale_x_log10(
#       breaks = c(10),
#       labels = scales::trans_format("log10", scales::math_format(10^.x)),
#       expand = c(0, 0)
#     ) +
#     scale_y_log10(
#       breaks = c(0.1, 1),
#       labels = scales::trans_format("log10", scales::math_format(10^.x)),
#       expand = c(0, 0)
#     ) +
#     annotation_logticks(sides = 'lb', size = 0.25, outside = TRUE,
#                         short = unit(0.1, "cm"),
#                         mid = unit(0.15, "cm"),
#                         long = unit(0.2, "cm")) +
#     scale_fill_viridis_c(limits = c(-15, 100), breaks = c( -15, 0, 25, 50, 75, 100)) +
#     coord_fixed(clip = "off") +
#     facet_grid(cols = vars(genotype)) +
#     labs(x = 'Albendazole Sulfoxide (µM)', y = 'Levamisole (µM)', fill = 'ZIP fit\n(% inhibition)') +
#     theme_dark() +
#     theme(
#       legend.position = 'right',
#       legend.title = element_text(hjust = 0.5),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       strip.text.x = element_markdown()
#     ) +
#     NULL)

# save_plot('Desktop/A_L_fit_norm_DMSO.pdf',
#           azs_lev_fit_plot_norm_DMSO ,
#           base_height = 4, base_width = 14)

(azs_lev_syn_plot_norm_DMSO <- azs_lev_fit_norm_DMSO %>% 
    filter(ZIP_metric == 'ZIP_synergy') %>% 
    ggplot() +
    metR::geom_contour_fill(aes(x = conc1, y = conc2, z = value), breaks = MakeBreaks(binwidth = 2)
    ) +
    #stat_subset(aes(x = conc1, y = conc2, subset = value > 10 | value < -10), geom = "point", size = 0.1) +
    geom_contour(aes(x = conc1, y = conc2, z = value), color = "black", size = 0.25, breaks = c(-10, 10)) +
    scale_x_log10(
      breaks = c(1, 10, 100),
      labels = c(1, 10, 100),
      expand = c(0, 0)
    ) +
    scale_y_log10(
      breaks = c(0.1, 1),
      labels = c(0.1, 1),
      expand = c(0, 0)
    ) +
    annotation_logticks(sides = 'lb', size = 0.25, outside = TRUE,
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scico::scale_fill_scico(palette = 'vik', limits = c(-60, 60), breaks = c( -60, -30, -10, 0, 10, 30, 60)) +
    coord_fixed(clip = "off") +
    facet_grid(cols = vars(genotype)) +
    labs(x = 'Albendazole sulfoxide (µM)', y = 'Levamisole (µM)', fill = 'ZIP synergy\n(% expected\nresponse)', tag = "A") +
    theme(
      legend.position = 'right',
      legend.title = element_text(hjust = 0.5),
      legend.text = element_text(angle = 45),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x = element_markdown(size = 10),
      axis.text.x = element_text(angle = , vjust = -1),
      strip.background =element_rect(fill="white"),
      axis.text.y = element_text(hjust = 0.75)
    ) +
    NULL)

# save_plot('Desktop/A_L_synergy_norm_DMSO.pdf',
#           azs_lev_syn_plot_norm_DMSO, 
#           base_height = 4, base_width = 14)

# all AZS x LEV
# save_plot('Desktop/A_L_all_norm_DMSO.pdf',
#           plot_grid(azs_lev_norm_DMSO, azs_lev_fit_plot_norm_DMSO, azs_lev_syn_plot_norm_DMSO,
#                     nrow = 3, align = 'v', axis = 'lr'),
#           base_height = 10, base_width = 14)


# azs vs ivm --------------------------------------------------------------

# raw
(azs_ivm_norm_DMSO <- all_strains %>% 
   filter(
     treatment1 == 'AlbendazoleSulfoxide',
     treatment2 == 'Ivermectin',
     conc1 > 0.78124, 
     conc1 < 100.1,
     conc2 > 0.00015624, 
     conc2 < 0.021
   ) %>% 
   ggplot(aes(x = conc1, y = conc2)) +
   metR::geom_contour_fill(aes(z = mean_well_norm_DMSO_length), breaks = MakeBreaks(binwidth = 0.03)) +
   scale_x_log10(
     breaks = c(10, 100),
     labels = c(10, 100),
     expand = c(0, 0)
   ) +
   scale_y_log10(
     breaks = c(.01, .001),
     labels = c(.01, .001),
     expand = c(0, 0)
   ) +
   annotation_logticks(sides = 'lb', size = 0.25, outside = TRUE,
                       short = unit(0.1, "cm"),
                       mid = unit(0.15, "cm"),
                       long = unit(0.2, "cm")) +
   scale_fill_viridis_c(option = 'magma', limits = c(0, 1.5), breaks = c(0, 0.75, 1.5)) +
   coord_fixed(clip = "off") +
   facet_wrap(facets = vars(genotype), nrow = 1) +
   labs(x = 'Albendazole sulfoxide (µM)', y = 'Ivermectin (µM)', fill = 'Worm length', tag = "B") +
   theme_dark() +
   theme(
     legend.position = 'right',
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     strip.text.x = element_blank(),
     strip.background =element_rect(fill="white"),
     axis.text.y = element_text(margin = margin(r = 5))
   ) +
   NULL)

# save_plot(here('Desktop/A_I_raw_mean_norm_DMSO.pdf'),
#           azs_ivm_norm_DMSO,
#           base_height = 4, base_width = 14)

# synergy
azs_ivm_s_norm_DMSO <- synergy %>%
  ungroup() %>% 
  filter(
    treatment1 %in% c('DMSO', 'AlbendazoleSulfoxide'),
    treatment2 %in% c(NA, 'Ivermectin'),
    conc1 <  200,
    conc2 < 0.04
  ) %>% 
  # keep the DMSO data but recode the treatments
  mutate(treatment1 = 'AlbendazoleSulfoxide',
         treatment2 = 'Ivermectin') %>% 
  group_by(strains) %>% 
  mutate(block_id = cur_group_id()) %>% 
  ungroup() %>% 
  select(block_id, drug1 = treatment1, drug2 = treatment2,
         conc1, conc2, response = mean_inhibition_norm_DMSO) %>%
  mutate(conc_unit1 = 'µM', 
         conc_unit2 = 'µM') 

azs_ivm_s_res_norm_DMSO <- ReshapeData(
  data = azs_ivm_s_norm_DMSO,
  data_type = "inhibition",
  impute = TRUE,
  impute_method = NULL,
  noise = TRUE,
  iteration = 10, # Number of iterations for bootstrapping
  seed = 1)

azs_ivm_s_res_norm_DMSO <- CalculateSynergy(
  data = azs_ivm_s_res_norm_DMSO, 
  method = 'ZIP'
)

azs_ivm_s_res_norm_DMSO <- CalculateSensitivity(
  data = azs_ivm_s_res_norm_DMSO,
  correct_baseline = "non",
  iteration = 10)

azs_ivm_fit_norm_DMSO <- azs_ivm_s_res_norm_DMSO$synergy_scores %>% 
  mutate(genotype = case_when(
    block_id == 1 ~ '*nhr-8(ok186)*',
    block_id == 2 ~ '*eat-2(ad453)*',
    block_id == 3 ~ '*bus-5(br19)*',
    block_id == 4 ~ '*agmo-1(e3016)*',
    block_id == 5 ~ 'Wild type',
    block_id == 6 ~ '*che-1(p672)*',
    block_id == 7 ~ '*dyf-2(m160)*',
  )) %>% 
  mutate(genotype = factor(genotype, levels = c('Wild type', '*eat-2(ad453)*', '*nhr-8(ok186)*', '*agmo-1(e3016)*', '*bus-5(br19)*', '*che-1(p672)*', '*dyf-2(m160)*'))) %>%
  filter(
    conc1 > 0.78124, 
    conc1 < 100.1,
    conc2 > 0.00015624, 
    conc2 < 0.021
  ) %>% 
  pivot_longer(cols = c(ZIP_fit, ZIP_synergy), names_to = 'ZIP_metric', values_to = 'value')

azs_ivm_output <- data.frame(azs_ivm_s_res_norm_DMSO$drug_pairs) %>%
  mutate(genotype = case_when(
    block_id == 1 ~ '*nhr-8(ok186)*',
    block_id == 2 ~ '*eat-2(ad453)*',
    block_id == 3 ~ '*bus-5(br19)*',
    block_id == 4 ~ '*agmo-1(e3016)*',
    block_id == 5 ~ 'Wild type',
    block_id == 6 ~ '*che-1(p672)*',
    block_id == 7 ~ '*dyf-2(m160)*',
  )) 

azs_ivm_responses <- data.frame(azs_ivm_s_res_norm_DMSO$response) %>%
  mutate(genotype = case_when(
    block_id == 1 ~ '*nhr-8(ok186)*',
    block_id == 2 ~ '*eat-2(ad453)*',
    block_id == 3 ~ '*bus-5(br19)*',
    block_id == 4 ~ '*agmo-1(e3016)*',
    block_id == 5 ~ 'Wild type',
    block_id == 6 ~ '*che-1(p672)*',
    block_id == 7 ~ '*dyf-2(m160)*',
  )) %>%
  mutate(treatment = 'azs_ivm')

# (azs_ivm_fit_plot_norm_DMSO <- azs_ivm_fit_norm_DMSO %>% 
#     filter(ZIP_metric == 'ZIP_fit') %>% 
#     ggplot() +
#     metR::geom_contour_fill(aes(x = conc1, y = conc2, z = value),
#     ) +
#     scale_x_log10(
#       breaks = c(10),
#       labels = scales::trans_format("log10", scales::math_format(10^.x)),
#       expand = c(0, 0)
#     ) +
#     scale_y_log10(
#       breaks = c(.01, .001),
#       labels = scales::trans_format("log10", scales::math_format(10^.x)),
#       expand = c(0, 0)
#     ) +
#     annotation_logticks(sides = 'lb', size = 0.25, outside = TRUE,
#                         short = unit(0.1, "cm"),
#                         mid = unit(0.15, "cm"),
#                         long = unit(0.2, "cm")) +
#     scale_fill_viridis_c(limits = c(-15, 100), breaks = c( -15, 0, 25, 50, 75, 100)) +
#     coord_fixed(clip = "off") +
#     facet_grid(cols = vars(genotype)) +
#     labs(x = 'Albendazole Sulfoxide (µM)', y = 'Ivermectin (µM)', fill = 'ZIP fit\n(% inhibition)') +
#     theme_dark() +
#     theme(
#       legend.position = 'right',
#       legend.title = element_text(hjust = 0.5),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       strip.text.x = element_markdown()
#     ) +
#     NULL)

# save_plot(here('Desktop/A_I_fit_norm_DMSO.pdf'),
#           azs_ivm_fit_plot_norm_DMSO,
#           base_height = 4, base_width = 14)

(azs_ivm_syn_plot_norm_DMSO <- azs_ivm_fit_norm_DMSO %>% 
    filter(ZIP_metric == 'ZIP_synergy') %>% 
    ggplot() +
    metR::geom_contour_fill(aes(x = conc1, y = conc2, z = value), breaks = MakeBreaks(binwidth = 2)
    ) +
    #stat_subset(aes(x = conc1, y = conc2, subset = value > 10 | value < -10), geom = "point", size = 0.1) +
    geom_contour(aes(x = conc1, y = conc2, z = value), color = "black", size = 0.25, breaks = c(-10, 10)) +
    scale_x_log10(
      breaks = c(1, 10, 100),
      labels = c(1, 10, 100),
      expand = c(0, 0)
    ) +
    scale_y_log10(
      breaks = c(0.001, 0.01),
      labels = c(0.001, 0.01),
      expand = c(0, 0)
    ) +
    annotation_logticks(sides = 'lb', size = 0.25, outside = TRUE,
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scico::scale_fill_scico(palette = 'vik', limits = c(-60, 60), breaks = c(-60, -10, 10, 60)) +
    coord_fixed(clip = "off") +
    facet_grid(cols = vars(genotype)) +
    labs(x = 'Albendazole sulfoxide (µM)', y = 'Ivermectin (µM)', fill = 'ZIP synergy\n(% expected\nresponse)', tag = "B") +
    theme_dark() +
    theme(
      legend.position = 'right',
      legend.title = element_text(hjust = 0.5),
      legend.text = element_text(angle = 45),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 0, vjust = -1),
      strip.background =element_rect(fill="white"),
      axis.text.y = element_text(hjust = 0.75)
    ) +
    NULL)

# save_plot(here('Desktop/A_I_synergy_norm_DMSO.pdf'),
#           azs_ivm_syn_plot_norm_DMSO,
#           base_height = 4, base_width = 14)

# all AZS x IVM
# save_plot(here('Desktop/A_I_all_norm_DMSO.pdf'),
#           plot_grid(azs_ivm_norm_DMSO, azs_ivm_fit_plot_norm_DMSO, azs_ivm_syn_plot_norm_DMSO,
#                     nrow = 3, align = 'v', axis = 'lr'),
#           base_height = 10, base_width = 14)

# ivm vs lev --------------------------------------------------------------

# raw
(ivm_lev_norm_DMSO <- all_strains %>% 
   filter(
     treatment1 == 'Ivermectin',
     treatment2 == 'Levamisole',
     conc1 > 0.00015624, 
     conc1 < 0.021,
     conc2 > 0.0390624, 
     conc2 < 5.1
   ) %>% 
   ggplot(aes(x = conc2, y = conc1)) +
   metR::geom_contour_fill(aes(z = mean_well_norm_DMSO_length), breaks = MakeBreaks(binwidth = 0.03)) +
   scale_x_log10(
     breaks = c(1, .1),
     labels = c(1, .1),
     expand = c(0, 0)
   ) +
   scale_y_log10(
     breaks = c(0.001, 0.01),
     labels = c(0.001, 0.01),
     expand = c(0, 0)
   ) +
   annotation_logticks(sides = 'lb', size = 0.25, outside = TRUE,
                       short = unit(0.1, "cm"),
                       mid = unit(0.15, "cm"),
                       long = unit(0.2, "cm")) +
   scale_fill_viridis_c(option = 'magma', limits = c(0, 1.5), breaks = c(0, 0.75, 1.5)) +
   coord_fixed(clip = "off") +
   facet_wrap(facets = vars(genotype), nrow = 1) +
   labs(x = 'Levamisole (µM)', y = 'Ivermectin (µM)', fill = 'Worm length\n(% control)', tag = "C") +
   #theme_dark() +
   theme(
     legend.position = 'right',
     legend.title = element_text(hjust = 0.5),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     strip.text.x = element_blank(),
     strip.background =element_rect(fill="white"),
     axis.text.y = element_text(margin = margin(r = 5))
   ) +
   NULL)

# save_plot(here('Desktop/I_L_raw_norm_DMSO.pdf'),
#           ivm_lev_norm_DMSO,
#           base_height = 4, base_width = 14)

# synergy
ivm_lev_s_norm_DMSO <- synergy %>%
  ungroup() %>% 
  filter(
    treatment1 %in% c('DMSO', 'Ivermectin'),
    treatment2 %in% c(NA, 'Levamisole'),
    conc1 < 0.04, 
    conc2 < 10) %>% 
  # keep the DMSO data but recode the treatments
  mutate(treatment1 = 'Ivermectin',
         treatment2 = 'Levamisole') %>% 
  group_by(strains) %>% 
  mutate(block_id = cur_group_id()) %>% 
  ungroup() %>% 
  select(block_id, drug1 = treatment1, drug2 = treatment2,
         conc1, conc2, response = mean_inhibition_norm_DMSO) %>%
  mutate(conc_unit1 = 'µM', 
         conc_unit2 = 'µM')

ivm_lev_s_res_norm_DMSO <- ReshapeData(
  data = ivm_lev_s_norm_DMSO,
  data_type = "inhibition",
  impute = TRUE,
  impute_method = NULL,
  noise = TRUE,
  iteration = 10, 
  seed = 1)

ivm_lev_s_res_norm_DMSO <- CalculateSynergy(
  data = ivm_lev_s_res_norm_DMSO, 
  method = 'ZIP')

ivm_lev_s_res_norm_DMSO <- CalculateSensitivity(
  data = ivm_lev_s_res_norm_DMSO,
  correct_baseline = "non",
  iteration = 10)

ivm_lev_fit_norm_DMSO <- ivm_lev_s_res_norm_DMSO$synergy_scores %>% 
  mutate(genotype = case_when(
    block_id == 1 ~ '*nhr-8(ok186)*',
    block_id == 2 ~ '*eat-2(ad453)*',
    block_id == 3 ~ '*bus-5(br19)*',
    block_id == 4 ~ '*agmo-1(e3016)*',
    block_id == 5 ~ 'Wild type',
    block_id == 6 ~ '*che-1(p672)*',
    block_id == 7 ~ '*dyf-2(m160)*',
  )) %>% 
  mutate(genotype = factor(genotype, levels = c('Wild type', '*eat-2(ad453)*', '*nhr-8(ok186)*', '*agmo-1(e3016)*', '*bus-5(br19)*', '*che-1(p672)*', '*dyf-2(m160)*'))) %>%
  filter(
    conc1 > 0.00015624, 
    conc1 < 0.021,
    conc2 > 0.0390624, 
    conc2 < 5.1
  ) %>% 
  pivot_longer(cols = c(ZIP_fit, ZIP_synergy), names_to = 'ZIP_metric', values_to = 'value')

ivm_lev_output <- data.frame(ivm_lev_s_res_norm_DMSO$drug_pairs) %>%
  mutate(genotype = case_when(
    block_id == 1 ~ '*nhr-8(ok186)*',
    block_id == 2 ~ '*eat-2(ad453)*',
    block_id == 3 ~ '*bus-5(br19)*',
    block_id == 4 ~ '*agmo-1(e3016)*',
    block_id == 5 ~ 'Wild type',
    block_id == 6 ~ '*che-1(p672)*',
    block_id == 7 ~ '*dyf-2(m160)*',
  )) 

ivm_lev_responses <- data.frame(ivm_lev_s_res_norm_DMSO$response) %>%
  mutate(genotype = case_when(
    block_id == 1 ~ '*nhr-8(ok186)*',
    block_id == 2 ~ '*eat-2(ad453)*',
    block_id == 3 ~ '*bus-5(br19)*',
    block_id == 4 ~ '*agmo-1(e3016)*',
    block_id == 5 ~ 'Wild type',
    block_id == 6 ~ '*che-1(p672)*',
    block_id == 7 ~ '*dyf-2(m160)*',
  )) %>%
  mutate(treatment = 'ivm_lev')

  
# (ivm_lev_fit_plot_norm_DMSO <- ivm_lev_fit_norm_DMSO %>% 
#     filter(ZIP_metric == 'ZIP_fit') %>% 
#     ggplot() +
#     metR::geom_contour_fill(aes(x = conc1, y = conc2, z = value),
#     ) +
#     scale_x_log10(
#       breaks = c(0.001, 0.01),
#       labels = scales::trans_format("log10", scales::math_format(10^.x)),
#       expand = c(0, 0)
#     ) +
#     scale_y_log10(
#       breaks = c(1, .1),
#       labels = scales::trans_format("log10", scales::math_format(10^.x)),
#       expand = c(0, 0)
#     ) +
#     annotation_logticks(sides = 'lb', size = 0.25, outside = TRUE,
#                         short = unit(0.1, "cm"),
#                         mid = unit(0.15, "cm"),
#                         long = unit(0.2, "cm")) +
#     scale_fill_viridis_c(limits = c(-15, 100), breaks = c( -15, 0, 25, 50, 75, 100)) +
#     coord_fixed(clip = "off") +
#     facet_grid(cols = vars(genotype)) +
#     labs(x = 'Ivermectin (µM)', y = 'Levamisole (µM)', fill = 'ZIP fit\n(% inhibition)') +
#     theme_dark() +
#     theme(
#       legend.position = 'right',
#       legend.title = element_text(hjust = 0.5),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       strip.text.x = element_markdown()
#     ) +
#     NULL)

# save_plot(here('Desktop/I_L_fit_norm_DMSO.pdf'),
#           ivm_lev_fit_plot_norm_DMSO,
#           base_height = 4, base_width = 14)


(ivm_lev_syn_plot_norm_DMSO <- ivm_lev_fit_norm_DMSO %>% 
    filter(ZIP_metric == 'ZIP_synergy') %>% 
    ggplot() +
    metR::geom_contour_fill(aes(x = conc2, y = conc1, z = value), breaks = MakeBreaks(binwidth = 2)
    ) +
    geom_contour(aes(x = conc2, y = conc1, z = value), color = "black", size = 0.25, breaks = c(-10, 10)) +
    #geom_text_contour(aes(x = conc1, y = conc2, z = value), stroke = 0.1) +
    #geom_contour(aes(x = conc1, y = conc2, z = value), color = "white", size = 0.25, breaks = c( 0)) +
    #stat_subset(aes(x = conc1, y = conc2, subset = value > 10 | value < -10), geom = "point", size = 0.7) +
    scale_x_log10(
      breaks = c(0.1, 1),
      labels = c(0.1, 1),
      expand = c(0, 0)
    ) +
    scale_y_log10(
      breaks = c(0.001, 0.01),
      labels = c(0.001, 0.01),
      expand = c(0, 0)
    ) +
    annotation_logticks(sides = 'lb', size = 0.25, outside = TRUE,
                        short = unit(0.1, "cm"),
                        mid = unit(0.15, "cm"),
                        long = unit(0.2, "cm")) +
    scico::scale_fill_scico(palette = 'vik', limits = c(-60, 60), breaks = c(-60, -10, 10, 60)) +
    coord_fixed(clip = "off") +
    facet_grid(cols = vars(genotype)) +
    labs(x = 'Levamisole (µM)', y = 'Ivermectin (µM)', fill = 'ZIP synergy\n(% expected\nresponse)', tag = "C") +
    #theme_dark() +
    theme(
      legend.position = 'right',
      legend.title = element_text(hjust = 0.5),
      legend.text = element_text(angle = 20),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x = element_blank(),
      axis.text.x = element_text(angle = 0, vjust = -1),
      strip.background =element_rect(fill="white"),
      axis.text.y = element_text(hjust = 0.75)
      #panel.border = element_rect(colour = "black", fill=NA, size=0.4)
    ) +
    NULL)

# save_plot(here('Desktop/I_L_synergy_norm_DMSO.pdf'),
#           ivm_lev_syn_plot_norm_DMSO,
#           base_height = 4, base_width = 14)


# all IVM x LEV
# save_plot(here('Desktop/I_L_all_norm_DMSO.pdf'),
#           plot_grid(ivm_lev_norm_DMSO, ivm_lev_fit_plot_norm_DMSO, ivm_lev_syn_plot_norm_DMSO,
#                     nrow = 3, align = 'v', axis = 'lr'),
#           base_height = 10, base_width = 14)




# all_outputs <- rbind(azs_ivm_output, azs_lev_output, ivm_lev_output) %>%
#   select(-block_id, -conc_unit1, -conc_unit2) %>%
#   select(drug1, drug2, response_p_value, response_origin_p_value, ZIP_synergy_p_value, genotype)
#   # mutate(significance = case_when(
#   #   ZIP_synergy_p_value > 0.05 ~ as.character(ZIP_synergy_p_value),
#   #   ZIP_synergy_p_value < 0.0001 ~ paste(as.character(ZIP_synergy_p_value), '****', sep=""),
#   #   ZIP_synergy_p_value < 0.001 ~ paste(as.character(ZIP_synergy_p_value), '***', sep=""),
#   #   ZIP_synergy_p_value < 0.01 ~ paste(as.character(ZIP_synergy_p_value), '**', sep=""),
#   #   ZIP_synergy_p_value < 0.05 ~ paste(as.character(ZIP_synergy_p_value), '*', sep=""),
#   # ))
# 
# all_responses <- rbind(azs_ivm_responses, azs_lev_responses, ivm_lev_responses) %>%
#   group_by(genotype, conc1, conc2, treatment) %>%
#   mutate(response_per_well = mean(response),
#          conc = paste(conc1, "_", conc2)) %>%
#   ungroup() %>%
#   select(genotype, treatment, response_per_well, conc1, conc2) %>%
#   unique()
# 
# (heatmap_ivm_lev <- Plot2DrugHeatmap(
#   data = ivm_lev_s_res_norm_DMSO,
#   plot_block = 1,
#   drugs = c(1, 2),
#   plot_value = "ZIP_synergy",
#   statistic = "ci",
#   dynamic = FALSE,
#   summary_statistic = c("quantile_25", "quantile_75")
# ))
# 
# (heatmap_azs_ivm <- Plot2DrugHeatmap(
#   data = azs_ivm_s_res_norm_DMSO,
#   plot_block = 1,
#   drugs = c(1, 2),
#   plot_value = "ZIP_synergy",
#   statistic = "ci",
#   dynamic = FALSE,
#   summary_statistic = c("quantile_25", "quantile_75")
# ))
# 
# (heatmap_azs_lev <- Plot2DrugHeatmap(
#   data = azs_lev_s_res_norm_DMSO,
#   plot_block = 1,
#   drugs = c(1, 2),
#   plot_value = "ZIP_synergy",
#   statistic = "ci",
#   dynamic = FALSE,
#   summary_statistic = c("quantile_25", "quantile_75")
# ))
# 
# #looking at variance to explain p values for the interactions
# variances <- trans_data %>%
#   group_by(strains, treatment1, treatment2, conc1, conc2) %>%
#   mutate(variance_raw = var(raw_length)) %>%
#   ungroup() %>%
#   group_by(treatment1, treatment2, strains) %>%
#   mutate(variance_avg = mean(variance_raw)) %>%
#   ungroup() %>%
#   select(strains, treatment1, treatment2, variance_avg) %>%
#   unique() %>%
#   na.omit()

save_plot(here('plots_EGR/Figure_3/Figure_3_updated.pdf'),
          plot_grid(azs_lev_norm_DMSO + remove_legend(),
                    azs_ivm_norm_DMSO + remove_legend(),
                    ivm_lev_norm_DMSO + theme(legend.position = 'bottom'),
                    nrow = 3, align = 'h', axis = 'l', rel_heights = c(1, 1, 1.2)), base_height = 10, base_width = 14)

save_plot(here('plots_EGR/Figure_3/Figure_3_updated.png'),
          plot_grid(azs_lev_norm_DMSO + remove_legend(),
                    azs_ivm_norm_DMSO + remove_legend(),
                    ivm_lev_norm_DMSO + theme(legend.position = 'bottom'),
                    nrow = 3, align = 'h', axis = 'l', rel_heights = c(1, 1, 1.2)), base_height = 10, base_width = 14)

save_plot(here('plots_EGR/Figure_4/Figure_4_updated.pdf'),
          plot_grid(azs_lev_syn_plot_norm_DMSO + remove_legend(),
                    azs_ivm_syn_plot_norm_DMSO + remove_legend(),
                    ivm_lev_syn_plot_norm_DMSO + theme(legend.position = 'bottom'),
                    nrow = 3, align = 'h', axis = 'l', rel_heights = c(1, 1, 1.2)), base_height = 10, base_width = 14)

save_plot(here('plots_EGR/Figure_4/Figure_4_updated.png'),
          plot_grid(azs_lev_syn_plot_norm_DMSO + remove_legend(),
                    azs_ivm_syn_plot_norm_DMSO + remove_legend(),
                    ivm_lev_syn_plot_norm_DMSO + theme(legend.position = 'bottom'),
                    nrow = 3, align = 'h', axis = 'l', rel_heights = c(1, 1, 1.2)), base_height = 10, base_width = 14)
# combine plots -----------------------------------------------------------
save_plot(here('plots_EGR/Figure_3/Figure_3_updated.pdf'),
          plot_grid(azs_lev_norm_DMSO + remove_legend(),
                    azs_ivm_norm_DMSO + remove_legend(),
                    ivm_lev_norm_DMSO + theme(legend.position = 'bottom'),
                    nrow = 3, align = 'h', axis = 'l', rel_heights = c(1, 1, 1.2)), base_height = 10, base_width = 14)

save_plot(here('plots_EGR/Figure_3/Figure_3_updated.png'),
          plot_grid(azs_lev_norm_DMSO + remove_legend(),
                    azs_ivm_norm_DMSO + remove_legend(),
                    ivm_lev_norm_DMSO + theme(legend.position = 'bottom'),
                    nrow = 3, align = 'h', axis = 'l', rel_heights = c(1, 1, 1.2)), base_height = 10, base_width = 14)

save_plot(here('Desktop/all_fit_mean_norm_DMSO_updated.pdf'),
          plot_grid(azs_lev_fit_plot_norm_DMSO + remove_legend(),
                    azs_ivm_fit_plot_norm_DMSO + remove_legend(),
                    ivm_lev_fit_plot_norm_DMSO + theme(legend.position = 'bottom'),
                    nrow = 3, align = 'h', axis = 'l', rel_heights = c(1, 1, 1.2)), base_height = 10, base_width = 14)

save_plot(here('plots_EGR/Figure_4/Figure_4_updated.pdf'),
          plot_grid(azs_lev_syn_plot_norm_DMSO + remove_legend(),
                    azs_ivm_syn_plot_norm_DMSO + remove_legend(),
                    ivm_lev_syn_plot_norm_DMSO + theme(legend.position = 'bottom'),
                    nrow = 3, align = 'h', axis = 'l', rel_heights = c(1, 1, 1.2)), base_height = 10, base_width = 14)

save_plot(here('plots_EGR/Figure_4/Figure_4_updated.png'),
          plot_grid(azs_lev_syn_plot_norm_DMSO + remove_legend(),
                    azs_ivm_syn_plot_norm_DMSO + remove_legend(),
                    ivm_lev_syn_plot_norm_DMSO + theme(legend.position = 'bottom'),
                    nrow = 3, align = 'h', axis = 'l', rel_heights = c(1, 1, 1.2)), base_height = 10, base_width = 14)
