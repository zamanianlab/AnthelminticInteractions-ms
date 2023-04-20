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

# data import -------------------------------------------------------------
print(here())
data_files <- tibble(path = list.files(path = here('data_EGR'), pattern = "*_tidy.csv", recursive = TRUE)) %>%
  filter(!str_detect(path, 'pipe')) %>%
  mutate(path = str_c('data_EGR', path, sep = '/'))

get_data <- function(...) {

  df <- tibble(...)

  data <- read_csv(here(df$path), col_types = c('ccccccccccddccccccccccdddddddddddddddddddddddd'))

}

data <- data_files %>%
  pmap_dfr(get_data) %>%
  janitor::clean_names() %>%
  # fix metadata
  mutate(
    treatment = case_when(
      is.na(treatment) ~ 'Untreated',
      conc == '1p' ~ 'DMSO',
      TRUE ~ treatment),
    genotype = case_when(
      strains == 'AE501' ~ '*nhr-8(ok186)*',
      strains == 'LC144' ~ '*agmo-1(e3016)*',
      strains == 'PR672' ~ '*che-1(p672)*',
      strains == 'SP1234' ~ '*dyf-2(m160)*',
      strains == 'DC19' ~ '*bus-5(br19)*',
      TRUE ~ strains))
  

# prune -------------------------------------------------------------------

# see roundup.csv for the list of good/bad plates
 iso_data <- data %>%
  filter(case_when(
    metadata_date == '20211115'  ~ FALSE,
    metadata_date == '20211118'  ~ FALSE,
    metadata_date == '20211118'  ~ FALSE,
    metadata_date == '20211208'  ~ FALSE,
    metadata_date == '20211216'  ~ FALSE,
    metadata_date == '20220203' & metadata_plate %in% c('p05', 'p06', 'p07', 'p08', 'p09', 'p10', 'p11', 'p12') ~ TRUE,
    metadata_date == '20220210' & metadata_plate %in% c('p02', 'p03', 'p04', 'p05', 'p06', 'p07', 
                                                        'p09', 'p10', 'p11', 'p12', 'p13', 'p14',
                                                        'p16', 'p18', 'p19', 'p20', 'p21', 'p24', 
                                                        'p25', 'p26', 'p27') ~ TRUE,
    metadata_date == '20220311' & metadata_plate %in% c('p02', 'p03') ~ TRUE,
    metadata_date == '20220324' & metadata_plate %in% c('p01', 'p02', 'p03', 'p04', 'p06', 'p07', 'p08', 'p09') ~ TRUE,
    metadata_date == '20220325' & metadata_plate %in% c('p02', 'p03', 'p04', 'p05') ~ TRUE,
    metadata_date == '20220331' & metadata_plate %in% c('p03', 'p04', 'p06', 'p07') ~ TRUE,
    metadata_date == '20220407' & metadata_plate %in% c('p01', 'p02', 'p03', 'p04') ~ TRUE,
    metadata_date == '20220422' & metadata_plate %in% c('p01', 'p02', 'p03', 'p04', 'p05', 'p06') ~ TRUE,
  )) %>%
  mutate(treatment = ifelse(treatment == 'AlbendazoleSulfoxide', 'AlbendazoleSulfoxide_Levamisole', treatment))%>%
  mutate(treatment = ifelse(treatment %in% c('Ivermectin', 'Levamisole'), 'Ivermectin_Levamisole', treatment))%>%
  separate(treatment, c('treatment1', 'treatment2'), sep = '_') %>%
  separate(conc, c('conc1', 'conc2'), sep = '_') %>%
  mutate(
    assay_date = str_extract(plate, '202[1-2][0-9]{4}'),
    conc1 = case_when(
      treatment1 == 'DMSO' ~ 0.01,
      treatment1 != 'DMSO' ~ as.numeric(str_remove(conc1, 'uM'))
    ),
    conc2 = case_when(
      treatment2 == 'DMSO' ~ 0.01,
      treatment2 != 'DMSO' ~ as.numeric(str_remove(conc2, 'uM'))
    )
  ) %>%
  mutate(treatment2 = ifelse(is.na(treatment2), treatment1,treatment2))
  
 # iso_invest <- iso_data %>%
 #  select(plate, well, treatment, conc) %>%
 #   filter(treatment %in% c('AlbendazoleSulfoxide_Levamisole')) %>%
 #  unique() %>%
 #   count(conc)
 # 
 # iso_invest <- iso_data %>%
 #   select(plate, treatment) %>%
 #   filter(treatment == 'Ivermectin') %>%
 #   unique() %>%
 #   count(conc)

### trim outliers
coef <- 1.5

outliers <- iso_data %>%
  group_by(treatment1, treatment2, conc1, conc2, strains) %>%
  group_nest() %>%
  mutate(quant = map(data, ~ as.numeric(quantile(.x$area_shape_major_axis_length, c(0, 0.25, 0.5, 0.75, 1)))),
         iqr = unlist(map(quant, ~ diff(.x[c(2, 4)]))),
         top_cutoff = unlist(map(quant, ~ pluck(.x, 4))) + coef * iqr,
         bottom_cutoff = unlist(map(quant, ~ pluck(.x, 2))) - coef * iqr) %>% 
  select(treatment1, treatment2, conc1, conc2, strains, contains('cutoff'))

trimmed_data <- iso_data %>% 
  left_join(outliers) %>% 
  mutate(outlier = case_when(
    area_shape_major_axis_length > top_cutoff | area_shape_major_axis_length < bottom_cutoff ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  filter(outlier == FALSE) %>% 
  select(-contains('cutoff'), -outlier)

print(here())
#write_rds(trimmed_data, here("iso_EGR_data.rds"))
