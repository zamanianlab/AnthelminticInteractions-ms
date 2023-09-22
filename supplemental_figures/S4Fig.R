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
library(patchwork)


# tables
library(gt)
library(gtExtras)

# stats
library(drc)
library(broom)

# misc
library(conflicted)
library(here)
print(here())
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")

#load iso data
trans_data <- read_rds('/Users/elenagr/Desktop/iso_data_EGR.rds')%>% #trans_data <- read_rds('/Users/elenagarncarz/Library/CloudStorage/Box-Box/ZamanianLab/LabMembers/Nic/project-drug_access/data_EGR/iso_data_EGR.rds')%>%
  rename(raw_length = area_shape_major_axis_length) %>%
  select(plate, well, conc1, conc2, strains, treatment1, treatment2, metadata_date, raw_length)

# find DMSO avg length for normalization
DMSO_mean <- trans_data %>%
  filter(treatment1 %in% c('DMSO', 'Untreated')) %>%
  group_by(metadata_date, strains) %>% 
  summarise(DMSO_mean = mean(raw_length, na.rm = FALSE)) %>%
  ungroup()

# normalize by dividing by DMSO treated avg
normalized_data <- trans_data %>% 
  left_join(DMSO_mean) %>%
  mutate(norm_DMSO = raw_length / DMSO_mean)

# filter for data we want - N2 and dyf-2 strains, IVM_AZS plates only, and at conc where we see a major change (IVM - 0.01, AZS - 50) and (IVM - 0.005, AZS - 3.125)
# antagonism in dyf-2, no significant ant or syn in N2
data_dyf_antshift1 <- normalized_data %>% 
  filter(strains %in% c('SP1234', 'N2'), treatment1 %in% c('DMSO', 'AlbendazoleSulfoxide'), 
         treatment2 %in% c('DMSO', 'Ivermectin'), conc1 %in% c(0, 0.01, 50), conc2 %in% c(NA, 0, 0.01)) %>%
  mutate(treatment = case_when(
    treatment1 == 'DMSO' ~ 'DMSO',
    treatment1 == 'AlbendazoleSulfoxide' & treatment2 == 'Ivermectin' & conc1 != 0 & conc2 != 0 ~ 'AZS_IVM',
    treatment1 == 'AlbendazoleSulfoxide' & treatment2 == 'Ivermectin' & conc1 == 0  ~ 'IVM',
    treatment1 == 'AlbendazoleSulfoxide' & treatment2 == 'Ivermectin' & conc2 == 0  ~ 'AZS',
  )) %>%
  mutate(genotype = case_when(
    strains == 'SP1234' ~ '*dyf-2(m160)*',
    strains == 'N2' ~ 'Wild type'
  )) %>% mutate(
    plate_well = paste(plate, '_', well)
  )

# antagonism in N2, no significant ant or syn in dyf-2
data_dyf_antshift2 <- normalized_data %>% 
  filter(strains %in% c('SP1234', 'N2'), treatment1 %in% c('DMSO', 'AlbendazoleSulfoxide'), 
         treatment2 %in% c('DMSO', 'Ivermectin'), conc1 %in% c(0, 0.01, 3.125), conc2 %in% c(NA, 0, 0.005)) %>%
  mutate(treatment = case_when(
    treatment1 == 'DMSO' ~ 'DMSO',
    treatment1 == 'AlbendazoleSulfoxide' & treatment2 == 'Ivermectin' & conc1 != 0 & conc2 != 0 ~ 'AZS_IVM',
    treatment1 == 'AlbendazoleSulfoxide' & treatment2 == 'Ivermectin' & conc1 == 0  ~ 'IVM',
    treatment1 == 'AlbendazoleSulfoxide' & treatment2 == 'Ivermectin' & conc2 == 0  ~ 'AZS',
  )) %>%
  mutate(genotype = case_when(
    strains == 'SP1234' ~ '*dyf-2(m160)*',
    strains == 'N2' ~ 'Wild type'
  )) %>% mutate(
    plate_well = paste(plate, '_', well)
  )

means_raw1 <- data_dyf_antshift1 %>%
  group_by(treatment, strains) %>%
  summarise( mean_raw = mean(raw_length)) 

means_norm1 <- data_dyf_antshift1 %>%
  group_by(treatment, strains) %>%
  summarise( mean_norm = mean(norm_DMSO))

means_raw2 <- data_dyf_antshift2 %>%
  group_by(treatment, strains) %>%
  summarise( mean_raw = mean(raw_length)) 

means_norm2 <- data_dyf_antshift2 %>%
  group_by(treatment, strains) %>%
  summarise( mean_norm = mean(norm_DMSO))

ann_text <- data.frame(lab = "X")

# make heatmaps from Fig 4, have to run Fig 4 script first before running this method
(heatmap_WT <- azs_ivm_fit_norm_DMSO %>% 
    filter(ZIP_metric == 'ZIP_synergy', genotype %in% c('Wild type')) %>% 
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

    labs(x = 'Albendazole sulfoxide (µM)', y = 'Ivermectin (µM)', fill = 'ZIP synergy\n(% expected\nresponse)', tag = "A") +
    theme_dark() +
    geom_richtext(data = ann_text,label = "C", aes(x = 3.125, y = 0.005), fill = NA, label.color = NA, size = 4) +
    geom_richtext(data = ann_text,label = "D", aes(x = 50, y = 0.0105), fill = NA, label.color = NA, size = 4) +
    ggtitle('Wild type') +
    theme(
      plot.title = element_text(face="italic", hjust = 0.5),
      legend.position = 'none',
      legend.title = element_text(hjust = 0.5),
      legend.text = element_text(angle = 45),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x = element_markdown(size = 10, color = 'black'),
      axis.text.x = element_text(angle = 0, vjust = -1),
      strip.background =element_rect(fill="white"),
      axis.text.y = element_text(hjust = 0.75)
    ) +
    NULL)

(heatmap_dyf <- azs_ivm_fit_norm_DMSO %>% 
    filter(ZIP_metric == 'ZIP_synergy', genotype %in% c('*dyf-2(m160)*')) %>% 
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
    labs(x = 'Albendazole sulfoxide (µM)', y = 'Ivermectin (µM)', fill = 'ZIP synergy\n(% expected\nresponse)', tag = "B") +
    theme_dark() +
    geom_richtext(data = ann_text,label = "C", aes(x = 3.125, y = 0.005), fill = NA, label.color = NA, size = 4) +
    geom_richtext(data = ann_text,label = "D", aes(x = 50, y = 0.0105), fill = NA, label.color = NA, size = 4) +
    ggtitle('dyf-2(m160)') +
    theme(
      plot.title = element_text(face="italic", hjust = 0.5),
      legend.position = 'none',
      legend.title = element_text(hjust = 0.5),
      legend.text = element_text(angle = 45),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x = element_markdown(size = 10, color = 'black'),
      axis.text.x = element_text(angle = 0, vjust = -1),
      strip.background =element_rect(fill="white"),
      axis.text.y = element_text(hjust = 0.75)
    ) +
    NULL)

# antagonism in dyf-2, no significant ant or syn in N2
(raw_plot1 <- data_dyf_antshift1 %>%
    mutate(treatment_geno = paste(treatment, genotype)) %>%
    replace(is.na(.), 0) %>%
    ggplot(aes(x = treatment, y = raw_length)) +
    facet_wrap(~genotype) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 350), breaks = seq(0, 320, 60)) +
    #geom_violin(aes(), alpha = 1, size = 0.5) +
    # geom_point(aes(color = plate), alpha = 1, size = 0.5) +
    stat_summary(aes(color = treatment_geno), geom = 'point', fun = mean, size = 4) +
    stat_summary(aes(color = plate_well), geom = 'point', fun = mean, size = 2, alpha = 0.5) +
    geom_text(aes(label = ifelse(metadata_date %in% c('0131'),as.character(metadata_date),'')), hjust=10, vjust=0, size = 1) +
    scale_color_manual(values = c('black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 
                                  'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', '#008aff', 
                                  'black', '#008aff', 'black' , 'black', 'black', '#008aff' , 'black', 'black', '#008aff', 
                                  'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black')) +
    scale_fill_manual(values = c('black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black',
                                 'black', 'black', 'black', 'black', 'black', 'black', 'black', '#008aff', 'black', '#008aff',
                                 'black', 'black', 'black', '#008aff' , 'black', 'black', '#008aff', 'black', 'black', 
                                 'black', 'black', 'black', 'black', 'black', 'black', 'black')) +
    scale_x_discrete(limits = c('DMSO', "AZS", "IVM", 'AZS_IVM'), guide = guide_axis(angle = 25)) +
    coord_fixed(ratio=0.03) +
    ggtitle('AZS [50 µM]; IVM [0.01 µM]' ) +
    facet_grid(~factor(genotype, levels=c('Wild type', '*dyf-2(m160)*'))) +
    theme_nw2() +
    theme(legend.position = 'none', plot.tag = element_text(size = 13), plot.title = element_text(size = 11, hjust = 0.5)) +
    guides(fill='none', color='none') +
    labs(tag = "D") +
    ylab('') +
    xlab('Treatment') +
    NULL)

# antagonism in N2, no significant ant or syn in dyf-2
(raw_plot2 <- data_dyf_antshift2 %>%
    mutate(treatment_geno = paste(treatment, genotype)) %>%
    replace(is.na(.), 0) %>%
    ggplot(aes(x = treatment, y = raw_length)) +
    facet_wrap(~genotype) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 350), breaks = seq(0, 320, 60)) +
    #geom_violin(aes(), alpha = 1, size = 0.5) +
    # geom_point(aes(color = plate), alpha = 1, size = 0.5) +
    stat_summary(aes(color = treatment_geno), geom = 'point', fun = mean, size = 4) +
    stat_summary(aes(color = plate_well), geom = 'point', fun = mean, size = 2, alpha = 0.5) +
    geom_text(aes(label = ifelse(metadata_date %in% c('0131'),as.character(metadata_date),'')), hjust=10, vjust=0, size = 1) +
    scale_color_manual(values = c('black', 'black', 'black', 'black', '#008aff', 'black', 'black', 'black', 'black', 
                                  '#008aff', 'black', '#008aff', 'black', 'black', 'black', 'black', 'black', 'black', 
                                  'black', 'black', 'black' , 'black', 'black', 'black' , 'black', 'black', 'black', 
                                  '#008aff', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black')) +
    scale_fill_manual(values = c('black', 'black', 'black', 'black', '#008aff', 'black', 'black', 'black', 'black', 
                                 '#008aff', 'black', '#008aff', 'black', 'black', 'black', 'black', 'black', 'black', 
                                 'black', 'black', 'black' , 'black', 'black', 'black' , 'black', 'black', 'black', 
                                 '#008aff', 'black', 'black', 'black', 'black', 'black', 'black', 'black', 'black')) +
    scale_x_discrete(limits = c('DMSO', "AZS", "IVM", 'AZS_IVM'), guide = guide_axis(angle = 25)) +
    coord_fixed(ratio=0.03) +
    facet_grid(~factor(genotype, levels=c('Wild type', '*dyf-2(m160)*'))) +
    theme_nw2() +
    ggtitle('AZS [3.125 µM]; IVM [0.005 µM]') +
    theme(legend.position = 'none', plot.tag = element_text(size = 13), plot.title = element_text(size = 11, hjust = 0.5)) +
    guides(fill='none', color='none') +
    labs(tag = "C") +
    ylab('Raw length') +
    xlab('Treatment') +
    NULL)

# antagonism in dyf-2, no significant ant or syn in N2
# (norm_plot1 <- data_dyf_antshift1 %>%
#     mutate(treatment_geno = paste(treatment, genotype)) %>%
#     replace(is.na(.), 0) %>%
#     ggplot(aes(x = treatment, y = norm_DMSO)) +
#     facet_wrap(~genotype) +
#     scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5), breaks = seq(0, 1.5, 0.25)) +
#     #geom_violin(aes(), alpha = 1, size = 0.5) +
#     # geom_point(aes(color = plate), alpha = 1, size = 0.5) +
#     stat_summary(aes(color = treatment_geno), geom = 'point', fun = mean, size = 1.5) +
#     geom_text(aes(label = ifelse(metadata_date %in% c('0131'),as.character(metadata_date),'')), hjust=10, vjust=0, size = 1) +
#     scale_color_manual(values = c('black', 'black', '#008aff', 'black', 'black', 'black', 'black', 'black')) +
#     scale_fill_manual(values = c('black', 'black', '#008aff', 'black', 'black', 'black', 'black', 'black')) +
#     scale_x_discrete(limits = c('DMSO', "AZS", "IVM", 'AZS_IVM'), guide = guide_axis(angle = 25)) +
#     theme_nw2() +
#     theme(legend.position = 'none', plot.tag = element_text(size = 10)) +
#     guides(fill='none', color='none') +
#     labs(tag = "C") +
#     ylab('Normalized length') +
#     xlab('Treatment') +
#     NULL)

# antagonism in N2, no significant ant or syn in dyf-2
# (norm_plot2 <- data_dyf_antshift2 %>%
#     mutate(treatment_geno = paste(treatment, genotype)) %>%
#     replace(is.na(.), 0) %>%
#     ggplot(aes(x = treatment, y = norm_DMSO)) +
#     facet_wrap(~genotype) +
#     scale_y_continuous(expand = c(0, 0), limits = c(0, 1.5), breaks = seq(0, 1.5, 0.25)) +
#     #geom_violin(aes(), alpha = 1, size = 0.5) +
#     # geom_point(aes(color = plate), alpha = 1, size = 0.5) +
#     stat_summary(aes(color = treatment_geno), geom = 'point', fun = mean, size = 1.5) +
#     geom_text(aes(label = ifelse(metadata_date %in% c('0131'),as.character(metadata_date),'')), hjust=10, vjust=0, size = 1) +
#     # geom_segment(data = data_dyf_antshift2 %>% filter(treatment %in% c('AZS', 'IVM'), genotype = 'Wild type'),aes(x = treatment, xend  = treatment, y = 12, 
#     #                  yend = Temp), color = "blue", lwd = 1) +
#     scale_color_manual(values = c('black', 'black', 'black', '#008aff', 'black', 'black', 'black', 'black')) +
#     scale_fill_manual(values = c('black', 'black', 'black', '#008aff', 'black', 'black', 'black', 'black')) +
#     scale_x_discrete(limits = c('DMSO', "AZS", "IVM", 'AZS_IVM'), guide = guide_axis(angle = 25)) +
#     theme_nw2() +
#     theme(legend.position = 'none', plot.tag = element_text(size = 10)) +
#     guides(fill='none', color='none') +
#     labs(tag = "D") +
#     ylab('') +
#     xlab('Treatment') +
#     NULL)


#(curve_plots_combined <- plot_grid(raw_plot1, raw_plot2, norm_plot1, norm_plot2, align = "v", rel_heights = c(1, 1, 1, 1), nrow = 2))

(curve_plots_bottom <- plot_grid(raw_plot2, raw_plot1, align = "v", rel_heights = c(1, 1), nrow = 1))
(curve_plots_top <- plot_grid(heatmap_WT, heatmap_dyf, align = "v", rel_heights = c(1, 1), nrow = 1))
(SFig3 <- plot_grid(curve_plots_top, curve_plots_bottom, align = 'h', rel_heights = c(0.75, 1),  nrow = 2))


ggsave('/Users/elenagr/Desktop/Supp_Fig_4.tiff', SFig3, width = 200, height = 200, units = 'mm')
ggsave('/Users/elenagr/Desktop/Supp_Fig_4.pdf', SFig3, width = 200, height = 200, units = 'mm')





