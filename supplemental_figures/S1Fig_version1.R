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
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  
  annotation_logticks(sides = 'b', size = 0.25, outside = TRUE, 
                      short = unit(0.1, "cm"),
                      mid = unit(0.15, "cm"),
                      long = unit(0.2, "cm")) +
  
  scale_y_continuous(limits = c(0.2, 1.37), expand = c(0, 0), breaks = seq(0.4, 1.2, 0.2)) +


  scale_shape(guide = 'none') +
  
  coord_cartesian(clip = "off") +
  
  labs(x = 'Concentration (µM)', y = "Normalized length",
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
    axis.title.y = element_markdown(face = 'plain')) +
  # remove_legend() +
  NULL)

ggsave(here('/Users/elenagr/Desktop/Supp_2_bus5.pdf'), dc19_plot, width = 200, height = 80, units = 'mm')
ggsave(here('/Users/elenagr/Desktop/Supp_2_bus5.png'), dc19_plot, width = 200, height = 80, units = 'mm')

