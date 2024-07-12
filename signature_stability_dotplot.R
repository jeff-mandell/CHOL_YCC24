# Purpose: Plot mean contribution and detection frequency of mutational signatures in CCA

library(data.table)
library(ggplot2)


for_plot = fread('output/signature_comparisons.txt')
for_plot[, short_signature := sub('SBS', '', signature)]
for_plot = rbind(for_plot[, .(comparison, short_signature, cca_type = grp1, detection_frequency = grp1.detection, contribution = grp1.mean, signif_label)],
                        for_plot[, .(comparison, short_signature, cca_type = grp2, detection_frequency = grp2.detection, contribution = grp2.mean, signif_label)])

# Remove signatures that are infrequently detected and have low contributions.
minor_sigs = for_plot[, .(highest_detection = max(detection_frequency), highest_contrib = max(contribution)),
                                   by = 'short_signature'][highest_detection < .2 & highest_contrib < .05, short_signature]
for_plot = for_plot[! short_signature %in% minor_sigs]

# Get sample counts
signature_attributions = fread('output/final_unblended_signature_weights.txt')
sample_counts_by_type = signature_attributions[, .N, by = 'cca_type']
ehc_count = sample_counts_by_type[cca_type != 'IHC', sum(N)]
sample_counts_by_type[cca_type == 'EHC', N := ehc_count]
sample_counts_by_type[, count_label := paste0('(n = ', N, ')')]
for_plot[sample_counts_by_type, count_label := count_label, on = 'cca_type' ]

# Build subtype labels
for_plot[, type_label := fcase(cca_type == 'IHC', paste0('iCCA\n', count_label),
                             cca_type == 'PHC', paste0('pCCA\n', count_label),
                             cca_type == 'DCC', paste0('dCCA\n', count_label),
                             cca_type == 'EHC', paste0('eCCA\n', count_label))]

for_type_labeling = setNames(unique(for_plot$type_label), unique(for_plot$cca_type))
for_plot[, cca_type := factor(cca_type, levels = c('IHC', 'EHC', 'PHC', 'DCC'))]

# Will denote significance for two sets of group comparisons
signif_labels = unique(for_plot[comparison %in% c('IHC-EHC', 'DCC-PHC'), .(comparison, short_signature, signif_label)])
signif_labels[, x := fcase(comparison == 'IHC-EHC', 1.5,
                           comparison == 'DCC-PHC', 3.5)]

# Sort by (mostly) numerical signature
for_plot[, sig_num := as.numeric(gsub('[^0-9]', '', short_signature))] # remove all non-numeric
for_plot[, sig_suffix := sub('^\\d+', '', short_signature)]
for_plot = for_plot[order(-sig_num, -sig_suffix)]
for_plot[, c('signif_label', 'sig_num', 'sig_suffix', 'comparison') := NULL]
for_plot[, short_signature := factor(short_signature, levels = unique(for_plot$short_signature))]

dotplot = ggplot(data = for_plot, 
                 aes( y = short_signature, x = cca_type, fill = detection_frequency, size = contribution)) +
  geom_point(shape = 21) +
  geom_text(data = signif_labels, aes(label = signif_label, x = x, y = short_signature), 
            nudge_y = -0.15, size = 7, na.rm = TRUE, inherit.aes = FALSE) +
  scale_x_discrete(labels = for_type_labeling) + scale_y_discrete() +
  scale_fill_distiller(palette = "RdYlBu", limits = c(0, 1)) +
  scale_size(range = c(1, 8)) +
  labs(x = "Subtype", y = "SBS signature", fill = "Detection frequency", size = "Mean contribution") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 9), 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(linewidth = 0.05, color = "grey"),
        panel.border = element_blank(), axis.ticks.y = element_blank()) +
  guides(size = guide_legend(reverse = TRUE))

# save plot
ggsave(file = 'figures/figure3.png', plot = dotplot, height = 4, width = 9, units = 'in')

