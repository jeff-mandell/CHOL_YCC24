# Purpose: Perform group comparisons on signature attributions between CCA subtypes
# Used by: figures/signature_stability_dotplot.R

# Load in required libraries
library(data.table)

# Load in bootstrapped attribution data (see signature_attribution_analysis.R):
# Contains the proportional attribution each signature contributes
# to each sample's total attribution, averaged across boots (370 samples)
signature_attributions = fread('output/final_unblended_signature_weights.txt')

sig_cols = names(signature_attributions)[names(signature_attributions) %like% 'SBS']

# Exclude signatures with no attributions (that is, artifact signatures and non-CCA signatures) from testing
all_zero_sigs = names(which(signature_attributions[, sapply(.SD, function(x) all(x == 0)), 
                                                   .SDcols = sig_cols]))
sig_cols = setdiff(sig_cols, all_zero_sigs)


# Function to run Mann-Whitney U on two groups subsetted from signature_attributions
compare_groups = function(grps) {
  stopifnot(is.list(grps), length(grps) == 2, uniqueN(names(grps)) == 2)
  results = sapply(sig_cols, function(x) {
    wilcox.test(
      grps[[1]][[x]],
      grps[[2]][[x]],
      alternative = 'two.sided',
      exact = FALSE
    )$p.value
  })
  
  quantile_and_mean = function(x) c(quantile(x), c(mean = mean(x)))
  grp1_quantiles = as.data.table(t(grps[[1]][, sapply(.SD, quantile_and_mean), .SDcols = sig_cols]), 
                                 keep.rownames = 'signature')
  grp2_quantiles = as.data.table(t(grps[[2]][, sapply(.SD, quantile_and_mean), .SDcols = sig_cols]),
                                 keep.rownames = 'signature')
  quantile_names = c('0%', '25%', '50%', '75%', '100%', 'mean')
  setnames(grp1_quantiles, quantile_names, paste0('grp1.', c('min', '25', 'median', '75', 'max', 'mean')))
  setnames(grp2_quantiles, quantile_names, paste0('grp2.', c('min', '25', 'median', '75', 'max', 'mean')))
  comparison_name = paste(names(grps), collapse = '-')
  all_quantiles = merge.data.table(grp1_quantiles, grp2_quantiles, by = 'signature')
  output = data.table(comparison = comparison_name, signature = sig_cols, 
                      grp1 = names(grps)[1], grp2 = names(grps)[2], pval = results)
  output = merge.data.table(output, all_quantiles, by = 'signature')
  return(output)
}

IHC = signature_attributions[cca_type == 'IHC']
PHC = signature_attributions[cca_type == 'PHC']
DCC = signature_attributions[cca_type == 'DCC']
EHC = signature_attributions[cca_type %in% c('EHC', 'PHC', 'DCC')]

signature_comparisons = rbindlist(list(compare_groups(list(IHC = IHC, EHC = EHC)),
                                       compare_groups(list(IHC = IHC, PHC = PHC)),
                                       compare_groups(list(IHC = IHC, DCC = DCC)),
                                       compare_groups(list(DCC = DCC, PHC = PHC))))

# Merge in signature detection frequency from bootstraps
raw_mp_out = fread('output/bootstrapped_mp_out.txt')
raw_mp_out = raw_mp_out[Unique_Patient_Identifier %in% signature_attributions$Unique_Patient_Identifier]
raw_mp_out[signature_attributions, cca_type := cca_type, on = 'Unique_Patient_Identifier']

signature_names = names(raw_mp_out)[names(raw_mp_out) %like% 'SBS']
ihc_phc_dcc_detection = raw_mp_out[cca_type %in% c('IHC', 'PHC', 'DCC'), 
                                   lapply(.SD, function(x) mean(x > 0)), 
                                   .SDcols = signature_names, by = 'cca_type']
ehc_detection = raw_mp_out[cca_type %in% c('PHC', 'DCC', 'EHC'),
                           lapply(.SD, function(x) mean(x > 0)),
                           .SDcols = signature_names]
ehc_detection$cca_type = 'EHC'
detection = melt(rbind(ihc_phc_dcc_detection, ehc_detection),
                 id.vars = 'cca_type', value.factor = FALSE, variable.name = 'signature')

signature_comparisons[detection, grp1.detection := value, on = c(grp1 = 'cca_type', 'signature')]
signature_comparisons[detection, grp2.detection := value, on = c(grp2 = 'cca_type', 'signature')]

# Significance labels for plots
signature_comparisons[pval < .05, signif_label := '*']
signature_comparisons[pval < .01, signif_label := '**']
signature_comparisons[pval < .001, signif_label := '***']

fwrite(signature_comparisons, 'output/signature_comparisons.txt', sep = "\t")






