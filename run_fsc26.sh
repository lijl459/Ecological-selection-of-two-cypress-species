# Perl scripts were downloaded from two sources:
# 1. GitHub repository: https://github.com/wk8910/bio_tools/
# 2. Ru D, Sun Y, Wang D, et al. Population genomic analysis reveals that homoploid hybrid speciation can be a lengthy process. Mol Ecol. 2018; 27: 4875–4887. https://doi.org/10.1111/mec.14909


# get_0fold-4fold_sites.pl extracts 0-fold and 4-fold degenerate sites from a VCF file.
# These sites are useful for analyzing selective pressure and neutral evolution.

# vcf2maf.pl generates the two-dimensional site frequency spectrum (SFS) from the VCF file.
# The SFS is important for studying allele frequency distributions across different populations or species.

# 26.fastsimcoal.pl is used to run fastsimcoal26, performing 50 replications to simulate the demographic model.

# 01.best.pl selects the best model based on the lowest AIC (Akaike Information Criterion) value.
# The AIC is used for model selection, helping identify the model that best fits the observed data while penalizing for model complexity.
