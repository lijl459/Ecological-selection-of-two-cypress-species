# Convert SNPs in VCF format to FASTA alignments for phylogenetic analysis using vcf2phylip
vcf2phylip.py -i SNPs.vcf -f
# This command uses `vcf2phylip.py` to convert a VCF file containing SNPs to a FASTA format alignment file, 
# which can be used for subsequent phylogenetic analysis.

# Run RAxML for phylogenetic analysis
/data/soft/raxml/8.2.12/raxmlHPC-PTHREADS-AVX -T 8 -f a -p 12345 -x 12345 -m GTRGAMMA -N 200 -n $outname -s $input.fa
# This step runs RAxML to infer a phylogenetic tree using maximum likelihood.
# -T 8: Uses 8 threads for parallel processing.
# -f a: Performs a rapid bootstrap analysis followed by the search for the best-scoring tree.
# -p and -x: Random seeds for the analysis to ensure reproducibility.
# -m GTRGAMMA: Specifies the GTR model with Gamma-distributed rate heterogeneity.
# -N 200: Performs 200 bootstrap replicates.
# -n $outname: Specifies the name of the output files.
# -s $input.fa: Input alignment file in FASTA format.

# Load the environment variables
source ~/.bashrc
# This command loads any necessary environment variables and configurations from the `.bashrc` file.

# Run vcftools to prepare the data for plink
vcftools --gzvcf $1 --plink --out plink.1
# `vcftools` is used to convert the VCF file into PLINK format, which is required for further genetic analysis.

# Modify the PLINK map file to correct the allele coding
sed -i 's/^0/1/' plink.1.map
# This step modifies the `plink.1.map` file, changing allele coding from 0 to 1, as needed for the analysis.

# Perform LD pruning to remove correlated SNPs using PLINK
plink --noweb --file plink.1 --indep-pairwise 50 10 0.4 --out plink.2
# This command performs LD (linkage disequilibrium) pruning, removing SNPs that are highly correlated, with a window size of 50 SNPs, a step size of 10, and a r2 threshold of 0.4.

# Create a pruned dataset with the selected SNPs
plink --noweb --file plink.1 --extract plink.2.prune.in --make-bed --out plink.LD
# This creates a new dataset (`plink.LD`) based on the pruned SNPs from the previous step, saving the data in PLINK binary format.

# Convert the pruned data into a ped file for further analysis
plink --noweb --bfile plink.LD --recode12 --out final
# This converts the binary PLINK dataset (`plink.LD`) into a `.ped` file (`final.ped`) for use in downstream genetic analysis.

# Run Admixture for population structure analysis
for i in {1..10};
do
  admixture -j8 --cv final.ped $i | tee log$i.out
done
# This loop runs Admixture 10 times with different values of K (population clusters), where `-j8` specifies 8 CPU threads for parallel processing.
# The `--cv` option computes cross-validation errors for each K, helping to determine the optimal number of clusters.
# The `tee` command writes the output to both the console and log files (`log$i.out`).

# Run smartpca for principal component analysis (PCA)
smartpca.perl -i final.ped -a final.map -b final.pedind -o final.PCA -p final.PCA.plot -e final.PCA.eigenvalues -l final.PCA.log
# `smartpca` performs PCA on the genotype data to identify patterns of genetic variation.
# -i final.ped: Input PED file containing the genotype data.
# -a final.map: MAP file with SNP information.
# -b final.pedind: Individual information file for the PED file.
# -o final.PCA: Output PCA file.
# -p final.PCA.plot: Plot file for visualizing PCA results.
# -e final.PCA.eigenvalues: Eigenvalue output for the PCA analysis.
# -l final.PCA.log: Log file that records the steps and results of the PCA.
