# Environment setup
source activate eggnog  # Activate the eggnog environment
export EGGNOG_DATA_DIR=/data/blastdb/eggnog.db  # Set the path for the EGGNOG database

# Create Diamond BLAST database for Viridiplantae
create_dbs.py -m diamond --dbname Viridiplantae --taxa Viridiplantae  # Create a Diamond BLAST database for Viridiplantae

# Run eggnog-mapper to annotate the peptides for *Cupressus chengiana*, which were used to designed the Probes
emapper.py --cpu 20 -m diamond --dmnd_db /data/blastdb/eggnog.db/Viridiplantae.dmnd -i LJL17009_target.pep -o LJL17009_target --override 

# Extract GO terms in the format: gene_name GO_term
python anno2go.py LJL17009_target.emapper.annotations LJL17009_target.emapper.go.txt  # Extract and format GO annotations

# Perform GO enrichment analysis
Rscript topGo.R

# Use ggplot2 to plot the enrichment analysis results
Rscript plot.R
