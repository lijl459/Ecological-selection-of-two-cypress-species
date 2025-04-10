# Run Trimmomatic to filter reads
cd /data/rawdata/Hybseq/Cupressus/torulosa/trim
java -jar /data/soft/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 -phred33 /data/rawdata/Hybseq/Cupressus/torulosa/clean/Lars-E5A/FCH2KTKCCX2_L5_CUPpbqTAADKAAA-133_1.fq.gz /data/rawdata/Hybseq/Cupressus/torulosa/clean/Lars-E5A/FCH2KTKCCX2_L5_CUPpbqTAADKAAA-133_2.fq.gz Lars-E5A_trimmed-paired.1.fq.gz Lars-E5A_trimmed-UNpaired.1.fq.gz Lars-E5A_trimmed-paired.2.fq.gz Lars-E5A_trmmed-UNpaired.2.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20
# Trimmomatic is used to filter raw reads by removing low-quality bases and trimming short reads.
# The command trims paired-end reads by removing bases at the 3' end with a quality score below 3.
# It also uses a sliding window approach (4 bases) and requires a minimum average quality score of 20 for each window.

# Run HybPiper to assemble the reference

# Run the main HybPiper script for 20 CPUs
while read i
do
  /usr/bin/time -vo $i.time python2 /data/soft/HybPiper/reads_first.py -r /data/user016/fastq/$i*fq -b LJL17009_target_cds.fas --prefix $i --bwa --cpu 20
done < namelist.txt
# This step runs the main HybPiper script for each sample in the namelist.txt file. It uses 20 CPU cores for parallel processing and aligns reads using BWA.

# Generate the seq_lengths.txt file
python2 /data/soft/HybPiper/get_seq_lengths.py LJL17009_target_cds.fas namelist.txt dna > seq_lengths.txt
# This script calculates the sequence lengths for the assembled CDS sequences and stores the result in seq_lengths.txt.

# Test for paralogs in the assembly
while read i
do
  python2 /data/soft/HybPiper/paralog_investigator.py $i
done < namelist.txt
# This script tests for potential paralogs in the assemblies by analyzing sequence similarity and clustering.

# Run IntronRate for intron analysis
while read i
do
  python2 /data/soft/HybPiper/intronerate.py --prefix $i
done < namelist.txt
# IntronRate is used to assess the intron-exon structure of the assembly and calculate the intron-exon rate for each sequence.

# Retrieve supercontigs
python retrieve_sequences_multiprocessing.py LJL17009_target_cds.fas retrieve_sequences_multiprocessing.py supercontig 20
# This script retrieves supercontigs from the assembled sequences, generating a reference supercontig fasta file (`Lars-E5_supercontig.fasta`).
# The process is parallelized with 20 threads to speed up sequence retrieval.

#### Map reads to the reference
# Run bwa mem to map the trimmed paired-end reads to the reference supercontig
bwa mem -t 5 -R '@RG\\tID:Lars-E5A\\tSM:Lars-E5A\\tLB:Illumina' Lars-E5A Lars-E5A_trimmed-paired.1.fq.gz Lars-E5A_trimmed-paired.2.fq.gz > ./0bwa/Lars-E5A.sam && touch ./ok/Lars-E5A.mem.ok
# BWA MEM is used to map the paired-end reads to the reference genome (`Lars-E5_supercontig.fasta`), with read group information for downstream processing.
# The `-t 5` flag specifies the number of threads to use (5 in this case).

# Convert SAM to BAM, filter for mapped reads, and create a BAM file
samtools view -@ 5 -q 20 -F 4 -bS ./0bwa/Lars-E5A.sam -o ./0bwa/Lars-E5A.bam && touch ./ok/Lars-E5A.view.ok
# This step converts the SAM file to BAM format, filtering out unmapped reads (with the `-F 4` flag) and ensuring that only high-quality reads (quality score >= 20) are kept.
# The `-@ 5` flag uses 5 threads for faster processing.

# Sort the BAM file
samtools sort -@ 5 ./0bwa/Lars-E5A.bam -o ./0bwa/Lars-E5A.sort.bam && touch ./ok/Lars-E5A.sort.ok
# Sort the BAM file by the coordinates of the reads for downstream processing.

# Mark and remove PCR duplicates with Picard
java -jar $picard MarkDuplicates REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 INPUT=./0bwa/Lars-E5A.sort.bam OUTPUT=1picard/Lars-E5A.picard.bam METRICS_FILE=1picard/Lars-E5A.picard.bam.metrics && touch ./ok/Lars-E5A.picard.ok
# Picard's MarkDuplicates tool is used to identify and mark duplicate reads resulting from PCR amplification.
# The `REMOVE_DUPLICATES=true` flag ensures that duplicates are removed, improving the quality of downstream variant calling.

# Index the BAM file
samtools index 1picard/Lars-E5A.picard.bam && touch ./ok/Lars-E5A.index.ok
# This step creates an index file for the BAM file, allowing for faster random access to the file in subsequent analyses.

# Run GATK RealignerTargetCreator to find realignment regions around indels
java -jar /data/soft/gatk/gatk-3.7/GenomeAnalysisTK.jar -R Lars-E5_supercontig.fasta -T RealignerTargetCreator -I 1picard/Lars-E5.picard.bam -o 2RealigneRegion/Lars-E5.picard.bam.intervals && touch ./ok/Lars-E5.RealigneRegion.ok
# GATK RealignerTargetCreator identifies regions of the genome around indels (insertion-deletions) that need realignment for more accurate variant calling.
# The resulting interval file is used to specify the regions for the next realignment step.

# Run GATK IndelRealigner to realign reads around the indels
java -jar /data/soft/gatk/gatk-3.7/GenomeAnalysisTK.jar -R Lars-E5_supercontig.fasta -T IndelRealigner -targetIntervals 2RealigneRegion/Lars-E5.picard.bam.intervals -I 1picard/Lars-E5.picard.bam -o 3finalbam/Lars-E5.final.bam && touch ./ok/Lars-E5.indelRealine.ok
# GATK IndelRealigner realigns the reads around the indel regions identified earlier, improving the accuracy of variant calling.
# The final realigned BAM file is saved in the `3finalbam` directory, and the `touch` command creates a flag file indicating successful completion of the step.


# Run samtools mpileup and bcftools to generate variants
samtools mpileup -t AD,ADF,ADR,DP,SP -Q 20 -q 20 -ugf Lars-E5_supercontig.fasta -b bam1.list | bcftools call -m > variants.vcf
# This step uses `samtools mpileup` to generate a pileup of the mapped reads with specific tags (AD, ADF, ADR, DP, SP), filtering for reads with a quality score >= 20.
# The `bcftools call -m` command is used to call variants from the pileup and output the results in VCF format (`variants.vcf`).





