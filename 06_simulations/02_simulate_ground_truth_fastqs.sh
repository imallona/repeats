#!/bin/bash
##
## gets a randomly picked subsequence inside each line of a GTF record (repeats file)
##  (maybe it should do it on beds, because gtfs lines from GTFs are nested but anyway)
## and mimics a fastq record with uniform base quality values
##
##
##
# define read length (fixed)
# LENGTH=92
# NTHREADS=20
# QUAL_SYMBOL=F
# # get rmsk coordinates
# GTF=~/repeats_sc/annotation/mm10_rmsk_TE.gtf.gz

# GENOME_FASTA=~/repeats_sc/annotation/Mus_musculus.GRCm38.dna.primary_assembly.fa

# BEDTOOLS=~/soft/bedtools/bedtools-2.29.2/bin/bedtools


# sample one read length per rmsk unit

# if the repeat is shorter than the specified read length,
#    report the full coordinates of the repeat
# else
#    get a random substring of the repeat matching the read length
#


simulate_cdna () {
    local LENGTH="$1"
    local QUAL_SYMBOL="$2"
    local GTF="$3"
    local GENOME_FASTA="$4"
    local NTHREADS="$5"
    local BEDTOOLS="$6"
    local OUTPUT="$7"

    zcat $GTF |   awk -v RL="$LENGTH" '{OFS=FS="\t"}{
$0=$0;
L=$6
chrom=$1
start=$4
end=$5
strand=$7
name=$9
srand(NR);
if (RL >= L)
   print chrom,start,end,name,".",strand
else
   min=$4
   max=$5-$RL
   R=int(min+rand()*(max-min+1))
   print chrom,R,R+RL,name,".",strand
}' > sampled_"$LENGTH"_intervals.bed

    sed -i 's/gene_id//g ; s/transcript_id//g; s/class_id// ; s/family_id// ; s/"//g ; s/ //g' sampled_"$LENGTH"_intervals.bed

    pigz --decompress -p $NTHREADS "$GENOME_FASTA".gz

    $BEDTOOLS getfasta  -fi $GENOME_FASTA \
              -bed sampled_"$LENGTH"_intervals.bed \
              -tab -name > sampled_"$LENGTH"_intervals.fa

    pigz -p $NTHREADS "$GENOME_FASTA"

    # fastq-ify

    awk -v qual="$QUAL_SYMBOL" 'OFS=FS="\t" {
$0=$0;
L=length($2)
fa=$2
name=$1

printf name"\n"; 
printf fa"\n";
printf "+\n"
for(c=0;c<L;c++) printf qual; printf "\n"
}' sampled_"$LENGTH"_intervals.fa > "$OUTPUT"

    rm sampled_"$LENGTH"_intervals.bed
    rm sampled_"$LENGTH"_intervals.fa
    # now, CBs and UMIs must be simulated too
}

# simulate_cdna $LENGTH $QUAL_SYMBOL $GTF $GENOME_FASTA $NTHREADS $BEDTOOLS foo.fastq
