#!/bin/sh
#$ -S /bin/sh
#$ -cwd

ccSAGdir=$(dirname $0)

# Default values
config_file=""
Outdir=""
Rawreaddir=""
cleaned_contig=""
outcontig=""
spades_option="-m 256 --sc --disable-rr --careful --disable-gzip-output "

# Help function
show_help() {
    cat << EOF
Usage: sh ccSAG_clamping.sh -o outdir -r rawreaddir -c cleaned_contig -C outcontig [options]

Description:
  Contig clamping pipeline that merges raw reads, performs assembly filtering, and connects contigs.

Options:
  -c, --config <file>           Load settings from config file (optional)
  -o, --outdir <dir>            Output directory (required)
  -r, --rawreaddir <dir>         Raw read directory (required)
  -c, --cleaned-contig <file>   Cleaned contig input file (required)
  -C, --outcontig <file>        Output contig filename (required)
  -O, --spades-option '<opts>'  SPAdes assembly options (default: "-m 256 --sc --disable-rr --careful --disable-gzip-output ")
  -h, --help                    Show this help message

Examples:
  # Using short options with required parameters
  sh ccSAG_clamping.sh -o /path/to/output -r /path/to/rawreads -c cleaned.fasta -C output.fasta

  # Using config file
  sh ccSAG_clamping.sh -c example_ccSAG_clamping.config

  # Custom SPAdes options
  sh ccSAG_clamping.sh -o /path/to/output -r /path/to/rawreads -c cleaned.fasta -C output.fasta -O "-m 512 --careful"

EOF
}

# Parse options
while [ $# -gt 0 ]; do
    case "$1" in
        --help|-h)
            show_help
            exit 0
            ;;
        --config)
            config_file="$2"
            shift 2
            ;;
        --outdir|-o)
            Outdir="$2"
            shift 2
            ;;
        --rawreaddir|-r)
            Rawreaddir="$2"
            shift 2
            ;;
        --cleaned-contig|-c)
            cleaned_contig="$2"
            shift 2
            ;;
        --outcontig|-C)
            outcontig="$2"
            shift 2
            ;;
        --spades-option|-O)
            spades_option="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Load config file if provided
if [ -n "$config_file" ]; then
    source "$config_file"
fi

# Check required parameters
if [ -z "$Outdir" ]; then
    echo "Error: --outdir option is required"
    echo "Use -h or --help for usage information"
    exit 1
fi
if [ -z "$Rawreaddir" ]; then
    echo "Error: --rawreaddir option is required"
    echo "Use -h or --help for usage information"
    exit 1
fi
if [ -z "$cleaned_contig" ]; then
    echo "Error: --cleaned-contig option is required"
    echo "Use -h or --help for usage information"
    exit 1
fi
if [ -z "$outcontig" ]; then
    echo "Error: --outcontig option is required"
    echo "Use -h or --help for usage information"
    exit 1
fi

# make contig clamps

if [ -s $Outdir/raw_SAG_merge_contigs.fasta ]; then
    echo "Skipped raw-merged-contig assembly."
else

echo -n > $Outdir/raw_merge_R1_001.fastq
echo -n > $Outdir/raw_merge_R2_001.fastq
for raw_reads in `ls $Rawreaddir | grep "_R1_001.fastq$" | sed 's/_R1_001.fastq/\t/g'`;
do
cat $Rawreaddir/${raw_reads}_R1_001.fastq >> $Outdir/raw_merge_R1_001.fastq
cat $Rawreaddir/${raw_reads}_R2_001.fastq >> $Outdir/raw_merge_R2_001.fastq
done

spades.py $spades_option -1 $Outdir/raw_merge_R1_001.fastq -2 $Outdir/raw_merge_R2_001.fastq -o $Outdir/raw_merge_SPAdes
cp $Outdir/raw_merge_SPAdes/contigs.fasta $Outdir/raw_SAG_merge_contigs.fasta
rm $Outdir/raw_merge_SPAdes/ -r
rm $Outdir/raw_merge_R1_001.fastq $Outdir/raw_merge_R2_001.fastq

fi

python $ccSAGdir/bin/longercontig.py $cleaned_contig $Outdir/cleaned_500.fasta 500
python $ccSAGdir/bin/longercontig.py $Outdir/raw_SAG_merge_contigs.fasta $Outdir/noncleaned_500.fasta 500

# make first_fished_contigs
makeblastdb -in $Outdir/cleaned_500.fasta -dbtype nucl -out $Outdir/fish_index -hash_index
blastn -query $Outdir/noncleaned_500.fasta -out $Outdir/fish_blastresult -db $Outdir/fish_index -perc_identity 99 -outfmt 6 -max_target_seqs 1
python $ccSAGdir/bin/chimera_cross_reference_third_uniq_blastresult.py $Outdir/fish_blastresult $Outdir/fish_uniq_blastresult 250 on
python $ccSAGdir/bin/chimera_cross_reference_third_fishing_from_noncleaned_merge.py $Outdir/fish_uniq_blastresult $Outdir/noncleaned_500.fasta $Outdir/fished_contigs.fasta 250
rm $Outdir/fish_index* $Outdir/fish_blastresult $Outdir/fish_uniq_blastresult

# remove short contigs from cleaned_merge_contigs
makeblastdb -in $Outdir/fished_contigs.fasta -dbtype nucl -out $Outdir/fished_index -hash_index
blastn -query $Outdir/cleaned_500.fasta -out $Outdir/cleaned_blastresult -db $Outdir/fished_index -perc_identity 99 -outfmt 6 -max_target_seqs 1
python $ccSAGdir/bin/chimera_cross_reference_third_uniq_blastresult.py $Outdir/cleaned_blastresult $Outdir/cleaned_uniq_blastresult 250 off
python $ccSAGdir/bin/chimera_cross_reference_third_remove_short_contig_in_cleaned.py $Outdir/cleaned_uniq_blastresult $Outdir/cleaned_500.fasta $Outdir/cleaned_long.fasta $Outdir/cleaned_keep.fasta 250
rm $Outdir/cleaned_blastresult $Outdir/cleaned_uniq_blastresult

# connect contigs
makeblastdb -in $Outdir/cleaned_long.fasta -dbtype nucl -out $Outdir/cleaned_long_index -hash_index
blastn -query $Outdir/cleaned_long.fasta -out $Outdir/cleaned_to_fished_blastresult -db $Outdir/fished_index -perc_identity 99 -outfmt 6 -max_target_seqs 2
python $ccSAGdir/bin/chimera_cross_reference_third_connect_contigs.py $Outdir/cleaned_to_fished_blastresult $Outdir/cleaned_long.fasta $Outdir/fished_contigs.fasta $Outdir/connected_contigs.fasta 150
cat $Outdir/cleaned_keep.fasta $Outdir/connected_contigs.fasta > $Outdir/$outcontig
rm $Outdir/*index.* $Outdir/cleaned_to_fished_blastresult $Outdir/cleaned_long.fasta $Outdir/fished_contigs.fasta $Outdir/noncleaned_500.fasta $Outdir/cleaned_500.fasta $Outdir/connected_contigs.fasta $Outdir/cleaned_keep.fasta
