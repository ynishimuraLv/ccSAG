#!/bin/sh
#$ -S /bin/sh
#$ -cwd

ccSAGdir=$(dirname $0)

# Default values
num_threads=1
config_file=""
Outdir=""
Seqdir=""
cleanedcontig=""
spades_option="-m 256 --sc --disable-rr --careful --disable-gzip-output "

# Help function
show_help() {
    cat << EOF
Usage: bash ccSAG_cross_reference_cleaning.sh -o outdir -s seqdir -C outputcontig [options]

Description:
  Cross-reference genome assembly cleaning pipeline with quality control, assembly, and mapping.

Options:
  -c, --config <file>           Load settings from config file (optional)
  -o, --outdir <dir>            Output directory for results (required)
  -s, --seqdir <dir>            Input sequence directory (required)
  -C, --cleaned-contig <file>   Output cleaned contig filename (required)
  -t, --threads <num>           Number of threads to use (default: 1)
  -O, --spades-option '<opts>'  SPAdes assembly options (default: "-m 256 --sc --disable-rr --careful --disable-gzip-output ")
  -h, --help                    Show this help message

Examples:
  # Using short options with required parameters
  bash ccSAG_cross_reference_cleaning.sh -o /path/to/output -s /path/to/raw/data -C cleaned.fasta

  # With thread specification
  bash ccSAG_cross_reference_cleaning.sh -o /path/to/output -s /path/to/raw/data -C cleaned.fasta -t 4

  # Using config file
  bash ccSAG_cross_reference_cleaning.sh -c example_ccSAG_cross_reference_cleaning.config -t 4

  # Custom SPAdes options
  bash ccSAG_cross_reference_cleaning.sh -o /path/to/output -s /path/to/raw/data -C cleaned.fasta -O "-m 512 --careful"

EOF
}

# Parse options
while [ $# -gt 0 ]; do
    case "$1" in
        --help|-h)
            show_help
            exit 0
            ;;
        --config|-c)
            config_file="$2"
            shift 2
            ;;
        --threads|-t)
            num_threads="$2"
            shift 2
            ;;
        --outdir|-o)
            Outdir="$2"
            shift 2
            ;;
        --seqdir|-s)
            Seqdir="$2"
            shift 2
            ;;
        --cleaned-contig|-C)
            cleanedcontig="$2"
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
if [ -z "$Seqdir" ]; then
    echo "Error: --seqdir option is required"
    echo "Use -h or --help for usage information"
    exit 1
fi
if [ -z "$cleanedcontig" ]; then
    echo "Error: --cleaned-contig option is required"
    echo "Use -h or --help for usage information"
    exit 1
fi

mkdir $Outdir -p
mkdir $Outdir/QC -p
mkdir $Outdir/Assemble -p
mkdir $Outdir/Mapping -p
mkdir $Outdir/Mapping_index -p

######
# QC
######

for r1 in "$Seqdir"/*_R1*.f*q*; do
    r2=${r1/_R1/_R2}
    sample=$(basename "$r1" | sed -E 's/_R1(_001)?\.(fastq|fq)(\.gz)?$//')
    ######
    # QC
    ######
    if [ -s $Outdir/QC/${sample}.html ]; then
        echo "Skipped ${sample} QC."
    else
    fastp -q 25 -u 50 -3 -Q 20 \
          -i $r1 -I $r2 \
          -o $Outdir/QC/${sample}_QC_R1_001.fastq -O $Outdir/QC/${sample}_QC_R2_001.fastq \
          -h $Outdir/QC/${sample}.html -j $Outdir/QC/${sample}.json \
          --thread $num_threads

    fi

    ######
    # Assemble
    ######
    if [ -s $Outdir/Assemble/${sample}_QC_contigs.fasta ]; then
        echo "Skipped ${sample} assembly."
    else

        spades.py $spades_option --threads $num_threads \
                  -1 $Outdir/QC/${sample}_QC_R1_001.fastq \
                  -2 $Outdir/QC/${sample}_QC_R2_001.fastq \
                  -o $Outdir/Assemble/${sample}_QC_SPAdes
        cp $Outdir/Assemble/${sample}_QC_SPAdes/contigs.fasta $Outdir/Assemble/${sample}_QC_contigs.fasta
        rm $Outdir/Assemble/${sample}_QC_SPAdes/ -r

    fi

    ######
    # make mapping index for cross reference
    #####a
    if [ -s $Outdir/Mapping_index/${sample}_QC_contigs_500_index.amb ]; then
        echo "Skipped ${sample} bwa indexing."
    else

        python $ccSAGdir/bin/longercontig.py $Outdir/Assemble/${sample}_QC_contigs.fasta $Outdir/Assemble/${sample}_QC_contigs_500.fasta 500
        bwa index -p $Outdir/Mapping_index/${sample}_QC_contigs_500_index $Outdir/Assemble/${sample}_QC_contigs_500.fasta

    fi

done


for r1 in "$Outdir/QC"/*_QC_R1.fq.gz; do
    file=$(basename "$r1" | sed -E 's/_QC_R1.fq.gz$//')
    r2=${r1/_R1/_R2}

    ######
    # cross reference mapping
    ######

    for file1 in "$Outdir"/QC/*_QC_R1.fq.gz; do
        index=$(basename "$file1" | sed -E 's/_QC_R1.fq.gz$//')
        classify="$Outdir/Mapping/${file}_${index}_uniq_classify.sam"
        if [ -s "$classify" ]; then
            echo "Skipped ${file} mapping to ${index}."
        else
            if [[ "$file" != "$index" ]]; then
                samfile="$Outdir/Mapping/${file}_${index}.sam"
                samfile_uniq="$Outdir/Mapping/${file}_${index}_uniq.sam"
                echo "Mapping ${file}_${index}"
                bwa mem -t $num_threads $Outdir/Mapping_index/${index}_QC_contigs_500_index ${r1} ${r2} > "$samfile"
                python $ccSAGdir/bin/get_primary_result_from_sam.py -i $samfile -o $samfile_uniq
                grep -e "^@" -v $samfile_uniq > $classify
                rm $samfile
            fi
        fi
    done

    python $ccSAGdir/bin/classify_chimera_read.py $Outdir/Mapping ${file}

    mkdir -p $Outdir/${file}_chimera
    mkdir -p $Outdir/${file}_chimera/QC
    mkdir -p $Outdir/${file}_chimera/Mapping

    mv $Outdir/Mapping/${file}_*.fastq $Outdir/QC
    mv $Outdir/QC/${file}_cut_chimera.fastq $Outdir/${file}_chimera/QC
    echo -n > $Outdir/QC/${file}_multicut_chimera.fastq

    while test $(wc -l < $Outdir/${file}_chimera/QC/${file}_cut_chimera.fastq) != 0
        do
            for file2 in "$Outdir"/QC/*_QC_R1.fq.gz; do
                index=$(basename "$file2" | sed -E 's/_QC_R1.fq.gz$//')
                    if test ${file} != ${index}
                    then
                        echo "Map ${index} ${file}_chimera"
                        bwa mem -t $num_threads $Outdir/Mapping_index/${index}_QC_contigs_500_index $Outdir/${file}_chimera/QC/${file}_cut_chimera.fastq > $Outdir/${file}_chimera/Mapping/${file}_${index}.sam
                        python $ccSAGdir/bin/get_primary_result_from_sam.py -i $Outdir/${file}_chimera/Mapping/${file}_${index}.sam -o $Outdir/${file}_chimera/Mapping/${file}_${index}_uniq.sam
                        grep -e "^@" -v $Outdir/${file}_chimera/Mapping/${file}_${index}_uniq.sam > $Outdir/${file}_chimera/Mapping/${file}_${index}_uniq_classify.sam
                    fi
            done

        python $ccSAGdir/bin/cut_chimera_read.py $Outdir/${file}_chimera/Mapping ${file}
        cat $Outdir/${file}_chimera/Mapping/${file}_normal_001.fastq >> $Outdir/QC/${file}_multicut_chimera.fastq
        mv $Outdir/${file}_chimera/Mapping/${file}_cut_chimera.fastq $Outdir/${file}_chimera/QC
    done
    rm $Outdir/${file}_chimera -r

    cat $Outdir/QC/${file}_multicut_chimera.fastq $Outdir/QC/${file}_normal_R1_001.fastq $Outdir/QC/${file}_normal_R2_001.fastq > $Outdir/QC/${file}_cleaned.fastq

done


######
# Assemble using cleaned fastq
######

echo -n > $Outdir/QC/cleaned_merge.fastq
for cleaned_reads in `\ls $Outdir/QC | grep "_QC_cleaned.fastq$"`;
    do
    cat $Outdir/QC/$cleaned_reads >> $Outdir/QC/cleaned_merge.fastq
done

spades.py $spades_option --threads $num_threads -s $Outdir/QC/cleaned_merge.fastq -o $Outdir/Assemble/cleaned_merge_SPAdes
cp $Outdir/Assemble/cleaned_merge_SPAdes/contigs.fasta $Outdir/$cleanedcontig

