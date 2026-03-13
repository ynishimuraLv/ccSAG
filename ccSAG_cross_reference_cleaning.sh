#!/bin/sh
#$ -S /bin/sh
#$ -cwd

ccSAGdir=$(dirname $0)
source $ccSAGdir/path.txt
source $1

# Number of threads (default to 1 if not specified)
num_threads=${2:-1}

mkdir $Outdir -p
mkdir $Outdir/QC -p
mkdir $Outdir/Assemble -p
mkdir $Outdir/Mapping -p
mkdir $Outdir/Mapping_index -p

for file in `\ls $Seqdir | grep "R1_001.fastq$" | sed 's/_R1_001.fastq/\t/g'`;
do

######
# QC
######
if [ -s $Outdir/QC/${file}_QC_R1_001.fastq ]; then
    echo "Skipped ${file} QC."
else

    fastp -q 25 -u 50 -3 -Q 20 \
                -i $Seqdir/${file}_R1_001.fastq -I $Seqdir/${file}_R2_001.fastq \
                -o $Outdir/QC/${file}_QC_R1_001.fastq -O $Outdir/QC/${file}_QC_R2_001.fastq \
                -h $Outdir/QC/${file}.html -j $Outdir/QC/${file}.json \
                --thread $num_threads

fi

######
# Assemble
######
if [ -s $Outdir/Assemble/${file}_QC_contigs.fasta ]; then
    echo "Skipped ${file} assembly."
else

$spades_path $spades_option --threads $num_threads -1 $Outdir/QC/${file}_QC_R1_001.fastq -2 $Outdir/QC/${file}_QC_R2_001.fastq -o $Outdir/Assemble/${file}_QC_SPAdes
cp $Outdir/Assemble/${file}_QC_SPAdes/contigs.fasta $Outdir/Assemble/${file}_QC_contigs.fasta
rm $Outdir/Assemble/${file}_QC_SPAdes/ -r

fi

######
# make mapping index for cross reference
#####a
if [ -s $Outdir/Mapping_index/${file}_QC_contigs_500_index.amb ]; then
    echo "Skipped ${file} bwa indexing."
else

python $ccSAGdir/bin/longercontig.py $Outdir/Assemble/${file}_QC_contigs.fasta $Outdir/Assemble/${file}_QC_contigs_500.fasta 500
bwa index -p $Outdir/Mapping_index/${file}_QC_contigs_500_index $Outdir/Assemble/${file}_QC_contigs_500.fasta

fi

done


for file in `\ls $Outdir/QC | grep "_QC_R1_001.fastq$" | sed 's/_R1_001.fastq/\t/g'`;
do

######
# cross reference mapping
######

for index in `\ls $Outdir/QC | grep "_QC_R1_001.fastq$" | sed 's/_R1_001.fastq/\t/g'`;
do

if [ -s $Outdir/Mapping/${file}_${index}_uniq_classify.sam ]; then
    echo "Skipped ${file} mapping to ${index}."
else

if test $file != $index
then
 bwa mem -t $num_threads $Outdir/Mapping_index/${index}_contigs_500_index $Outdir/QC/${file}_R1_001.fastq $Outdir/QC/${file}_R2_001.fastq > $Outdir/Mapping/${file}_${index}.sam
 python $ccSAGdir/bin/get_primary_result_from_sam.py -i $Outdir/Mapping/${file}_${index}.sam -o $Outdir/Mapping/${file}_${index}_uniq.sam
 grep -e "^@" -v $Outdir/Mapping/${file}_${index}_uniq.sam > $Outdir/Mapping/${file}_${index}_uniq_classify.sam
 rm $Outdir/Mapping/${file}_${index}.sam
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
for index in `\ls $Outdir/QC | grep "_QC_R1_001.fastq$" | sed 's/_R1_001.fastq/\t/g'`;
do
if test ${file} != ${index}
then
 bwa mem -t $num_threads $Outdir/Mapping_index/${index}_contigs_500_index $Outdir/${file}_chimera/QC/${file}_cut_chimera.fastq > $Outdir/${file}_chimera/Mapping/${file}_${index}.sam
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

$spades_path $spades_option --threads $num_threads -s $Outdir/QC/cleaned_merge.fastq -o $Outdir/Assemble/cleaned_merge_SPAdes
cp $Outdir/Assemble/cleaned_merge_SPAdes/contigs.fasta $Outdir/$cleanedcontig

