# Run Stacks for single-end RADseq data:
|  |  |
|---|---|
|[Stacks V1](#v1) |[Stacks V2](#v2)|
## Stacks V1 <a name="v1"></a>
<pre>
## load stacks
module load Stacks/1.48-intel-2017A

## create a folder for stacks output
out_folder="stacks_out"
mkdir $out_folder

## format the barcode file
perl -ne 's/(\S+)\s+(\S+)/$2\t$1/; print $_' barcodes.txt  >barcodes_for-stacks.txt 

## step 1.  process_radtags
### remove low quality and adapter-containning reads
process_radtags -p fastq -i gzfastq -b barcodes_for-stacks.txt -o $out_folder -e pstI -E phred33 -r -c -q

## step 2. ustacks, generate unique stcks
## get sample names
cut -f 2 barcodes_for-stacks.txt >sample_names.txt 

names=`cat sample_names.txt`

counter=0
for name in $names;
do
  echo "counter: $counter"
  ustacks -t gzfastq -f $out_folder/${name}.R1.1.fq.gz -r -o $out_folder -i $counter -m 5 -M 3 -p 1
  let "counter+=1"
done
echo "ustacks is done"

## step 3. cstack, build catelog
samp=""
for name in $names;
do
  samp+="-s $out_folder/${name} "
done

batch_id=1
cstacks -b $batch_id $samp -o $out_folder -n 3  -p 15 
echo "cstacks is done"

## step 4. sstacks, map tags of each sample to catalog
for name in $names;
do
  sstacks -b $batch_id -c $out_folder/batch_${batch_id} -s $out_folder/$name -o $out_folder  
done
echo "sstacks is done"

## step 5. population, Calculate population statistics and export several output files
# generate a population map file
perl -ne 'chomp; print $_, "\t", substr($_, 0, 1), "\n"' sample_names.txt >pop_map.txt

populations -b $batch_id -P $out_folder -M pop_map.txt -r 2 -m 5 -e pstI -t 15 --genomic --fasta --vcf --structure --phylip

##
</pre>

##bsub file:

<pre>
#!/bin/bash
#BSUB -W 70:05                    # wall-clock time (hrs:mins)
#BSUB -L /bin/bash                # login shell    
#BSUB -n 1                        # number of tasks in job
#BSUB -R "span[ptile=1]"          # run one MPI task per node
#BSUB -R "rusage[mem=8000]"     # memory to reserve, in MB
#BSUB -J myjob                    # job name
#BSUB -o myjob.%J.%I.out             # output file name in which %J is replaced by the job ID
#BSUB -e myjob.%J.%I.err             # error file name in which %J is replaced by the job ID

date
## load stacks
module load Stacks/1.48-intel-2017A

## create a folder for stacks output
$out_folder="stacks_out"
mkdir $out_folder

## format the barcode file
perl -ne 's/(\S+)\s+(\S+)/$2\t$1/; print $_' barcodes.txt  >barcodes_for-stacks.txt 

## step 1.  process_radtags
### remove low quality and adapter-containning reads
process_radtags -p fastq -i gzfastq -b barcodes_for-stacks.txt -o $out_folder -e pstI -E phred33 -r -c -q

## step 2. ustacks, generate unique stcks
## get sample names
cut -f 2 barcodes_for-stacks.txt >sample_names.txt 

names=`cat sample_names.txt`

counter=0
for name in $names;
do
  echo "counter: $counter"
  ustacks -t gzfastq -f $out_folder/${name}.fq.gz -r -o $out_folder -i $counter -m 5 -M 3 -p 1
  let "counter+=1"
done
echo "ustacks is done!"

## step 3. cstack, build catelog
samp=""
for name in $names;
do
  samp+="-s $out_folder/${name} "
done

batch_id=1
cstacks -b $batch_id $samp -o $out_folder -n 3  -p 15 
echo "cstacks is done!"

## step 4. sstacks, map tags of each sample to catalog
for name in $names;
do
  sstacks -b $batch_id -c $out_folder/batch_${batch_id} -s $out_folder/$name -o $out_folder  
done
echo "sstacks is done!"

## step 5. population, Calculate population statistics and export several output files
# generate a population map file
perl -ne 'chomp; print $_, "\t", substr($_, 0, 1), "\n"' sample_names.txt >pop_map.txt

populations -b $batch_id -P $out_folder -M pop_map.txt -r 2 -m 5 -e pstI -t 15 --genomic --fasta --vcf --structure --phylip

</pre>

# Stacks V2 <a name="v2"></a>

## Load modules
```shell
module load parallel/20151222-intel-2015B
module load Stacks


INPUT_DIR=./fastq/

SAMPLE_DIR=samples
if [ ! -d $SAMPLE_DIR ]; then mkdir $SAMPLE_DIR; fi
OUTPUT_DIR=stacks2
if [ ! -d $OUTPUT_DIR ]; then mkdir $OUTPUT_DIR; fi


## format the barcode file
perl -ne 's/(\S+)\s+(\S+)/$2\t$1/; print $_' barcodes.txt  >barcodes_for-stacks.txt 

cut -f 2 barcodes_for-stacks.txt >accession_names.txt 
names=`cat accession_names.txt`
```

## step 1.  process_radtags, remove low quality and adapter-containning reads
### If demuxing has not been done yet
```shell
process_radtags -p $INPUT_DIR -i gzfastq -b barcodes_for-stacks.txt -o $SAMPLE_DIR -e pstI -E phred33 -r -c -q
```
### If demultiplexing was already done
```shell
  cmds=""
  for name in $names;
  do
    cmd="process_radtags -1 $INPUT_DIR/${name}_R1.fastq.gz  -2 $INPUT_DIR/${name}_R2.fastq.gz -o $SAMPLE_DIR -c -q disable_rad_check -i gzfastq"
    cmds+="$cmd;"
  done

  echo $cmds |tr ";" "\n" |parallel -j 20
```

## step 2, ustacks, generate unique stcks
```shell
cmd_file="ustacks.cmds"
if [ ! -s $cmd_file ]; then
  counter=1
  for name in $names;
  do
    echo "counter: $counter"
    cmd="ustacks -t fastq.gz -f $SAMPLE_DIR/${name}.fq.gz -r -o $OUTPUT_DIR -i $counter -m 5 -M 2 -p 5 --name $name"
    echo $cmd >>$cmd_file
    let "counter+=1"
  done
fi

cat $cmd_file | parallel -j 20 
```
 
## step 3, cstack, build catelog
```shell
names=`cat accession_names.txt`
samp=""
for name in $names;
do
  samp+="-s $OUTPUT_DIR/${name} "
done

date
echo "cstacks is started"

## Build the catalog of loci available in the metapopulation from the samples contained
## in the population map. To build the catalog from a subset of individuals, supply
## a separate population map only containing those samples.

cstacks  -n 2 -P $OUTPUT_DIR -M ./popmap -p 20

echo "cstacks is done"
date
```


## Step 4, Run sstacks. Match all samples supplied in the population map against the catalog.
```shell
date
echo "start sstacks"
sstacks  -P $OUTPUT_DIR  -M popmap -p 20
date
```

## Step 5, tsvbam and gstacks
### Run tsv2bam to transpose the data so it is stored by locus, instead of by sample. We will include paired-end reads using tsv2bam. tsv2bam expects the paired read files to be in the samples directory and they should be named consistently with the single-end reads,e.g. sample_01.1.fq.gz and sample_01.2.fq.gz, which is how process_radtags will output them.

Single-end reads:
```shell
tsv2bam -P $OUTPUT_DIR -M popmap -t 20
```
Pair-end reads:
```shell
tsv2bam -P $OUTPUT_DIR -M popmap --pe-reads-dir $SAMPLE_DIR -t 20
```

### Run gstacks: build a paired-end contig from the metapopulation data (if paired-reads provided), align reads per sample, call variant sites in the population, genotypes in each individual.

Gstacks:
```shell
gstacks -P $OUTPUT_DIR  -M popmap -t 20
```

## Step 6, population: Calculate population statistics and export several output files
### Run populations. Calculate Hardy-Weinberg deviation, population statistics, f-statistics, export several output files.

```shell
populations -P $OUTPUT_DIR -M popmap -r 0.80 --vcf --genepop --structure --fstats --hwe --phylip --fasta --write_single_snp  -t 20
```

