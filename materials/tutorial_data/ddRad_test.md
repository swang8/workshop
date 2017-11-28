#Run Stacks for pairend ddRAD data:

<pre>
## load stacks
module load Stacks

## create a folder for stacks output
out_folder="stacks_out"
mkdir $out_folder

## format the barcode file
perl -ne 's/(\S+)\s+(\S+)/$2\t$1/; print $_' ddRad_barcodes.txt  >ddRad_barcodes_for-stacks.txt 
head ddRad_barcodes.txt  ddRad_barcodes_for-stacks.txt

## step 1.  process_radtags
### remove low quality and adapter-containning reads
process_radtags -p fastq -P -i gzfastq -b ddRad_barcodes_for-stacks.txt -o $out_folder --renz_1 pstI --renz_2 ecoRI -E phred33 -r -c -q

## step 2. ustacks, generate unique stcks
## get sample names
cut -f 2 ddRad_barcodes_for-stacks.txt >sample_names.txt 

names=`cat sample_names.txt`

counter=0
for name in $names;
do
  echo "counter: $counter"
  ## stacks treat pairend reads as single loci, so let's merge fastq files for each sample
  cat $out_folder/${name}*.[12].fq.gz >$out_folder/${name}.fq.gz
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

## load stacks
module load Stacks

date
echo "Job started"

## create a folder for stacks output
out_folder="stacks_out"
mkdir $out_folder

## format the barcode file
perl -ne 's/(\S+)\s+(\S+)/$2\t$1/; print $_' ddRad_barcodes.txt  >ddRad_barcodes_for-stacks.txt 
head ddRad_barcodes.txt  ddRad_barcodes_for-stacks.txt

## step 1.  process_radtags
### remove low quality and adapter-containning reads
process_radtags -p fastq -P -i gzfastq -b ddRad_barcodes_for-stacks.txt -o $out_folder --renz_1 pstI --renz_2 ecoRI -E phred33 -r -c -q

## step 2. ustacks, generate unique stcks
## get sample names
cut -f 2 ddRad_barcodes_for-stacks.txt >sample_names.txt 

names=`cat sample_names.txt`

counter=0
for name in $names;
do
  echo "counter: $counter"
  ## stacks treat pairend reads as single loci, so let's merge fastq files for each sample
  cat $out_folder/${name}*.[12].fq.gz >$out_folder/${name}.fq.gz
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

echo "All done!"
date
</pre>




