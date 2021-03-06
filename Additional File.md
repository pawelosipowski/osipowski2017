# Additional file. Commands.
## Script 1. Reference genome sequence assembly and processing commands and parameters used
```
PBcR -sensitive -length 500 -partitions 200 -l -s -fastq genomeSize=367000000
```
```
PBcR v8.3rc2 specification file content:
ovlMemory=32
ovlStoreMemory=32000
merylMemory=32000
threads=32
ovlConcurrency=1
cnsConcurrency=32
merylThreads=32
frgCorrThreads=32
frgCorrBatchSize=100000
ovlCorrBatchSize=100000
asmOvlErrorRate=0.1
asmUtgErrorRate=0.06
asmCgwErrorRate=0.1
asmCnsErrorRate=0.1
asmOBT=1
asmObtErrorRate=0.08
asmObtErrorLimit=4.5
batOptions=-RS -NS -CS -SR -R
utgGraphErrorRate=0.05
utgMergeErrorRate=0.05
ovlHashBits=24
ovlHashLoad=0.8
ovlHashBlockLength=300000000
ovlRefBlockLength=0
ovlRefBlockSize=2000000
```
## Script 2. Reference genome correction and quality assessment

```
pbalign --nproc 22
quiver --referenceFilename= -j22 --annotateGFF -o 
trimmomatic-0.35.jar PE ILLUMINACLIP:TruSeq3-PE.fa:2:30:15 TRAILING:30 MINLEN:50
ecc.sh in= in2= out1= out2= hist= histout=
bowtie2 -x -1 -2 -X 300 -p 7
samblaster -i -o 
pilon-1.20.jar --genome --bam --output --vcf --chunksize 12000000 --diploid
run_BUSCO.py -i -o -l -m genome -c 5 -sp arabidopsis --long
```
## Script 3. Yang et al. 2013 primer data alignment
```
bwa mem -k 5 -t 20 -r 0.1 -c 1000000
```
## Script 4. Short read alignment for variant calling
```
trimmomatic-0.35.jar PE  ILLUMINACLIP:Trimmomatic-0.35/adapters/TruSeq3-PE.fa:2:30:15 TRAILING:30 MINLEN:50
bfc -s 367m -t16 -k 55 
bwa index
bwa mem -t 7 -R 
samblaster -i -o
```
## Script 5. Freebayses variant calling, genotyping and filtering commands
```
mdust -w 6 -v 20 -c 
samtools faidx 
freebayes -f -b
vcffilter -f "QUAL > 30" 
lcr_filter.py - our own script
vcffilter -f "DP > 30"
avg_depth=$(cut -f 8 .vcf | awk -F \; '{print $8}'| sed 's/DP=//' | sed '/^$/d' | awk '{a+=$1} END{print a/NR}')
depth_score=$(echo "$avg_depth + sqrt($avg_depth) * 3" | bc -l)
vcffilter -f "DP < $depth_score" 
vcffilter -f "SAF > 5" 
vcffilter -f "SAR > 5" 

```
## Script 6. DeepVariant variant calling shell script template
```
SMPL=
REF=
BAM=
DV_FOLDER=
MODEL_NAME=DeepVariant-inception_v3-0.4.0+cl-174375304.data-wgs_standard

OUTPUT_DIR="${SMPL}_output"
mkdir -p "${OUTPUT_DIR}"
MODEL="${DV_FOLDER}/${MODEL_NAME}/model.ckpt"

LOGDIR="./logs_${SMPL}"
N_SHARDS=6

#mkdir -p "${LOGDIR}"
#time seq 0 $((N_SHARDS-1)) | \
  parallel --eta --halt 2 --joblog "${LOGDIR}/log" --res "${LOGDIR}" \
  python "${DV_FOLDER}/bazel-bin/deepvariant/make_examples.zip" \
    --mode calling \
    --ref "${REF}" \
    --reads "${BAM}" \
    --examples "${OUTPUT_DIR}/examples.tfrecord@${N_SHARDS}.gz" \
    --task {}
CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/call_variants_output.tfrecord.gz"

python "${DV_FOLDER}/bazel-bin/deepvariant/call_variants.zip" \
 --outfile "${CALL_VARIANTS_OUTPUT}" \
 --examples "${OUTPUT_DIR}/examples.tfrecord@${N_SHARDS}.gz" \
 --checkpoint "${MODEL}"

FINAL_OUTPUT_VCF="${OUTPUT_DIR}/output_${SMPL}.vcf.gz"

python  "${DV_FOLDER}/bazel-bin/deepvariant/postprocess_variants.zip" \
  --ref "${REF}" \
  --infile "${CALL_VARIANTS_OUTPUT}" \
  --outfile "${FINAL_OUTPUT_VCF}"
```
## Script 7. GRIDSS SV calling
```
gridss O= ASSEMBLY= REFERENCE_SEQUENCE= I= TMP_DIR= WORKING_DIR= THREADS=7 
```
## Script 8. Transposable element genome annotation
```
BuildDatabase -name 
RepeatModeler -pa 16 -database 
RepeatMasker -pa 22 -s -no_is -nolow -norna -lib -xm -ace -u -gff -excln -dir 
```
## Script 9. Trancriptome-based genome structural annotation (main approach)
```
RepeatMasker -pa 12 -species viridiplantae -xsmall -gff -dir 
CpGcluster.pl CpG 50 1e-5 
barrnap --kingdom euk --threads 8 --reject 0.8 
barrnap --kingdom mito --threads 8 --reject 0.8 
aragorn -v -t -i -l -w -o 
cmscan --cpu 8 -Z 684.576320 --cut_ga --rfam --nohmmonly --tblout --fmt 2 --clanin 
makeblastdb -in  -out  -dbtype nucl -title 
blastn -query  -db -task  -num_threads 8 -evalue 1e-3 -max_target_seqs 1 -outfmt "6 qseqid sseqid sstart send pident qcovhsp evalue" -out
bbduk2.sh -Xmx2g threads=8 in= in2= out= out2= qtrim=w trimq=5 rref= k=23 mink=11 hdist=1 tbo tpe forcetrimleft=10 minlength=50 removeifeitherbad=t overwrite=t stats= 
bedtools getfasta -fi -s -bed
bbduk2.sh -Xmx4g threads=8 in= in2= out= out2= fref= k=31 hdist=1 removeifeitherbad=t overwrite=t stats= 
Trinity --seqType fq --max_memory 60G --CPU 8 --SS_lib_type RF --left  --right --output 
salmon index -p 8 -t -i 
salmon quant -p 8 -i --useVBOpt --seqBias --gcBias --libType ISR -1 -2 -o
rnaspades.py -t 8 -m 60 --pe1-1 --pe1-2 -o
salmon index -p 8 -t -i 
salmon quant -p 8 -i --useVBOpt --seqBias --gcBias --libType IU -1 -2 -o
```
```
PASA configuration file content
MYSQLDB=pasa_cucumber
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=75
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=95
validate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY=3
validate_alignments_in_db.dbi:--MAX_INTRON_LENGTH=10000
validate_alignments_in_db.dbi:--MIN_INTRON_LENGTH=20
subcluster_builder.dbi:-m=50
```
```
Launch_PASA_pipeline.pl -c -C -R -g -t --stringent_alignment_overlap 30.0 --INVALIDATE_SINGLE_EXON_ESTS --ALIGNERS gmap --CPU 8
BUSCO.py --in --out --lineage busco/embryophyta_odb9 --mode transcriptome --cpu 8 --force
cufflinks_gtf_genome_to_cdna_fasta.pl  
cufflinks_gtf_to_alignment_gff3.pl
makeblastdb -in -out -dbtype prot -title
TransDecoder.LongOrfs -t --gene_trans_map -m 100 -S
hmmscan --cpu 8 --domtblout
blastp -query -db -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8
TransDecoder.Predict -t --single_best_orf --retain_pfam_hits --retain_blastp_hits --cpu 8
cdna_alignment_orf_to_genome_orf.pl
```
## Script 10. Ab initio genome structural annotation (complementary approach)
```
tophat2 --num-threads 12 -o
braker.pl --AUGUSTUS_CONFIG_PATH= --AUGUSTUS_BIN_PATH= --AUGUSTUS_SCRIPTS_PATH= --BAMTOOLS_PATH= --GENEMARK_PATH= --cores 20 --species= --gff3 --genome= --bam=
```
