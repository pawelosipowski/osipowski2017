# Script 1. Reference genome sequence assembly and processing commands and parameters used
### PBcR v8.3rc2:
```
PBcR -sensitive -length 500 -partitions 200 -l -s -fastq genomeSize=367000000
```

### PBcR v8.3rc2 specification file content:
```
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
### Pbalign v3.0:
```
pbalign --nproc 22
```
### Quiver v2.1.0: 
```
quiver --referenceFilename= -j22 --annotateGFF -o 
```
### Trimmomatic v0.35
```
trimmomatic-0.35.jar PE ILLUMINACLIP:TruSeq3-PE.fa:2:30:15 TRAILING:30 MINLEN:50
```
### BBTools ecc.sh script v35.82
```
ecc.sh in= in2= out1= out2= hist= histout=
```
### Bowtie2 v2.2.9
```
bowtie2 -x -1 -2 -X 300 -p 7
```
### Samblaster v0.1.24
```
samblaster -i -o 
```
### Pilon v1.20
```
pilon-1.20.jar --genome --bam --output --vcf --chunksize 12000000 --diploid
```
### ScaffMatch v0.9
```
scaffmatch -w ./scfm_greedy/ -c -1 -2 -i 135000,300 -p fr,fr -s 60000,50 -t 1 -g max_weight -l
```
# Script 2. SNV calling commands and parameters
