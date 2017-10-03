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
## Script 3. Short read alignment for variant calling
```
trimmomatic-0.35.jar PE  ILLUMINACLIP:Trimmomatic-0.35/adapters/TruSeq3-PE.fa:2:30:15 TRAILING:30 MINLEN:50
bfc -s 367m -t16 -k 55 
bwa index
bwa mem -t 7 -R 
samblaster -i -o
```
## Script 4. Freebayses variant calling, genotyping and filtering commands
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
## Script 5. GATK variant calling, joint genotyping and filtering commands
```
-T HaplotypeCaller -R -I --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o
-T BaseRecalibrator -R -I -knownSites -o # for known sites Freebayes results were used
-T BaseRecalibrator -R -I -knownSites -BQSR -o # second recallibration round
-T AnalyzeCovariates -R -before -after -plots
-T PrintReads -R -I-BQSR -o
-T HaplotypeCaller -R -I --emitRefConfidence -o
-T GenotypeGVCFs -R -V -V -o                                                                    
-T SelectVariants -R -V  -selectType SNP -o # SNP extraction
-T VariantFiltration -R - V -filter "QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -filterName -o # SNP filtration
-T SelectVariants -R -V -selectType INDEL -o # INDEL extraction
-T VariantFiltration -R -V -filter "QD < 2.0 || FS > 200.0 || SOR > 10.0 || ReadPosRankSum < -20.0" -filterName -o # INDEL filtration


```
## Script 4. Transposable element genome annotation
```
BuildDatabase -name 
RepeatModeler -pa 16 -database 
RepeatMasker -pa 22 -s -no_is -nolow -norna -lib -xm -ace -u -gff -excln -dir 
```
## Script 5. Genome structural annotation (main approach)
```
RepeatMasker -pa 12 -species viridiplantae -xsmall -gff -dir 
CpGcluster.pl CpG 50 1e-5 
barrnap --kingdom euk --threads 8 --reject 0.8 
barrnap --kingdom mito --threads 8 --reject 0.8 
aragorn -v -t -i -l -w -o 
cmscan --cpu 8 -Z 684.576320 --cut_ga --rfam --nohmmonly --tblout --fmt 2 --clanin 
makeblastdb -in  -out  -dbtype nucl -title 
blastn -query  -db -task  -num_threads 8 -evalue 1e-3 -max_target_seqs 1 -outfmt "6 qseqid sseqid sstart send pident qcovhsp evalue" -out



# Kontrola jakości krótkich odczytów RNA-Seq (usuwanie adapterów i rejonów o niskiej jakości) dla próbki $sample
bbduk2.sh -Xmx2g threads=8 in=fastq/$sample_R1.fastq.gz in2=fastq/$sample_R2.fastq.gz out=fastq_trim/$sample_R1.fastq.gz out2=fastq_trim/$sample_R2.fastq.gz  qtrim=w trimq=5 rref=illumina.fasta k=23 mink=11 hdist=1 tbo tpe forcetrimleft=10 minlength=50 removeifeitherbad=t overwrite=t stats=status/$sample.bbduk2_adapters.txt 2> status/$sample.bbduk2_trimming.txt

# Utworzenie pliku z sekwencjami rRNA do usunięcia z odczytów
bedtools getfasta -fi data/genome.fasta -s -bed rna_abinitio/barrnap_euk.gff3 > rRNA.fasta
bedtools getfasta -fi data/genome.fasta -s -bed rna_abinitio/barrnap_mito.gff3 >> rRNA.fasta

# Usuwanie sekwencji rRNA z odczytów RNA-Seq dla próbki $sample
bbduk2.sh -Xmx4g threads=8 in=fastq_trim/$sample_R1.fastq.gz in2=fastq_trim/$sample_R2.fastq.gz out=fastq_decon/$sample_R1.fastq.gz out2=fastq_decon/$sample_R2.fastq.gz fref=rRNA.fasta k=31 hdist=1 removeifeitherbad=t overwrite=t stats=status/$sample.bbduk2_contaminants.txt 2> status/$sample.bbduk2_decontamination.txt

# Składanie de novo transkrpytomu dla próbki $sample (dane strand-specific)
Trinity --seqType fq --max_memory 60G --CPU 8 --SS_lib_type RF --left fastq_decon/$sample_R1.fastq.gz --right fastq_decon/$sample_R2.fastq.gz --output trinity/$sample_trinity 
salmon index -p 8 -t trinity/$sample_trinity/Trinity.fasta -i trinity/$sample_trinity/salmon_index 
salmon quant -p 8 -i trinity/$sample_trinity/salmon_index --useVBOpt --seqBias --gcBias --libType ISR -1 fastq_decon/$sample_R1.fastq.gz -2 fastq_decon/$sample_R2.fastq.gz -o trinity/$sample_trinity/salmon_results
rm -rf trinity/$sample_trinity/chrysalis trinity/$sample_trinity/insilico_read_normalization trinity/$sample_trinity/read_partitions trinity/$sample_trinity/salmon_index

# Składanie de novo transkrpytomu dla próbki $sample (dane non strand-specific)
rnaspades.py -t 8 -m 60 --pe1-1 fastq_decon/$sample_R1.fastq.gz --pe1-2 fastq_decon/$sample_R2.fastq.gz -o rnaspades/$sample
rm -rf rnaspades/$sample/corrected
salmon index -p 8 -t rnaspades/$sample/transcripts.fasta -i rnaspades/$sample/salmon_index 
salmon quant -p 8 -i rnaspades/$sample/salmon_index --useVBOpt --seqBias --gcBias --libType IU -1 fastq_decon/$sample_R1.fastq.gz -2 fastq_decon/$sample_R2.fastq.gz -o rnaspades/$sample/salmon_results
rm -rf rnaspades/$sample/salmon_index

#################################################################################################################################################

# Instalacja PASA
wget -O PASA_v2.1.0.tar.gz  https://github.com/PASApipeline/PASApipeline/archive/v2.1.0.tar.gz
tar xfz PASA_v2.1.0.tar.gz
mv PASApipeline-2.1.0/ PASA/
rm -rf PASA_v2.1.0.tar.gz
cd PASA/
make -j 12
cd seqclean/
tar xfz seqclean.tar.gz
mv seqclean seqclean_dir
mv seqclean_dir/* .
rm -rf seqclean_dir
cd ..
cd pasa_conf
cp pasa.CONFIG.template conf.txt
nano conf.txt  #zmodyfikuj plik pod kątem lokalnych ustawień serwera MySQL
cd ../..

# Przykładowa zawartość pliku conf.txt
## MySQL settings:
# server actively running MySQL
MYSQLSERVER=127.0.0.1
# read-write username and password
MYSQL_RW_USER=PASA
MYSQL_RW_PASSWORD=some_password

# Instalacja BUSCO
git clone https://gitlab.com/ezlab/busco.git
cd busco
wget busco.ezlab.org/datasets/embryophyta_odb9.tar.gz
tar xfz embryophyta_odb9.tar.gz
cd ..

# Instalacja TransDecoder
wget -O TransDecoder_v3.0.1.tar.gz https://github.com/TransDecoder/TransDecoder/archive/v3.0.1.tar.gz
tar xfz TransDecoder_v3.0.1.tar.gz
rm -f TransDecoder_v3.0.1.tar.gz
mv TransDecoder-3.0.1 TransDecoder
cd TransDecoder/
make -j 8
cd ..

# Instalacja eggnogMapper
git clone https://github.com/jhcepas/eggnog-mapper.git
eggnog-mapper/download_eggnog_data.py virNOG

# Instalacja InterProScan
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.23-62.0/interproscan-5.23-62.0-64-bit.tar.gz
tar -pxvzf interproscan-5.23-62.0-64-bit.tar.gz
mv interproscan-5.23-62.0 interproscan
cd interproscan/data
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-11.1.tar.gz
tar -pxvzf panther-data-11.1.tar.gz
rm -f panther-data-11.1.tar.gz
cd ../..

#################################################################################################################################################

# Utwórz plik konfiguracyjny do analizy PASA
cp PASA/pasa_conf/pasa.alignAssembly.Template.txt alignAssembly.config
# Przykładowa zawartość
MYSQLDB=pasa_cucumber
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=75
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=95
validate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY=3
validate_alignments_in_db.dbi:--MAX_INTRON_LENGTH=10000
validate_alignments_in_db.dbi:--MIN_INTRON_LENGTH=20
subcluster_builder.dbi:-m=50

# Uruchom adnotację genomu programem PASA
# transcripts_all.fasta zawiera wszystkie złożone de novo transkrypty o poziomie ekspresji przynajmniej 1 TPM (liczone programem Salmon)
PASA/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g data/genome.fasta -t transcripts_all.fasta --stringent_alignment_overlap 30.0 --INVALIDATE_SINGLE_EXON_ESTS --ALIGNERS gmap --CPU 8

# PASA może się zawieszać około zadania 14, w takim wypadku:
rm -rf assemblies/
PASA/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -s 14 -R -g data/genome.fasta -t transcripts_all.fasta --stringent_alignment_overlap 30.0 --INVALIDATE_SINGLE_EXON_ESTS --ALT_SPLICE --ALIGNERS gmap --CPU 1

# Analiza adnotacji programem BUSCO
busco/BUSCO.py --in pasa_cucumber.assemblies.fasta --out pasa_cucumber --lineage busco/embryophyta_odb9 --mode transcriptome --cpu 8 --force

# Wykrywanie rejonów kodujących transkryptów
# pasa.gtf zawiera adnotacje z pliku pasa_cucumber.pasa_assemblies.gtf z poprawionymi identyfikatorami genów/transkryptów
# pasa.gtmap.txt zawiera dwie kolumny (identyfikator genu, identyfikator transkryptu) oddzielone znakiem tabulacji
TransDecoder/util/cufflinks_gtf_genome_to_cdna_fasta.pl pasa.gtf data/genome.fasta > pasa.fasta
TransDecoder/util/cufflinks_gtf_to_alignment_gff3.pl pasa.gtf > pasa.gff3
wget -O Pfam-A.hmm.gz ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -out uniprot_sprot -dbtype prot -title "Swiss-Prot"
TransDecoder/TransDecoder.LongOrfs -t pasa.fasta --gene_trans_map pasa.gtmap.txt -m 100 -S
hmmscan --cpu 8 --domtblout pfam.domtblout Pfam-A.hmm pasa.fasta.transdecoder_dir/longest_orfs.pep
blastp -query pasa.fasta.transdecoder_dir/longest_orfs.pep  -db uniprot_sprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 8 > blastp.outfmt6
TransDecoder/TransDecoder.Predict -t pasa.fasta --single_best_orf --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --cpu 8
rm -f uniprot_sprot.* Pfam-A.*
TransDecoder/util/cdna_alignment_orf_to_genome_orf.pl pasa.fasta.transdecoder.gff3 pasa.gff3 pasa.fasta > pasa.fasta.transdecoder.genome.gff3

# Adnotacja funkcjonalna programem eggnog-mapper
# proteins.fasta zawiera sekwencje białek z pliku pasa.fasta.transdecoder.pep z poprawionymi identyfikatorami 
eggnog-mapper/emapper.py -i proteins.fasta --output_dir eggnog --output pasa -d virNOG --usemem --cpu 8

# Adnotacja funkcjonalna programem InterProScan
# Plik proteins.fasta ajlepiej rozbić na pliki proteins_X.fasta (X=1,2,3...) zawierających po 1000 sekwencji
interproscan/interproscan.sh -i proteins_X.fasta -b interpro/pasa_X -f TSV,GFF3,HTML,SVG -iprlookup -goterms --pathways --disable-precalc -appl Pfam,TIGRFAM,SFLD,ProDom,Hamap,SMART,CDD,ProSiteProfiles,ProSitePatterns,SUPERFAMILY,PRINTS,PANTHER,PIRSF,Coils,MobiDBLite

# Połącz wyniki poszczególnych uruchomień programu InterProScan
# Pliki tekstowe
cat interpro/*.tsv | sort -k1,1 -n > interpro_results.xls
# Raporty SVG
mkdir -p interpro_svg
for FILENAME in interpro/*.svg.tar.gz
do
tar xfz $FILENAME -C interpro_svg --overwrite
done
zip -r interpro_svg.zip interpro_svg
rm -rf interpro_svg
# Raporty HTML
mkdir -p interpro_html
for FILENAME in interpro/*.html.tar.gz
do
tar xfz $FILENAME -C interpro_html --overwrite
done
cd summary
zip -r interpro_html.zip interpro_html
rm -rf interpro_html
```

## SNV calling commands and parameters
