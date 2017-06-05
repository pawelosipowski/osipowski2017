##Structural Annotation
### This is sequence of operations and linux bash commands to transfer structural annotation from one version of the reference genome to another
get promoters annotation (flanks 500n after coding start & 1000n before coding start)
```
bedtools flank -i genes_pasa.bed -g /storage3/B10v2/2_polishing/ncbi_filtered/previous_ecc_corrected/B10v2_nm_sizes.genome -l 1000 -r 0 -s > temp.bed
bedtools flank -i temp.bed -g /storage3/B10v2/2_polishing/ncbi_filtered/previous_ecc_corrected/B10v2_nm_sizes.genome -l 0 -r 500 -s > promoters_pasa.bed
```
get promoters fasta sequences
```
bedtools getfasta -fi /storage3/B10v2/2_polishing/ncbi_filtered/previous_ecc_corrected/B10v2_c_corr.fsa -bed promoters_pasa.bed -s -name -fo promoters_pasa.fasta
```
get exons fasta sequences
```
bedtools getfasta -fi /storage3/B10v2/2_polishing/ncbi_filtered/previous_ecc_corrected/B10v2_c_corr.fsa -bed exons_pasa.bed -s -name -fo exons_pasa.fasta
```
blast annotation on quiver-only_corrected_genome
```
blastn -query promoters_pasa.fasta -db /storage3/B10v2/2_polishing/ncbi_filtered/final_quiveronly_corrected/B10v2_quiver_nm_c_corr.fasta -task blastn -outfmt 6 -perc_identity 95 -best_hit_overhang 0.25 -num_threads 21 > promoters_blasted_temp
```
remove blast dups by first column leaving first position (best score)
```
awk '!a[$1]++' promoters_blasted_temp > promoters_blasted
```
set of file modifications: 1. split strand from name to column 2 by "(" 2. and print table in .bed order 3. swaping stop/start for lower value in 2nd column
```
awk -F"(" '$1=$1' OFS="\t" genes_blasted | awk 'OFS="\t" {gsub(/\)/,"");print $3,$10,$11,$1, ".",$2,"PASA","gene",".",$1}' | awk '$2 > $3 { temp = $3; $3 = $2; $2 = temp } 1' OFS='\t' > genes_new.be
```
