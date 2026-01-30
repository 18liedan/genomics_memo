The following commands are the ones I used manually to prepare my reference genomes. Theoretically, I can make this into a bash script but I purposely kept it as manual commands because each genome differs in nature. As examples, I will show what I did with my golden eagle and mountain hawk-eagle genomes. Most steps assume you are in the reference genome directory, so please customize accordingly for your own use.

### Golden Eagle
#### Download reference genome from NCBI
- bAquChr1.4 (Scottish individual): https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_900496995.4/
- I used the RefSeq (GCF_900496995.4) because it had already removed the mitochondria.
- It is a chromosome assembly, with 26 autosomal chromosomes, and W & Z chromosomes, as well as several unplaced scaffolds
- After download, check the md5sum to make sure it is not corrupted.
- To make it easier to use, I renamed it to ge_ref.fna

#### Check if it is already masked or not - this genome was already soft-masked!
`head -100 ge_ref.fna`

#### Sort contigs by length
`seqkit sort -l -r ge_ref.fna > ge_ref_sort.fa`

#### Remove short sequences (< 100 kb)
`seqkit seq -m 100000 ge_ref_sort.fa > ge_ref_sort_lr.fa`

#### Make list of new chromosome labels, excluding sex chromosomes, but including unplaced scaffolds which were longer than 100 kb, which I will rename as chr27-chr47.
```
cat > rename.tsv <<EOF
NC_044004.1	chr1
NC_044005.1	chr2
NC_044006.1	chr3
NC_044007.1	chr4
NC_044008.1	chr5
NC_044009.1	chr6
NC_044010.1	chr7
NC_044011.1	chr8
NC_044012.1	chr9
NC_044013.1	chr10
NC_044014.1	chr11
NC_044015.1	chr12
NC_044016.1	chr13
NC_044017.1	chr14
NC_044018.1	chr15
NC_044019.1	chr16
NC_044020.1	chr17
NC_044021.1	chr18
NC_044022.1	chr19
NC_044023.1	chr20
NC_044024.1	chr21
NC_044025.1	chr22
NC_044026.1	chr23
NC_044027.1	chr24
NC_044028.1	chr25
NC_044029.1	chr26
NW_024470381.1	chr27
NW_024470382.1	chr28
NW_024470383.1	chr29
NW_024470384.1	chr30
NW_024470385.1	chr31
NW_024470386.1	chr32
NW_024470387.1	chr33
NW_024470388.1	chr34
NW_024470389.1	chr35
NW_024470390.1	chr36
NW_024470391.1	chr37
NW_024470392.1	chr38
NW_024470393.1	chr39
NW_024470394.1	chr40
NW_024470395.1	chr41
NW_024470396.1	chr42
NW_024470397.1	chr43
NW_024470398.1	chr44
NW_024470399.1	chr45
NW_024470400.1	chr46
NW_024470401.1	chr47
EOF
```

#### Extract current IDs to a temporary list for filtering
`cut -f1 rename.tsv > keep_ids.txt`

#### Run seqkit to rename chromosomes
```
seqkit grep -f keep_ids.txt ge_ref_sort_lr.fa | \
	seqkit replace -p "^(\S+)" -r "{kv}" -k rename.tsv > ge_ref_softmasked_auto.fa
```

#### Output and check names of chromosomes
`seqkit seq -n ge_ref_softmasked_auto.fa > ge_autochr.txt`

#### Extract soft-masked sites using a [python script](./get_masked_regions.py) for downstream analyses where you want to remove repeated sites
`python3 get_masked_regions.py ge_ref_softmasked_auto.fa > ge_ref_masked_regions.bed`
gatk IndexFeatureFile -I ge_ref_masked_regions.bed
bgzip ge_ref_masked_regions.bed && tabix -p bed ge_ref_masked_regions.bed.gz

#### Based on these masked sites, make a NON-masked region sites file (e.g.to be included in analyses) - you will need this in ANGSD
```
samtools faidx ge_ref_softmasked_auto.fa

# 1. Generate genome sizes (same as before)
cut -f1,2 ge_ref_softmasked_auto.fa.fai > genome.sizes

# 2. Sort, Merge, Complement, and Convert to 1-based in one go
bedtools sort -g genome.sizes -i ge_ref_masked_regions.bed | \
bedtools merge -i - | \
bedtools complement -i - -g genome.sizes | \
awk 'BEGIN{FS="\t"; OFS="\t"} {$2=$2+1; print $1, $2, $3}' > ge_nonmasked_regions.txt
angsd sites index ge_nonmasked_regions.txt
```
#### Index final reference genome
```
bwa-mem2 index ge_ref_softmasked_auto.fa
samtools faidx ge_ref_softmasked_auto.fa
gatk CreateSequenceDictionary -R ge_ref_softmasked_auto.fa -O ge_ref_softmasked_auto.dict
```

### Mountain Hawk-Eagle
#### Not-yet-published genome
- It is a scaffold level assembly with unlabelled W & Z scaffolds.
- As with any other genome, check the md5sum to make sure it is not corrupted.
- To make it easier to use, I renamed it to mhe_ref.fna

#### Sort contigs by length
`seqkit sort -l -r mhe_ref.fna > mhe_ref_sort.fa`

#### Remove short sequences (< 100 kb)
`seqkit seq -m 100000 mhe_ref_sort.fa > mhe_ref_sort_lr.fa`

#### Relabel scaffold number and check stats
```
seqkit replace -p .+ -r "Scaffold{nr}" mhe_ref_sort_lr.fa > mhe_ref_sort_lr_scaff.fa
seqkit fx2tab -n -l -g -i mhe_ref_softmasked.fa > mhe_ref_stats.txt
```

#### Run RepeatMasker
```
RepeatMasker mhe_ref_sort_lr_scaff.fa -species "nisaetus nipalensis" -xsmall -no_is -pa 32 -e hmmer
#-no_is: skip bacterial insertion element check
#-pa: threads
#-e: search engine
#-xsmall: softmasking
#result was: mhe_ref_sort_lr_scaff.fa.masked
#renamed to: mhe_ref_softmasked.fa
```

#### Filter out sex chromosomes
Make local blast database to find scaffs with sex-linked genes by blasting golden eagle w/z chr sequences.
(chrW_ge_head.fa & chrZ_ge_head.fa taken from bAquChr1.4 genome, chd_bAquChr1.4.fna from NCBI database)
```
makeblastdb -in mhe_ref_softmasked.fa -out mhedb -dbtype nucl -parse_seqids
blastn -db mhedb -query chrW_ge_head.fa -out blastedseq_W_out.txt
blastn -db mhedb -query chrZ_ge_head.fa -out blastedseq_Z_out.txt
```
Scaffold1 and Scaffold42 had high identity and low E value.

#### Remove sex chromosomes from reference genome
```
seqkit grep -v -p "Scaffold1" mhe_ref_softmasked.fa > mhe_ref_softmasked_autotemp.fa
seqkit grep -v -p "Scaffold42" mhe_ref_softmasked_autotemp.fa >  mhe_ref_softmasked_auto.fa
rm  mhe_ref_softmasked_autotemp.fa #remove temporary file
seqkit seq -n gfb_ref_softmasked_auto.fa  #check if it is really removed
```
#### Extract soft-masked sites using a [python script](./get_masked_regions.py) for downstream analyses where you want to remove repeated sites
`python3 get_masked_regions.py mh_ref_softmasked_auto.fa > mhe_ref_masked_regions.bed`
gatk IndexFeatureFile -I mhe_ref_masked_regions.bed
bgzip mhe_ref_masked_regions.bed && tabix -p bed mhe_ref_masked_regions.bed.gz

#### Based on these masked sites, make a NON-masked region sites file (e.g.to be included in analyses) - you will need this in ANGSD
```
samtools faidx mhe_ref_softmasked_auto.fa

# 1. Generate genome sizes
cut -f1,2 mhe_ref_softmasked_auto.fa.fai > genome.sizes

# 2. Sort, Merge, Complement, and Convert to 1-based in one go
bedtools sort -g genome.sizes -i mhe_ref_masked_regions.bed | \
bedtools merge -i - | \
bedtools complement -i - -g genome.sizes | \
awk 'BEGIN{FS="\t"; OFS="\t"} {$2=$2+1; print $1, $2, $3}' > mhe_nonmasked_regions.txt
angsd sites index mhe_nonmasked_regions.txt
```
#### Index final reference genome
```
bwa-mem2 index ge_ref_softmasked_auto.fa
samtools faidx ge_ref_softmasked_auto.fa
gatk CreateSequenceDictionary -R ge_ref_softmasked_auto.fa -O ge_ref_softmasked_auto.dict
```
#### Index final reference genome
```
bwa-mem2 index mhe_ref_softmasked_auto.fa
samtools faidx mhe_ref_softmasked_auto.fa
gatk CreateSequenceDictionary -R mhe_ref_softmasked_auto.fa -O mhe_ref_softmasked_auto.dict
```
#### Extract soft-masked sites using a [python script](./get_masked_regions.py) for downstream analyses where you want to remove repeated sites
`python3 get_masked_regions.py mhe_ref_softmasked_auto.fa > mhe_ref_masked_regions.bed`
