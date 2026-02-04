#### 🧬 This is my repository for keeping track of my bioinformatics workflows. 🧬

Please see the [Wiki page](https://github.com/18liedan/genomics_memo/wiki) for various notes and basic tips/tricks for beginners like myself.

I have uploaded scripts for general bioinformatic workflows (they use loops so you can process multiple samples and populations at once).
Please note not all scripts are ready yet, and the uploaded ones may be updated/refined occasionally too.
They go in the following order, and in most cases the inputs for one script depend on the outputs of an earlier script.

📖 [00_refprep.md](./00_refprep.md) describes how I prepared my reference genomes.
- [get_masked_regions.py](./get_masked_regions.py) is a custom python script used for extracting coordinates of masked sites in the reference genome. This should be placed in the same directory as your reference genome.

🗺️ [01_trimming-mapping.sh](./01_trimming-mapping.sh) is a script to process all raw genomes.

🔡 [02_haplotypecalling.sh](./02_haplotypecalling.sh) is a script to call, merge, genotype, and filter reliable variants.

⚖️ [03_bqsr.sh](./03_bqsr.sh) is a script to perform BQSR on all samples.

🔡 [04_finalhaplotypecalling.sh](./04_finalhaplotypecalling.sh) is a script to call, merge, genotype, and filter final variants.
Further manual filtering, tailored for your species/population/downstream analyses, is recommended.

🪢 [05_angsdhet.sh](05_angsdhet.sh) is a script to run ANGSD to calculate heterozygosity from sfs files per sample.

🪢 06_angsdstats.sh is a script to run ANGSD to calculate nucleotide diversity, Tajima's D, etc using sfs.

🪢 [07_roh.sh](07_roh.sh) is a script to run PLINK to calculate runs of homozygosity for all samples.
- [roh_plot.R](./roh_plot.R) is a custom R script used to plot ROH results.

📍 [08_pca.md](./08_pca.md) is a series of command line to manually run PCA on samples.

📍 09_admixture.sh is a script to run ADMIXTURE on all samples, with options to subsample different datasets.

📍 [10_angsdfst.sh](10_angsdfst.sh) is a script to run ANGSD to calculate Fst between designated datasets.

🕰️ 11_psmc.sh is a script to run PSMC on all samples.

🕰️ [12_smcpp.sh](12_smcpp.sh) is a script to run SMC++ on all samples, with bootstrapping. I run the software through singularity because I don't have admin priviledges in my cluster, but docker works too if you can use sudo and install necessary files (i.e. you have root permission).
- the following python script was used for bootstrapping: https://github.com/popgenmethods/smcpp/files/1182555/bootstrap_smcpp.zip
- [smcpp_makecsv.py](./smcpp_makecsv.py) is a custom python script used for exporting results to a csv.
- [smcpp_plot.R](./smcpp_plot.R) is a custom R script used to plot SMC++ results.

For most of my processes, I use a species prefix that can be changed at the top of each script, so it should be fairly easy to implement the pipelines as long as you adhere to the naming system. Please feel free to download and adjust the scripts for your own study systems! 🦅 Let me know if you come across any issues.
