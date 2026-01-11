#### 🧬 This is my repository for keeping track of my bioinformatics workflow for my conseravtion genomics projects. 🧬

Please see the [Wiki page](https://github.com/18liedan/genomics_memo/wiki) for various random notes and tips/tricks that may help you.

Scripts for general bioinformatic workflows are as follows:

📖 [00_refprep.md](./00_refprep.md) describes how I prepared my reference genomes.
- [get_masked_regions.py](./get_masked_regions.py) is a custom python script used for extracting coordinate of masked sites in the reference genome.

🗺️ [01_trimming-mapping.sh](./01_trimming-mapping.sh) is a script to process all raw genomes

🔡 [02_haplotypecalling.sh](./02_haplotypecalling.sh) is a script to call, merge, genotype, and filter reliable variants

⚖️ 03_bqsr.sh is a script to perform BQSR on all samples

🔡 04_finalhaplotypecalling.sh is a script to call, merge, genotype, and filter final variants for downstream analyses

📏 05_angsdhet.sh is a script to run ANGSD to calculate heterozygosity from sfs files per sample

📏 06_roh.sh is a script to run PLINK to calculate runs of homozygosity for all samples
- [roh_plot.R](./roh_plot.R) is a custom R script used to plot ROH results.

📍 07_pca.sh is a script to run PCA on all samples, with options to subsample different datasets

🗂️ 08_admixture.sh is a script to run ADMIXTURE on all samples, with options to subsample different datasets

🕰️ 09_psmc.sh is a script to run PSMC on all samples

🕰️ 10_smcpp.sh is a script to run SMC++ on all samples, with bootstrapping
- the following python script was used for bootstrapping: https://github.com/popgenmethods/smcpp/files/1182555/bootstrap_smcpp.zip
- [smcpp_makecsv.py](./smcpp_makecsv.py) is a custom python script used for exporting results to a csv.
- [smcpp_plot.R](./smcpp_plot.R) is a custom R script used to plot SMC++ results.

For most of my processes, I use a species prefix that can be changed at the top of each script, so it should be fairly easy to implement the pipelines as long as you adhere to the naming system. Please feel free to download and adjust the scripts for your own study systems! 🦅
