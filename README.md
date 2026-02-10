# PacBio_depth_analysis

## Introduction

This repository contains the scripts and results used to analyze how sequencing depth affects the detection power in PacBio full-length 16S rRNA libraries. Raw reads will be deposited in public archives once the related work is published.

## Dataset

The dataset consists of 30 samples:
- 27 samples are hindguts of honeybees colonized with a synthetic community of 25 bacterial species
- 3 samples (HF, Mix_15_2, Mix_26) are the mixed cultures of these bacteria used to inoculate the bees.

DNA was extracted using a combination of mechanical and chemical lysis, then purified with magnetic beads. Libraries were prepared using a Kinnex protocol and were not normalized based on amount of amplicon. The full length 16S rRNA amplicons were then sequenced at the GTF of UNIL on a Revio machine.

## Method

I ran [this pipeline](https://github.com/momsane/16S_PacBio_dada2) without any pre-subsampling of the reads. Because samples consist of a SynCom and the number and sequence of all 16S rRNA copies of the members are known, I could infer genome equivalents from ASV counts. Therefore, plots are based on genome equivalents instead of reads. I then used the strain quantification matrix to run a simulation of the abundance of each strain at different depths, to see what would be the minimum depth to detect most taxa. The script used to perform this analysis and plot the results is `/scripts/compare_detection.R`.

## Results

- **Raw reads to final ASV table read retention**: on average, using maxEE=4, 73% of raw reads make it through.
- **Genome equivalents to reads relationship**: it seems that genome_depth=0.276*reads_after_tax. This was expected since most bacteria in the SynCom have 3-4 16S rRNA copies. This means that for instance with 20k raw reads, you would get about 4k genome equivalents at the end.
- **Detection of most abundant taxa**: 20k raw reads or 4k genome equivalents (first vertical dashed line in the curves) seems to be enough to detect the most abundant taxa.
- **Rare/low abundance taxa**: it also seems that some strains colonize at very low abundance, *i.e.* Commensalibacter, Bartonella, some Gilliamella and some Bifidobacteria. So if detecting those is needed, a lot more reads are required. However this will come with several downsides:
    - Having enough reads in the beginning: to get say 20k genome equivalents at the end you need at least 100k raw reads for each of your samples. You might end up having to discard many samples that don't meet this requirement
    - Computation then becomes much heavier, in particular for the denoising step which is the main bottleneck. For example, with 150 samples starting at around 60k reads each, denoising on 4 CPUs requires 64G of RAM and runs for 12h.
    - Even if you have those 100k reads, those taxa would still have <10 reads, so your estimation of relative abundance will be very noisy.

## Conclusion

Based on these results, I would recommend subsampling raw reads at 60k to keep a good balance between detection limit and computation costs. Of course this should be adapted to your research question and dataset. If half of your samples have <30k raw reads, you will anyways rarefy the other ones at the end to normalize depth, so it is maybe not needed to process so many reads in the first place.