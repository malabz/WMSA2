## WMSA 2: A multiple DNA/RNA sequence alignment tool implemented with accurate progressive mode and a fast win-win mode combining the center star and progressive strategies
[website](http://lab.malab.cn/soft/WMSA2) 

[CLUSTAR](https://github.com/Guuhua/CLUSTAR)

### Introduction

Multiple sequence alignment is widely used for sequence analysis, such as identifying important sites and phylogenetic analysis. Traditional methods, such as progressive alignment, are time-consuming. To address this issue, we introduce StarTree, a novel method to fast construct a guide tree by combining sequence clustering and hierarchical clustering. Furthermore, we develop a new heuristic similar region detection algorithm using the FM-index and apply the k-banded dynamic program to the profile alignment. We also introduce a win-win alignment algorithm that applies the central star strategy within the clusters to fast the alignment process, then uses the progressive strategy to align the central-aligned profiles, guaranteeing the final alignment's accuracy. We present WMSA 2 based on these improvements and compare the speed and accuracy with other popular methods. The results show that the guide tree made by the StarTree clustering method can lead to better accuracy than that of PartTree while consuming less time and memory than that of UPGMA and mBed methods on datasets with thousands of sequences. During the alignment of simulated data sets, WMSA 2 can consume less time and memory while ranking at the top of Q and TC scores. The WMSA 2 is still better at the time, and memory efficiency on the real datasets and ranks at the top on the average sum of pairs score. For the alignment of 1 million SARS-CoV-2 genomes, the win-win mode of WMSA 2 significantly decreased the consumption time than the former version. The source code and data are available at https://github.com/malabz/WMSA2.

### Installation

#### For OSX/Linux/Windows

1. Download and install JDK (version >= 11) for different systems from [here](https://www.oracle.com/java/technologies/downloads/#java11).
2. Download WMSA 2 from [relseases](https://github.com/malabz/WMSA2/releases).

#### Usages

```bash
java -jar WMSA2.jar [-m] mode [-s] similarity [-i] path [-o] path
```

```
necessary arguments: 
    -i  Input file path (nucleotide sequences in fasta format)
    -o  Output file path

optional arguments:
  -m  align option (default mode: Pro)
       1. Pro        progressive alignment with StarTree method
       2. Win        combine central and progressive alignment with CluStar method
       3. StarTree   only output the guidetree
  -t  tow guide tree option for MSA (default: StarTree)
       1. StarTree
       2. UPGMA
  -s  the similarity of the cluster (used in Win mode. default: 0.9)
```

#### Example

1. Download [data](https://github.com/malabz/WMSA2#data) and uncompress the data.
2. Align the data with WMSA 2.

Use Pro mode to align the data.

```bash
java -jar WMSA2.jar -i test.fasta -o test_algined_WMSA2.fasta
```

Use Win mode to align the data(the default similarity value is 0.9).

```bash
java -jar WMSA2.jar -m Win -s 0.95 -i test.fasta -o test_algined_WMSA2.fasta
```

### Download

#### Data

**Hierarchical tree simulated datasets**

> Ref: [HAlign 3](https://doi.org/10.1093/molbev/msac166): Fast Multiple Alignment of Ultra-Large Numbers of Similar DNA/RNA Sequences.

Download: ([more information](https://github.com/malabz/HAlign-3/tree/main/dataset#hierarchical-tree-simulated-datasets))

- [sars_cov_2_like_diff_similarity.tar.xz](http://lab.malab.cn/~tfr/HAlign3_testdata/sars_cov_2_like_diff_similarity.tar.xz) 
- [sars_cov_2_like_diff_treelength.tar.xz](http://lab.malab.cn/~tfr/HAlign3_testdata/sars_cov_2_like_diff_treelength.tar.xz)  
- [mt_like_diff_similarity.tar.xz](http://lab.malab.cn/~tfr/HAlign3_testdata/mt_like_diff_similarity.tar.xz)  
- [mt_like_diff_treelength.tar.xz](http://lab.malab.cn/~tfr/HAlign3_testdata/mt_like_diff_treelength.tar.xz)

**simluation RNA**

Download: ([more information](https://kim.bio.upenn.edu/software/csd.shtml)) ([all](http://lab.malab.cn/soft/WMSA2/data/RNA-all.zip))
   - [RNA-255](http://lab.malab.cn/soft/WMSA2/data/RNA-255.zip)
   - [RNA-511](http://lab.malab.cn/soft/WMSA2/data/RNA-511.zip)
   - [RNA-1023](http://lab.malab.cn/soft/WMSA2/data/RNA-1023.zip)
   - [RNA-2047](http://lab.malab.cn/soft/WMSA2/data/RNA-2047.zip)
   - [RNA-4095](http://lab.malab.cn/soft/WMSA2/data/RNA-4095.zip)
   - [RNA-8191](http://lab.malab.cn/soft/WMSA2/data/RNA-8191.zip)

**Ancient human partial mitochondrial DNA**

> Ref: Ancient DNA Reveals Key Stages in the Formation of Central European Mitochondrial Genetic Diversity. [DOI: 10.1126/science.1241844](https://www.science.org/doi/10.1126/science.1241844)

Download: 

- [mt.zip](http://lab.malab.cn/soft/WMSA2/data/mt.zip)

**SARS-CoV-2 genome**

> Ref: [HAlign 3](https://doi.org/10.1093/molbev/msac166): Fast Multiple Alignment of Ultra-Large Numbers of Similar DNA/RNA Sequences.

Download:

- [SARS-CoV-2_1020.zip](http://lab.malab.cn/soft/WMSA2/data/SARS-CoV-2_1020.zip)
- [SARS-CoV-2_1M.tar.xz](http://lab.malab.cn/~tfr/HAlign3_testdata/sars_cov_2_1Mseq.tar.xz) (more information see this [link](https://github.com/malabz/HAlign-3/tree/main/dataset#respiratory-syndrome-coronavirus-2-genomes))

#### Related Tools

**Score scripts**

- [SP score](https://github.com/malabz/MSATOOLS/tree/main/SPscore)
- [Q score](https://www.drive5.com/qscore/)

**MSA tools**

- [WMSA](https://github.com/malabz/WMSA)
- [HAlign3](https://github.com/malabz/HAlign-3/)
- [FMAlign](https://github.com/iliuh/FMAlign)
- [MAFFT](http://mafft.cbrc.jp/alignment/software/)
- [Kalign3](https://github.com/TimoLassmann/kalign)
- [Clustal-Omega](http://www.clustal.org/omega/)
- [MUSCLE](https://drive5.com/muscle/downloads_v3.htm)

### License

> WMSA 2 is a free software, License under [MIT](https://github.com/malabz/WMSA2/blob/main/LICENSE).

### Citation

If you use this software, please cite:

Chen J, Chao J, Liu H, Yang F, Zou Q, Tang F. WMSA 2: a multiple DNA/RNA sequence alignment tool implemented with accurate progressive mode and a fast win-win mode combining the center star and progressive strategies. Brief Bioinform. 2023 May 17:bbad190. doi: 10.1093/bib/bbad190. Epub ahead of print. PMID: 37200156.

### Contact us

The software tools are developed and maintained by [ZOU's lab](http://lab.malab.cn/~zq/en/index.html). If you find any bug, welcome to contact us on the [issues page](https://github.com/malabz/WMSA2/issues). More tools and infomation can visit our [github](https://github.com/malabz).
