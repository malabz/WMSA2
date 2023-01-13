## WMSA 2: A win-win multiple sequence alignment method combining central and progressive alignment strategies for the fast alignment of big data sets
[website](http://lab.malab.cn/soft/WMSA2) 

[CLUSTAR](https://github.com/Guuhua/CLUSTAR)

### Introduction

Multiple sequence alignment is essential for various sequence analyses, including identifying important sites and phylogenetic analysis. Currently, most tools adopt the progressive alignment strategy, which first needs to create a guide tree and perform profile-profile alignment; both steps have high time complexity. Thus, we introduce StarTree, a novel method to fast construct a guide tree by combining heuristic and hierarchical clustering. Furthermore, to overcome the high time consumption of profile-profile alignment, we develop a new heuristic similar region detection algorithm using the FM-index and apply the k-banded dynamic program to the profile alignments. Benefiting from the StarTree method, we also introduce a win-win alignment algorithm that combines the central and progressive alignment strategies. Apply the central strategy within the clusters to fast the alignment process, then use the progressive strategy to align the central-aligned profiles guaranteeing the accuracy of final alignment. Here, we present WMSA 2 based on these improvements and compare the speed and accuracy of WMSA 2 with other popular methods. The results show that the guide tree made by the StarTree clustering method can lead to better accuracy than that of PartTree while consuming less time and memory than that of UPGMA and mBed methods. During the alignment of simulated data sets, WMSA 2 can consume less time and memory while ranking at the top of Q and TC scores. The WMSA 2 is still better at the time and memory efficiency on the real data sets, with slightly less accuracy based on the average SP score. For the alignment of 1 million SARS-CoV-2 genomes, the win-win mode of WMSA 2 is better than the former version in efficiency and accuracy, which was mainly developed on the alignment tool MAFFT and CD-HIT clustering methods. Source code and data: https://github.com/malabz/WMSA2.

### Installation

#### For OSX/Linux/Windows

1. Download and install JDK (version >= 11) for different systems from [here](https://www.oracle.com/java/technologies/downloads/#java11).
2. Download WMSA2 from [relseases](https://github.com/malabz/WMSA2/releases).

#### Usages

```bash
java -jar WMSA2.jar [-m] mode [-s] similarity [-i] path [-o] path
```

```
necessary arguments: 
    -i  Input file path (nucleotide sequences in fasta format)
    -o  Output file path

optional arguments: 
    -m  three align option (default mode: Pro)
        1. Pro   more accurate but slower
        2. Win   less accurate but faster
    -s  the similarity of the cluster (used in Win mode. default: 0.9)
```

#### Example

1. Download [data](http://lab.malab.cn/soft/WMSA2/index.html#load) and uncompress the data.
2. Align the data with WMSA2.

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

**AmtDB**

> Ref: [AmtDB](https://amtdb.org/): a database of ancient human mitochondrial genomes.

Download: 

- [amtdb.zip](http://lab.malab.cn/soft/WMSA2/data/amtdb.zip)

**SARS-CoV-2 genome**

> Ref: [HAlign 3](https://doi.org/10.1093/molbev/msac166): Fast Multiple Alignment of Ultra-Large Numbers of Similar DNA/RNA Sequences.

Download:

- [SARS-CoV-2_1020.zip](http://lab.malab.cn/soft/WMSA2/data/SARS-CoV-2_1020.zip)
- [SARS-CoV-2_1M.tar.xz](http://lab.malab.cn/~tfr/HAlign3_testdata/sars_cov_2_1Mseq.tar.xz) (more information see this [link](https://github.com/malabz/HAlign-3/tree/main/dataset#respiratory-syndrome-coronavirus-2-genomes))

**more datasets **

Download:

- [site](http://lab.malab.cn/~cjt/MSA/datasets.html)

#### Related Tools

**Score scripts**

- [SP score](https://github.com/malabz/MSATOOLS/tree/main/SPscore)
- [Q score](https://www.drive5.com/qscore/) ([usage](http://lab.malab.cn/~cjt/MSA/measure.html))

**MSA tools**

- [WMSA](https://github.com/malabz/WMSA)
- [HAlign3](https://github.com/malabz/HAlign-3/)
- [FMAlign](https://github.com/iliuh/FMAlign)
- [MAFFT](http://mafft.cbrc.jp/alignment/software/)
- [Kalign3](https://github.com/TimoLassmann/kalign)
- [Clustal-Omega](http://www.clustal.org/omega/)
- [MUSCLE](https://drive5.com/muscle/downloads_v3.htm)

### License

> WMSA2 is a free software, License under [MIT](https://github.com/malabz/WMSA2/blob/main/LICENSE).

### Citation

If you use this software, please cite:

WMSA2: A win-win multiple sequence alignment method combining central and progressive alignment strategies for the fast alignment of big data sets(under review).

### Contact us

The software tools are developed and maintained by [MALAB](http://lab.malab.cn/~cjt/MSA). If you find any bug, welcome to contact us on the [issues page](https://github.com/malabz/WMSA2/issues). More tools and infomation can visit our [website](http://lab.malab.cn/~cjt/MSA) or [github](https://github.com/malabz).
