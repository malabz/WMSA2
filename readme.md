## WMSA 2: A win-win multiple sequence alignment method combining central and progressive alignment strategies for the fast alignment of big data sets
[website](http://lab.malab.cn/soft/MACF)

### Sequence clustering

[CLUSTAR](https://github.com/Guuhua/CLUSTAR)

### Usage

```shell
usage: java -jar WMSA2.jar [-m] mode [-i] path [-o] path

  necessary arguments: 
    -i  Input file path (nucleotide sequences in fasta format)
    -o  Output file path

  optional arguments: 
    -m  three align option (default mode: Pro)
         1. Pro   more accurate but slower
         2. Win   less accurate but faster
    -s  the similarity of the cluster (used in Win mode. default: 0.9)
```

### Change Log
---

- 11-01-2022, version 1.0
  
  inital version
  

### Dependencies

MACF requires JDK environment (version >= 1.11).

- [Download JDK](https://www.oracle.com/java/technologies/downloads/)

- [Installation Guide](https://docs.oracle.com/en/java/javase/17/install/overview-jdk-installation.html)
