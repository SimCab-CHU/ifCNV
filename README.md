# ifCNV : a novel isolation-forest-based package to detect copy number variations from various NGS datasets

## Installation

### From Conda (recommended)

1. Install Conda : [documentation here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
1. Create an envrionment : `conda create -n ifcnv`
1. Activate the environment : `conda activate ifcnv`
1. Install the package : `conda install -c conda-forge -c bioconda ifcnv`

### From PyPi

Make sure you have python >= 3.6

```sh
python --version
Python 3.X.X
```

Install pip : [documentation here](https://pip.pypa.io/en/stable/installation/).

```sh
python -m ensurepip --upgrade
```

Install [ifCNV](https://pypi.org/project/ifCNV/)

```sh
pip install ifCNV
```

## Users Guide

### Basic usage

ifCNV stands for **i**solation **f**orest based **C**opy **N**umber **V**ariation detection. 

Its usage is meant to be approachable for entry-level users. 

```sh
ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/
```

**ifCNV** input: 
- the aligned sequences (.bam or .cram files and their associated indexes)
- the genomic coordinates of the region of interest (a .bed file, see 
http://genome.cse.ucsc.edu/FAQ/FAQformat.html#format1 and below for more information).


**ifCNV** ouput:
- An html report

Some **BED file** informations:

The BED (Browser Extensible Data) format is a text file format used to store genomic regions as coordinates and associated annotations. The data are presented in the form of columns separated by spaces or tabs. This format was developed during the Human Genome Project and then adopted by other sequencing projects (source: wikipedia). For targeted sequencing, it is supposed to be provided by the vendor.

As it is of important for **ifCNV** let's see some key points:
- the BED file must have (at least) 4 columns, the three first are the cooridnates of the baited genomic regions and the forth one is a character describing these regions.
- The forth column is the name of the targeted region and it must be set carefuly as it will be used to compute the localization score (see the **Resolution paragraph**)

If you don’t have the capture regions BED file, but you do know which commercial exome capture kit was used to prepare your samples, you might find the file you need in Astra-Zeneca’s reference data repository. Otherwise, try searching the vendor’s website or contacting their customer support for the right file.

### Interpretation of the html report:

<img src="docs/img/output_example_2.png" alt="drawing" width="800"/>

- The 1st column is the run name (defult is ifCNV, but it can be specified with the -r flag).
- The 2nd column is the sample name (by default the table is sorted by this column).
- The 3rd column is the name of the region (here a gene).
- The 4th column is the reads ratio (1<reads_ratio<2 gain, reads_ratio>2 amplification, reads_ratio<1 loss).
- The 5th column is the localization score.

The user can click on each line to visualize a graph of the normalized log ratio for the region of interest where each point represents a target.

<img src="docs/img/output_example_4.png" alt="drawing" width="800"/>


### More specific commands

- if you want ifCNV to stop talking
```sh
ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/ -v ''
```
- if you don't want ifCNV to automatically open the report
```sh
ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/ -a ''
```
- if you don't want to re-create the reads matrix of your run at each
utilisation (ie. to rerun ifCNV with different parameters)
    - First, run ifCNV and write the reads matrix
    - Then, run ifCNV and tell it to take the reads matrix (it will skip its
    creation)

```sh
ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/ -rm /path/to/readsMatrix/file

ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/ -s /path/to/readsMatrix/file
```
- if you want to save the output in a .tsv file (it will be written in the
output directory)
```sh
ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/ -sv True
```
- changing mode from 'fast' (default) to 'extensive', will take a bit longer
but doing so ifCNV will check every sample of the run for altered targets
```sh
ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/ -m 'extensive'
```

### Resolution

The resolution of ifCNV is the smallest region covered by a target. It is set in
the .bed file. The 4th column of the .bed file must be the name of the targeted
region. ifCNV splits this name on the "\_" character and regroups whats on the
left of it as the region of interest.

Examples:
```tsv
1	65300120	65300271	JAK1_E25_STA100319	0	+	65300140	65300249	255,0,0
1	65300220	65300383	JAK1_E25_STA100320	0	+	65300246	65300353	255,0,0
1	65300278	65300427	JAK1_E25_STA100321	0	+	65300296	65300393	255,0,0
1	65301022	65301161	JAK1_E24_STA100322	0	+	65301046	65301143	255,0,0
1	65301065	65301228	JAK1_E24_STA100323	0	+	65301085	65301202	255,0,0
```
Using this .bed file, the resolution will be at the gene level, meaning the five
amplicons will be considered as belonging to the same region (JAK1) in the
calculus of the localisation score.

```tsv
1	65300120	65300271	JAK1-E25_STA100319	0	+	65300140	65300249	255,0,0
1	65300220	65300383	JAK1-E25_STA100320	0	+	65300246	65300353	255,0,0
1	65300278	65300427	JAK1-E25_STA100321	0	+	65300296	65300393	255,0,0
1	65301022	65301161	JAK1-E24_STA100322	0	+	65301046	65301143	255,0,0
1	65301065	65301228	JAK1-E24_STA100323	0	+	65301085	65301202	255,0,0
```

Using this .bed file, the resolution will be at the exon level, meaning the 3
first will be considered as belonging to the same region (JAK1-E25) and the 2
last to another region (JAK1-E24) in the calculus of the localisation score.

This implies a careful consideration to the localisation score threshold (-sT).
Indeed, the localisation score depends on the size of the region of interest.
For example, 3 altered targets on 3 targets of the region of interest will have
a smaller localisation score than 10 altered targets on 10 targets of the region
of interest (see image below).

<img src="docs/img/score_plot.png" alt="drawing" width="450"/>

### Contamination parameters

ifCNV uses 2 Isolation forests, one to detect the outlying samples (considered
as CNV positives) and another to detect the outlying targets. The
_contamination_ is a parameter of the isolation forest that defines the
proportion of outliers in the data set. It is set for both IF as "auto" by
default but can be changed by the user.

Changing the -ct parameter of ifCNV can be useful but a careful consideration
must be taken on the score threshold (-sT). For example, if the user sets the
-ct parameter to a small value (\~ \]0,0.01\]), less targets will be considered
as outliers and so the localisation scores will be lower. On the other hand, if
the user sets the -ct parameter to a high value (\~ \]0.1,0.5\]), more targets
will be considered as outliers and so the localisation scores will be higher.

An example is provided in the paper describing ifCNV using the ICR96 dataset,
and is summarized in the figure below:

<img src="docs/img/figure_5.jpg" alt="drawing" width="450"/>
