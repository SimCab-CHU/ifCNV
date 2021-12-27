# ifCNV : a novel isolation-forest-based package to detect copy number variations from various NGS datasets.

## Installation

### From PyPi (recommended)

Make sure you have python >= 3.6

```sh
$ python --version
Python 3.X.X
```

Install pip (https://pip.pypa.io/en/stable/installation/).

```sh
$ python -m ensurepip --upgrade
```

Install ifCNV (https://pypi.org/project/ifCNV/)

```sh
$ pip install ifCNV
```

## Users Guide

Basic usage, ifCNV creates the output directory in wich it stores the html report that will open automatically.

```sh
$ ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/
```

More specific commands:

- if you want ifCNV to stop talking
```sh
$ ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/ -v ''
```
- if you don't want ifCNV to automatically open the report
```sh
$ ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/ -a ''
```
- if you don't want to re-create the reads matrix of your run at each utilisation (ie. to rerun ifCNV with different parameters)
First, run ifCNV and write the reads matrix:
```sh
$ ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/ -rm /path/to/readsMatrix/file
```
Then, run ifCNV and tell it to take the reads matrix (it will skip its creation)
```sh
$ ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/ -s /path/to/readsMatrix/file
```
- if you want to save the output in a .tsv file
```sh
$ ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/ -sv True
```
- changing mode from 'fast' (default) to 'extensive', will take a bit longer but doing so ifCNV will check every sample of the run for altered targets
```sh
$ ifCNV -i /path/to/bam/directory/ -b /path/to/bed/file -o /path/to/output/directory/ -m 'extensive'
```






