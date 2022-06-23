# Files for user testing

The files in the present directory are from a NGS amplicon dataset. 

## Changing the resolution

### Chromosome arm resolution

For the chromosome arm resolution, the reads matrix was obtained using the chromArm.bed file and typing this command:

```sh
ifCNV -i /path/to/bam/directory/ -b /path/to/chromArm.bed -o /path/to/output/directory/ -rm readsMatrix.tsv
```

Then, multiple parameters can be tested, the following gives the result in the output_chromArm folder:

```sh
ifCNV -i /path/to/bam/directory/ -b /path/to/chromArm.bed -o /path/to/output_chromArm/ -s readsMatrix.tsv -ct 0.2 -sT 50
```

### Gene resolution

For the gene resolution, the reads matrix was obtained using the gene.bed file and typing this command:

```sh
ifCNV -i /path/to/bam/directory/ -b /path/to/gene.bed -o /path/to/output/directory/ -rm readsMatrix.tsv
```

Then, multiple parameters can be tested, the following gives the result in the output_gene folder:

```sh
ifCNV -i /path/to/bam/directory/ -b /path/to/chromArm.bed -o /path/to/output_chromArm/ -s readsMatrix.tsv -ct 0.2 -sT 25
```


