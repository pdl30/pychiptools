#pychiptools 

### Installation

Clone this repository and then:

```bash
$ cd pychiptools/
$ python setup.py install --user
```

This will install the scripts in the pychiptools/scripts directory. For more information on the individual scripts, use the --help command after each script. 

##Core Pipeline

- pychip_align.py -> Wrapper for FASTQC, cutadapt and bowtie1/bowtie2 sequence aligner. 
- pychip_ucsc.py -> Takes a sam file and converts it to a UCSC formatted bigWig
- pychip_peak_call.py -> Peak calling using MACS2 and SICER

##Additional Tools
- pychip_peak_anno.py -> Peak annotation using chipseqanno, homer and also a custom annotation
- pychip_motifs.py -> Denovo motif discovery implemented using MEME and HOMER. Also motif search using FIMO
- pychip_diff_bind.py -> Differential binding implemented using a variety of programs. 
- pychip_download.py -> Tool for downloading samples from GEO
