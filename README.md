# uORFs
This repository contains most of the data and all the scripts used during my final grade project (Apr-Jun 2020)

## Contents
- FGP_RESULTS:
Contains the custom track that can be uploaded to UCSC genome browser, and the code needed to create it and also the files with the scored uORFs using the 5 methods proposed.

- SCRIPTS:
Contains the scripts for iedntifying all the uORFs in a given genome and also for adding all the additional features. 

### uORFs_identifier

The uORFs identifier tool can be runned using a command like this:

`perl uORFs_identifier.pl -i ../../DATA/input_data/all_human_EnsemblGeneIDs_v34.txt -of ../../DATA/raw_data/raw_uORFs/allENSG_15-5-2020_at_17-12.tsv -m 100 -sp Human 2> error`

#### Requirements

Ensembl Core API installed
see http://www.ensembl.org/info/docs/api/api_installation.html for the installation procedure.
GitTools installation recommended in order to manage version updates easily.


### Adding additional features information to each uORF.

GENCODE annotation files available at:
www.gencodegenes.org

All conservations scores source files can be downloaded from these sites:
- phastCons: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.bw
- phyloP: http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw
- phyloCSF: https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg38_100/20170118/ all PhyloCSF*.bw files





