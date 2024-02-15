
# EANMD
<img src="https://github.com/dontkme/PersonalScripts/raw/master/Fig.logo.EANMD-02.png"  align="right" height="73" width="221"/>

[![releaseVersion](https://img.shields.io/badge/release%20version-1.42-green.svg?style=flat)](https://github.com/dontkme/EANMD/releases) [![Last-changedate](https://img.shields.io/badge/last%20change-2023--7--11-green.svg)](https://github.com/dontkme/EAHNMD/commit) ![perlVersion](https://img.shields.io/badge/perl-%3E%3D5.10-blue.svg?sytle=flat)

**Exon annotation for nonsense-mediated mRNA decay** (EANMD) is written in Perl to predict alternative splicing effects on potential to trigger NMD by 50-nt rules: premature stop-codon before last exon-exon junctions more than 50 nt.

## Features

- Support Alternative splicing events: Skipped exon (**SE**), Intron Retention (**IR**), Alternative 5' splicing site (**A5SS**) and Alternative 3' splicing site (**A3SS**)
- Support predict original isoform NMD type, new isoform NMD type, **NMD_in** and **NMD_ex** clustering.
- Support novel cassette exon, novel IR, A5SS and A3SS events NMD type prediction
- Support customize 50-nt rule
- Support multi-threaded


## Workflow

![EANMD workflow](https://github.com/dontkme/PersonalScripts/raw/master/Fig.workflow.202402-01.png)







## Installation

you can proceed to download the EANMD files from [here](https://github.com/dontkme/EANMD/archive/main.zip).

Simply unzip the downloaded zip file:


```bash
unzip EANMD-master.zip
```

Navigate to the extracted folder and run EANMD:

```
cd EANMD-main
perl EANMD_GENCODE -h
```

If the screen displays help and version information. It works.

**If need Perl Parallel::ForkManager.** You could install it by command: 

```bash
cpan Parallel::ForkManager
``` 
    
## Running Tests

To run tests, run the following command

0. Download reference genome and GTF annotation, if you don't had them before.
```bash
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz

  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
```
    
  Unzip them


  ```bash
  gunzip GRCm38.p6.genome.fa.gz
  gunzip gencode.vM25.primary_assembly.annotation.gtf.gz
  ```
  1. Run EANMD, test 28 mouse AS events.
  ```bash
  perl EANMD_GENOCDE -g gencode.vM25.primary_assembly.annotation.gtf -in TestMouseMM10_SE28.input.txt GRCm38.p6.genome.fa
  ```
  If it run, test pass.


## Usage/Examples

```bash
perl EANMD_GENCODE [options] -g <string|GTF annotation file> -in <string|input rMATs type list> <input FASTA file(s)>

    Options:
    [-o string|Output prefix. Default: getseqsOut]
    [-p int |Threads number. Default: 1]
    [-d int |Distance to the last exon-exon junction. Default: 50]
    [-s string|Specify attribute type in GTF annotation for sorting. Default: gene_id]
    [-f string|Specify feature type in GTF annotation. Default: '']
```


## Authors

- [@Kaining Hu](https://www.github.com/dontkme) - *Initial work* -


## Appendix

Any additional information goes here


## Acknowledgements

 - [Awesome Readme Templates](https://awesomeopensource.com/project/elangosundar/awesome-README-templates)
 - [Awesome README](https://github.com/matiassingers/awesome-readme)
 - [How to write a Good readme](https://bulldogjob.com/news/449-how-to-write-a-good-readme-for-your-github-project)


## 🔗 Links
[![portfolio](https://img.shields.io/badge/my_portfolio-000?style=for-the-badge&logo=ko-fi&logoColor=white)](https://katherineoelsner.com/)
[![linkedin](https://img.shields.io/badge/linkedin-0A66C2?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/)
[![twitter](https://img.shields.io/badge/twitter-1DA1F2?style=for-the-badge&logo=twitter&logoColor=white)](https://twitter.com/)

## License
![licence](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details