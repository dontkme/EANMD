
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
![EANMD main feature](https://github.com/dontkme/PersonalScripts/raw/master/Fig.workflow.202402.2.feature-02-02.png )


## Workflow

![EANMD workflow](https://github.com/dontkme/PersonalScripts/raw/master/Fig.workflow.202402.2.flow-03.png)







## Installation

You can proceed to download the EANMD files from [here](https://github.com/dontkme/EANMD/archive/main.zip).
<details>
<summary>Simply unzip the downloaded zip file:</summary>


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
</details>
    
## Running Tests

To run tests, run the following command

<details>
<summary>0. Download reference genome and GTF annotation, if you don't had them before.</summary>

```bash
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz

  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
```
    
  Unzip them


  ```bash
  gunzip GRCm38.p6.genome.fa.gz
  gunzip gencode.vM25.primary_assembly.annotation.gtf.gz
  ```
  </details>
  <details>
  <summary>1. Run EANMD, test 28 mouse AS events.</summary>

  ```bash
  perl EANMD_GENOCDE -g gencode.vM25.primary_assembly.annotation.gtf -in TestMouseMM10_SE28.input.txt GRCm38.p6.genome.fa
  ```
  If it run, test pass.
</details>

## Usage/Examples

```
perl EANMD_GENCODE [options] -g <string|GTF annotation file> -in <string|input rMATs type list> <input FASTA file(s)>

    Options:
    [-o string|Output prefix. Default: getseqsOut]
    [-p int |Threads number. Default: 1]
    [-d int |Distance to the last exon-exon junction. Default: 50]
    [-s string|Specify attribute type in GTF annotation for sorting. Default: gene_id]
    [-f string|Specify feature type in GTF annotation. Default: '']
```
### Example of main EANMD steps:
1. Run EANMD, input GTF, AS events and Genome, ouput an outCombined.txt file.
   ```
   perl EANMD_GENOCDE -g gencode.vM25.primary_assembly.annotation.gtf -p 4 -o TestMouseMM10_SE28.EANMD -in TestMouseMM10_SE28.input.txt GRCm38.p6.genome.fa
   ```
2. Filter non-ATG transcripts and MXE results.
   ```
   perl EANMDFilterOutSE.pl -o TestMouseMM10_SE28.EANMD.f TestMouseMM10_SE28.EANMD.outCombined.txt
   ```
   **Note**: **SE** events use **EANMDFilterOutSE.pl**; **A5SS** events use **EANMDFilterOutA5SS.pl**; **A3SS** and **IR** events use **EANMDFilterOutA3SS.pl**

3. Summary the unique NMD flag for each AS events. (Need R environment)
    ```
    Rscript EANMDflagcount_reverse.R -i TestMouseMM10_SE28.EANMD.f.filter.ATG.txt -o TestMouseMM10_SE28.EANMD.f.filter.ATG.txt.flag
    ```
### Optional Steps
1. Get an AS-evetns list from filtering rMATS results.
   ```
   perl GetSEinput.pl [options] <input rmats result>
   ```

   <details>
   <summary>Options</summary>

         [-o output prefix. default: rMATS_filtered.out]

         [-d int|min depth of average read counts. default: 20]

         [-m int|min count of inclusion events' UP|Downstream junction. default: 2]

         [-i float|delta PSI cutoff [0-1.0]. default: 0.15]

         [-f float|FDR cutoff [0-1.0]. default: 0.05]

         [-c1 int|The first sample numbers.  default: 2]

         [-c2 int|The second sample numbers. default: 2]

         [-mf float|The US and DS fold change cutoff. default: 0.05]
    </details>

    Example: 

    ```
    
    ```
     **Note**: **SE** events use **GetSEinput.pl**; **A5SS** events use **GetA5SSinput.pl**; **A3SS** events use **GetA3SSinput.pl**


2. Convert unique flag results to GTF
   ```
   perl TransFlag2GTF.pl [options] <input AS flag file>
   ```

## File format

### 1. AS events input: 3 Columns Tab delitimated file.

|Column|Description|Example|
|:------|:------------|:-------|
|Gene_name|Gene Symbol|Flna|
|AS|Alternative splicing event SEUSDS Coordinates: <br>SE@US@DS<br><br>SE\|US\|DS = chr:start0:end:strand<br>|chrX:74240815:74240872:-@chrX:74240412:74240550:-@chrX:74241102:74241303:-|
|Optional|Custom column, will remain in output, suggest use gene_id|ENSMUSG00000031328|
### 2. OutCombined Output file: 38 Columns
<details>
<summary>Descriptions of the 38 columns</summary>

|Column|Description|
|:---|:---|
|QueryCol1|Input Column 1|
|SEUSDSCoordinates|Input Column 2|
|QueryCol3|Input Column 3|
Transcript_id|Reference transcript id|
Strand|Transcript strand|
Exons|Total exon numbers of the reference transcript|
Start_exon|Reference start condon exon number|
Stop_exon|Reference stop condon exon number|
SE_exon_Number|Skipped exon number for reference transcript|
SE(US)_Pos|Skipped exon position for reference transcript|
SE_length|Skipped exon length|
Ori_CDS_length|Original CDS length|
Ori_Star_codon_to_exon_end_seq_len|Length of start-codon to exon end|
rm/add_SE_start_to_end_seq_len|Length of start-codon to exon end after remove or add the SE|
SEseq|Skipped exon sequence|
Ori_CDSexons_seq|Original Start codon to exons end sequence|
rm/add_SE_CDSexons_seq|Sequence of start codon to exons end after remove or add the SE|
Ori_last_junction_pos|Original last exon-exon junction position|
Ori_last_dj|Distance of original stop codon to the last exon-exon junction (EJ)|
Ori_NMD|Original or refernce transcript NMD type|
Start_codon|Star_codon sequence|
Ori_AA|Original amino acid sequence|
rm/add_SE_AA|Amino acid sequence after remove or add the SE|
AA_len+1|Original amino acid sequence length + 1 (stop codon)|
Ori_AA_1st_stop_pos|Original amino acid 1st stop codon positon|
Ori_AA_stop_pos|Original amino acid stop codon positons|
SEed_AA_1st_stop_pos|Amino acid 1st stop codon positon after remove or add the SE|
SEed_AA_stop_pos|Amino acid stop codon positons after remove or add the SE|
Frame_shift_flag|Frame shift flag|
New_1st_stop_pos_dj|Distance of the new 1st stop codon to the last EJ|
NMD_flag|NMD flag after remove or add the SE|
NMD_in/ex_flag|This AS event NMD type|
source|This record form which median file|
SEupstreamCDS|How many nt upstream of the SE|
SEupstreamAApos|How many AA upstream of the SE|
UStransexonnumber|US exon number in reference transcript|
DStransexonnumber|DS exon number in reference transcript|
innerExonsofUSandDS|inner exon(s) between US and DS|
</details>
   
   
## Authors

- [@Kaining Hu](https://www.github.com/dontkme) - *Initial work* -




## License
![licence](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details