
# EANMD
<img src="https://github.com/dontkme/PersonalScripts/raw/master/Fig.logo.EANMD-02.png"  align="right" height="73" width="221"/>

[![releaseVersion](https://img.shields.io/badge/release%20version-1.44-green.svg?style=flat)](https://github.com/dontkme/EANMD/releases) [![Last-changedate](https://img.shields.io/badge/last%20change-2024--11--06-green.svg)](https://github.com/dontkme/EANMD/commit) ![perlVersion](https://img.shields.io/badge/perl-%3E%3D5.10-blue.svg?style=flat)

**Easy Annotation of Alternative Splicing Events Coupled with Nonsense-Mediated mRNA Decay** (EANMD) is written in Perl to predict alternative splicing's potential to trigger NMD based on 50-nt rule: premature stop-codon before last exon-exon junctions (DJ) more than 50 nt. XGBoost is used to predict the NMD efficiency score based on transcripts factors.

## Features

- Support **long-read** and **short-read** RNA-seq sourced AS events.
- Support **wide range** **species**.
- Support predicting NMD types for original and new isoforms, as well as **NMD_in** and **NMD_ex** clustering.
- Support predicting NMD types for novel Skipped exon (**SE**), Intron Retention (**IR**), Alternative 5' splicing site (**A5SS**), and Alternative 3' splicing site (**A3SS**) events.
- Support **customizing the 50-nt rule**.
- Support **multi-threading**.
- Support filtering out non-ATG-started transcripts, MXE events
- Support sorting and scoring AS events using our trained **Xgboost model**.


<!-- ![EANMD main feature](https://github.com/dontkme/PersonalScripts/raw/master/Fig.workflow.202402.2.feature-02-02.png ) -->


## Workflow Overview
![EANMD workflow](https://github.com/dontkme/PersonalScripts/raw/master/Fig.workflow.202402.2.flow-03.png)

### 1). AS Events Detection from RNA-seq or GTF Annotation
<details> 
<summary>Click to see the AS detection steps from long-read and short-read RNA-seq</summary>

![AS Detection](https://github.com/dontkme/PersonalScripts/raw/master/1.AS_events.png)

</details>

### 2. AS-NMD Prediction

![EANMD main feature](https://github.com/dontkme/PersonalScripts/raw/master/Fig.workflow.202402.2.feature-02-02.png )

### 3). Filtering, Sorting and Scoring AS-NMD Events

<details>

<summary>Click to see the Filtering, Sorting and Scoring steps</summary>

![NMD Filtering, Sorting and Scoring](https://github.com/dontkme/PersonalScripts/raw/master/3.FilterOutCombined.png)

</details>

## Installation

You can download the EANMD from [here](https://github.com/dontkme/EANMD/archive/main.zip).
<details>
<summary>Simply unzip the downloaded zip file:</summary>


```bash
unzip EANMD-main.zip
```

Navigate to the extracted folder and run EANMD:

```
cd EANMD-main
perl EANMD_GENCODE -h
```

If the screen displays usage and version information. It works.

**If** it displays you need `Array::Utils` module, please install it by:
```
cpan Array::Utils
```

**If** you need Perl `Parallel::ForkManager` module. You could install it by command: 

```bash
cpan Parallel::ForkManager
``` 
</details>
    
## Running Tests

To run tests, run the following command

<details>
<summary>0. Download reference genome and GTF annotation, if you don't have them before.</summary>

```bash
  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz

  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
```
    
  Unzip them:


  ```bash
  gunzip GRCm38.p6.genome.fa.gz
  gunzip gencode.vM25.primary_assembly.annotation.gtf.gz
  ```
  </details>
  <details>
  <summary>1. Run EANMD, test the 28 mouse AS events.</summary>

  ```bash
  perl EANMD_GENCODE -g gencode.vM25.primary_assembly.annotation.gtf -in TestDataset/TestMouseMM10_SE28.input.txt GRCm38.p6.genome.fa
  ```
  If it runs, the test passes.
</details>

## Usage/Examples

```
perl EANMD_GENCODE [options] -g <string|GTF annotation file> -in <string|input AS-event list> <input genome FASTA file(s)>

    Options:
    [-o string|Output prefix. Default: getseqsOut]
    [-p int |Threads number. Default: 1]
    [-d int |Distance to the last exon-exon junction. Default: 50]
    [-s string|Specify attribute type in GTF annotation for sorting. Default: gene_id]
    [-f string|Specify feature type in GTF annotation. Default: '']
```
### Example of main EANMD steps:
1. Run EANMD, input GTF, AS events and Genome, output an outCombined.txt file.
   ```
   perl EANMD_GENCODE -g gencode.vM25.primary_assembly.annotation.gtf -p 4 -o TestMouseMM10_SE28.EANMD -in TestMouseMM10_SE28.input.txt GRCm38.p6.genome.fa
   ```
   Main output file: `TestMouseMM10_SE28.EANMD.outCombined.txt`

   **Note:** The default version of EANMD only used the US End and DS Start boundaries to fuzzy search the reference transcripts. If you want to search the exact coordinates of US and DS boundaries (contain both Start and End), you could use the `EANMD_GENCODE_strict` version.
2. Filter non-ATG-started transcripts and MXE results.
   ```
   perl EANMDFilterOutSE.pl -o TestMouseMM10_SE28.EANMD.f TestMouseMM10_SE28.EANMD.outCombined.txt
   ```
   Main output file: `TestMouseMM10_SE28.EANMD.f.filter.ATG.outCombined.txt`
   
   **Note:** **SE** events use **EANMDFilterOutSE.pl**; **A5SS** events use **EANMDFilterOutA5SS.pl**; **A3SS** and **IR** events use **EANMDFilterOutA3SS.pl**

3. Summarize the unique **NMD Flag** for each AS events. (Requires R environment)
    ```
    Rscript EANMDflagcount.R -i TestMouseMM10_SE28.EANMD.f.filter.ATG.outCombined.txt -o TestMouseMM10_SE28.EANMD.f
    ```
    Main output file: `TestMouseMM10_SE28.EANMD.f.FinalUniqNMDflag.txt`

    **Note:** Could use `EANMDflagcount_reverse.R` to get reverse ordered flags.  The 'no_info', 'no_NMD', 'No_stop_codon' > 'NMD_in' or 'NMD_ex', the final flag will be the most frequent flag.


4. **Optional**: If you want to get the **NMD** efficiency **score**, you could use the `EANMDflagcount_withUTRMFE.R` script. This script will calculate the NMD Score based on the 3'-UTR length, the maximum distance of stop codon to the last exon-exon junction (Max DJ), 3'-UTR MFE/nt, Minimal stop codon position fraction (Min stop Pos F) and other factors (See Paper Methods). [*RNAfold*](https://github.com/ViennaRNA/ViennaRNA), [*xgboost*](https://xgboost.readthedocs.io/en/stable/install.html#r) R package are required for this step. Higher the NMD Score, higher the NMD efficiency. We could filter the NMD events by the NMD Score, for example, NMD Score > 0.131.

    ```
    Rscript EANMDflagcount_withUTRMFE.R -i TestMouseMM10_SE28.EANMD.f.filter.ATG.outCombined.txt -o TestMouseMM10_SE28.EANMD.f
    ```

### Other Optional Steps
1. Get an AS-events list from filtering rMATS individual-counts results.
   ```
   perl GetSEinput.pl [options] <input rmats individual-counts result>
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
    perl GetSEinput.pl -o Test.SE -d 2 -m 2 -i 0.1 -f 0.05 -c1 2 -c2 2 -mf 0.05 SE.MATS.JCEC.txt
    ```
    **-d 2**: minimum depth of average read counts per sample is 2.
    **-m 2**: minimum junction count for SE inclusion is 2.
    **-i 0.1**: delta PSI cutoff is 0.1.
    **-f 0.05**: max FDR is 0.05.
    **-c1 2** and **-c2 2**: Group1 has 2 samples, Group2 has 2 samples.
    **-mf 0.05** The US and DS fold change cutoff is 0.05, meaning the fold of upstream exon or downstream exon to skipped exon junctions reads must be greater than 1:20.

    Output file: `Test.SE.EANMD.input.txt`. Can be used in `EANMD` main process.

     **Note**: **SE** events use **GetSEinput.pl**; **A5SS** events use **GetA5SSinput.pl**; **A3SS** events use **GetA3SSinput.pl**


2. Convert unique flag results to GTF
   ```
   perl TransFlag2GTF.pl [options] <input AS flag file>
   ```

## Main File Format

### 1. AS events input file: 3 Columns Tab delimited file.

|Column|Description|Example|
|:------|:------------|:-------|
|Gene_name|Gene Symbol|Flna|
|AS|Alternative splicing event SEUSDS Coordinates: <br>SE@US@DS<br><br>SE\|US\|DS = chr:start0\:end:strand<br>|chrX:74240815:74240872:-@chrX:74240412:74240550:-@chrX:74241102:74241303:-|
|Optional|Custom column, will remain in output, suggest use gene_id|ENSMUSG00000031328|
### 2. OutCombined Output file: 38 Columns Tab delimited file.
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
Start_exon|Reference start codon exon number|
Stop_exon|Reference stop codon exon number|
SE_exon_Number|Skipped exon number for reference transcript|
SE(US)_Pos|Skipped exon position for reference transcript|
SE_length|Skipped exon length|
Ori_CDS_length|Original CDS length|
Ori_Start_codon_to_exon_end_seq_len|Length of start-codon to exon end|
rm/add_SE_start_to_end_seq_len|Length of start-codon to exon end after remove or add the SE|
SEseq|Skipped exon sequence|
Ori_CDSexons_seq|Original Start codon to exons end sequence|
rm/add_SE_CDSexons_seq|Sequence of start codon to exons end after remove or add the SE|
Ori_last_junction_pos|Original last exon-exon junction position|
Ori_last_dj|Distance of original stop codon to the last exon-exon junction (EJ)|
Ori_NMD|Original or reference transcript NMD type|
Start_codon|Start_codon sequence|
Ori_AA|Original amino acid sequence|
rm/add_SE_AA|Amino acid sequence after remove or add the SE|
AA_len+1|Original amino acid sequence length + 1 (stop codon)|
Ori_AA_1st_stop_pos|Original amino acid 1st stop codon position|
Ori_AA_stop_pos|Original amino acid stop codon positions|
SEed_AA_1st_stop_pos|Amino acid 1st stop codon position after remove or add the SE|
SEed_AA_stop_pos|Amino acid stop codon positions after remove or add the SE|
Frame_shift_flag|Frame shift flag|
New_1st_stop_pos_dj|Distance of the new 1st stop codon to the last EJ|
NMD_flag|NMD flag after remove or add the SE|
NMD_in/ex_flag|This AS event NMD type|
source|This record from which intermediate file|
SEupstreamCDS|How many nt upstream of the SE|
SEupstreamAApos|How many AA upstream of the SE|
UStransexonnumber|US exon number in reference transcript|
DStransexonnumber|DS exon number in reference transcript|
innerExonsofUSandDS|inner exon(s) between US and DS|
</details>
   
   
## Author

- [@Kaining Hu](https://www.github.com/dontkme) - *Initial work* -




## License
![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details