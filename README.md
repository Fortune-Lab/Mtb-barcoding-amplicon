# Project Title
Data associated with the manuscript "CD4 T cells and CD8α+ lymphocytes are necessary for intravenous BCG-induced protection against tuberculosis in macaques" by Simonson et al. (2024, doi: https://doi.org/10.1101/2024.05.14.594183)

## Project Description
This project describes the computational workflow used to identify barcoded Mycobacterium tuberculosis (Mtb) in tissue extracts from animals infected with a library of barcoded strains. More specifically, the workflow is applied to FASTQ files generated as follows: Animals are infected with a genetically barcoded strain of Mycobacterium tuberculosis (Mtb), strain Erdman described in Martin CJ, Cadena AM, Leung VW, Lin PL, Maiello P, Hicks N, Chase MR, Flynn JL, Fortune SM. Digitally Barcoding Mycobacterium tuberculosis Reveals In Vivo Infection Dynamics in the Macaque Model of Tuberculosis. mBio. 2017 May 9;8(3):e00312-17. doi: 10.1128/mBio.00312-17. PMID: 28487426; PMCID: PMC5424202. Bacteria is cultured out of harvested tissues, followed by bacterial DNA extraction, and then the amplification of barcoded region using PCR. These amplicons are then pooled across tissues and subject to Illumina sequencing.

**Please don't hesitate to contact us if there are any issues.**

## File Guide
This project includes two perl scripts to process the fastq files in the analysis of this work. 

### Directory:
- ./scripts:
  - BarcodeReader.plx
  - FindThreshold.plx
    
### Tool versions:
  - Requires Perl, seqtk https://github.com/lh3/seqtk and the perl module Math::Derivative https://metacpan.org/pod/Math::Derivative to be installed on your path.

### Installation:
  - The tool versions can be manually installed in a unix like environment running bash. You can also install using a Anaconda environments.
  - `conda env create -f barcoding.yml`
  - In the script FindThreshold.plx change the path in `use lib '/n/boslfs02/LABS/sfortune_lab/Lab/envs/barcoding/lib/perl5/site_perl/';` to the path of your environment installation.

### Running the code:


### DOI of this repository: (Link to the Zenodo DOI, generate a unqiue DOI for this repo)
[![DOI](https://zenodo.org/)] **Zenodo will give you a link that you can paste here**



