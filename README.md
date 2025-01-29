# Mtb-barcoding-amplicon

## Project Description
This project describes the computational workflow used to identify barcoded Mycobacterium tuberculosis (Mtb) in tissue extracts from animals infected with a library of barcoded strains. More specifically, the workflow is applied to FASTQ files generated as follows: Animals are infected with a genetically barcoded strain of Mycobacterium tuberculosis (Mtb), strain Erdman described in Martin CJ, Cadena AM, Leung VW, Lin PL, Maiello P, Hicks N, Chase MR, Flynn JL, Fortune SM. Digitally Barcoding Mycobacterium tuberculosis Reveals In Vivo Infection Dynamics in the Macaque Model of Tuberculosis. mBio. 2017 May 9;8(3):e00312-17. doi: 10.1128/mBio.00312-17. PMID: 28487426; PMCID: PMC5424202. Bacteria is cultured out of harvested tissues, followed by bacterial DNA extraction, and then the amplification of barcoded region using PCR. These amplicons are then pooled across tissues and subject to Illumina sequencing.

**Please don't hesitate to contact us if there are any issues.**

### Directories:
- ./scripts:
  - BarcodeReader.plx
  - FindThreshold.plx
  - download.sh
  - parallel.sh
  - barcoding.yml
- ./data:
  - PRJNA1215398_sra_accession_list.txt
- barcoding.yml
    
### Prerequisites:
  - Requires Perl, seqtk https://github.com/lh3/seqtk and the perl module Math::Derivative https://metacpan.org/pod/Math::Derivative to be installed on your path.
  - To download the fastq files from SRA you will also need to install sra-tools: https://github.com/ncbi/sra-tools. To speed things up a little gnu parallel: https://www.gnu.org/software/parallel/

### Installation:
  - Download the repository and copy the scripts to a directory on your path and make them executable `chmod u+x` or just copy them to your directory with your fastq files.
  - The tool versions can be manually installed in a unix like environment running bash or can be install using a Anaconda environments.
  - `conda env create -f barcoding.yml`
  - If installed with conda, in the script FindThreshold.plx change the path in `use lib '/my_path_to_environment/barcoding/lib/perl5/site_perl/';` to the path of your environment installation.
  - On a system with a native perl module installation, delete or comment out '/my_path_to_environment/barcoding/lib/perl5/site_perl/' in the script FindThreshold.plx;

### Download the fastq files:
  - Copy download.sh && parallel.sh into your data analysis directory.
  - Modify the script parallel.sh to the number of cpus you want to use.
  - `chmod u+x download`
  - Copy 'PRJNA1215398_sra_accession_list.txt' into your directory.
  - `bash parallel.sh PRJNA1215398_sra_accession_list.txt`

### Running the code:
- `cp /path_to/BarcodeReader.plx .`
- `cp /path_to/FindThreshold.plx .`
- The script expects FASTQ files in the working directory as named in SRA: SRR32115076_1.fastq.gz.
- `perl BarcodeReader.plx`
- `perl FindThreshold.plx`
- outputs the file: directory_name_threshold_data.csv
  
### Troubleshooting:
- Verify read quality and counts.
- Ensure that the prefix of the FASTQ file conforms to the expected pattern. if not you will have to modify the subroutine get_file() in both scripts to match your file name convention.
- Accurate threshold is not always achieved and can require additional curration.

### Citation and associated data:
The script was used to analyze bacterial barcoding data associated with the manuscript "CD4 T cells and CD8α+ lymphocytes are necessary for intravenous BCG-induced protection against tuberculosis in macaques" by Simonson et al. (2024, doi: https://doi.org/10.1101/2024.05.14.594183)
- The fastqs are on SRA at https://www.ncbi.nlm.nih.gov/sra/PRJNA1215398
- The metadata is on Fairdomhub at https://fairdomhub.org/studies/1283

[![DOI](https://zenodo.org/)] **(https://doi.org/10.5281/zenodo.14766941)**



