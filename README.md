# Mtb-barcoding-amplicon

## Project Description
This project describes the computational workflow used to identify barcoded Mycobacterium tuberculosis (Mtb) in tissue extracts from animals infected with a library of barcoded strains. More specifically, the workflow is applied to FASTQ files generated as follows: Animals are infected with a genetically barcoded strain of Mycobacterium tuberculosis (Mtb), strain Erdman described in Martin CJ, Cadena AM, Leung VW, Lin PL, Maiello P, Hicks N, Chase MR, Flynn JL, Fortune SM. Digitally Barcoding Mycobacterium tuberculosis Reveals In Vivo Infection Dynamics in the Macaque Model of Tuberculosis. mBio. 2017 May 9;8(3):e00312-17. doi: 10.1128/mBio.00312-17. PMID: 28487426; PMCID: PMC5424202. Bacteria is cultured out of harvested tissues, followed by bacterial DNA extraction, and then the amplification of barcoded region using PCR. These amplicons are then pooled across tissues and subject to Illumina sequencing.

**Please don't hesitate to contact us if there are any issues.**

### Directories:
- ./scripts:
  - BarcodeReader.plx
  - FindThreshold.plx
  - barcoding.yml
- ./data:
  - SraRunInfo.csv
  - 
### Prerequisites:
  - Requires Perl, seqtk https://github.com/lh3/seqtk and the perl module Math::Derivative https://metacpan.org/pod/Math::Derivative to be installed on your path.
  - To download the fastq files from SRA you will also need to install sra-tools: https://github.com/ncbi/sra-tools.

### Installation:
  - Download the repository and copy the scripts to a directory on your path and make them executable `chmod u+x` or just copy them to your directory with your fastq files.
  - 
  - The tool versions can be manually installed in a unix like environment running bash or can be install using a Anaconda environments.
  - `conda env create -f barcoding.yml`
  - In the script FindThreshold.plx change the path in `use lib '/my_path_to_environment/barcoding/lib/perl5/site_perl/';` to the path of your environment installation.

### Running the code:
- `mkdir myFolder && cd myFolder`
- `conda activate barcoding`
- `cp /path_to/BarcodeReader.plx .`
- `cp /path_to/FindThreshold.plx .`
- `mv fastq_files . `
- The script expects fastqs as named in SRA: SRR32115076_1.fastq.gz.
- `perl BarcodeReader.plx`
- `perl FindThreshold.plx`
- output directory_name_threshold_data.csv
  
### Troubleshooting:
- Verify read quality and counts.
- Ensure that the prefix of the fastq file conforms to the expected pattern. if not you will have to modify the subroutine get_file() in both scripts to match your file name convention.
- Accurate threshold is not always achieved and can require additional curration.

### Data:
Data associated with the manuscript "CD4 T cells and CD8α+ lymphocytes are necessary for intravenous BCG-induced protection against tuberculosis in macaques" by Simonson et al. (2024, doi: https://doi.org/10.1101/2024.05.14.594183)
- SRA: https://www.ncbi.nlm.nih.gov/sra/PRJNA1215398
- Fairdomhub: https://fairdomhub.org/studies/1283

### DOI of this repository: (Link to the Zenodo DOI, generate a unqiue DOI for this repo)
[![DOI](https://zenodo.org/)] **Zenodo will give you a link that you can paste here**



