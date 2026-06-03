# Kivvi: HiFi genotyper for large-unit variable number tandem repeat

Kivvi is a HiFi-based tool that calls the copy number and sequence variants of large-unit variable number tandem repeats (VNTRs). With the large size of the full repeat and the highly polymorphic nature across individuals, these regions are challenging to resolve using currently available methods. Kivvi identifies unique copies of each repeat and assembles them into alleles. Kivvi has been applied to two medically important VNTRs and can be adapted to more large-unit VNTRs.
- The LPA Kringle IV-type 2 (KIV2) repeat (repeat unit 5.5kb). A short KIV2 allele is associated with a higher risk of cardiovascular diseases.
- The D4Z4 repeat (repeat unit 3.3kb). D4Z4 is involved in [Facioscapulohumeral Muscular Dystrophy (FSHD)](https://www.ncbi.nlm.nih.gov/books/NBK1443/), which is caused by chromatin relaxation (hypomethylation) and/or contraction of D4Z4.

For more details about Kivvi, please check out our latest [preprint](https://www.biorxiv.org/content/10.64898/2026.04.10.717730) on D4Z4.

# Table of Contents

- [Contact](#contact)
- [Installation](#installation)
- [Input](#input)
- [Running the program](#running-the-program)
- [Output](#output)
- [Demo and tutorials](#demo-and-tutorials)

## Contact

If you have suggestions or need assistance, please don't hesitate to reach out by email or open a GitHub issue.

Xiao Chen: xchen@pacificbiosciences.com

## Installation

```bash
# Specify the version
VERSION="v1.1.0"
# Download the release file
wget https://github.com/PacificBiosciences/kivvi/releases/download/${VERSION}/kivvi-${VERSION}-x86_64-unknown-linux-gnu.tar.gz
# Decompress the file
tar -xzvf kivvi-${VERSION}-x86_64-unknown-linux-gnu.tar.gz
cd kivvi-${VERSION}-x86_64-unknown-linux-gnu
# Check the md5 sum (optional)
md5sum -c kivvi.md5
# Execute help instructions
./kivvi -h
```

## Input

The input to Kivvi is a WGS bam (aligned to GRCh38). The WGS must be standard depth (20-30X or higher). Kivvi works better with higher coverage and longer reads. Targeted data is not supported due to the shorter read length.

Kivvi can take a bamlet of the WGS bam as input. The region needed is (GRCh38):
- KIV2: `chr6:160605000-160655000`.
- D4Z4: `chr4:190022510-190093263 chr4:190173122-190192666 chr10:133622567-133685491 chr10:133740609-133775186`. 

## Running the program

Kivvi requires a genome-aligned BAM, an output directory and a prefix to output files. The command ends in a preset (`kiv2` or `d4z4`) for specifying which target region to run.

```bash
kivvi -b $WGS_BAM -o $OUTPUT_DIRECTORY -p $OUTPUT_PREFIX kiv2
```
Or
```bash
kivvi -b $WGS_BAM -o $OUTPUT_DIRECTORY -p $OUTPUT_PREFIX d4z4
```

Kivvi is single-threaded and generally takes less than 5 minutes per sample.

## Output

- `$prefix.kivvi.$target.json`: detailed report of alleles, variants and other information for each sample. See more detailed documentation [here](docs/json.md).
- `$prefix.kivvi.$target.svg`: reads plotted onto assembled alleles for visualization. (Only produced when at least one allele is assembled).
- `$prefix.kivvi.$target.bam`: all reads realigned to one repeat unit. Can be loaded into IGV to view different repeat copies. (Group by HP tag).
- `$prefix.kivvi.$target.vcf`: small variant calls. See more detailed documentation [here](docs/vcf.md). 

For KIV2, the reference used for read alignment (BAM) and variant calling (VCF) is one copy of the repeat on GRCh38 (`chr6:160613619-160619170`). Please navigate to this region for visualizing the Kivvi produced BAM and VCF files. 

For D4Z4, as the region is noisy on GRCh38, a different [reference](data/d4z4/d4z4_ref.fa) is used for read alignment (BAM) and variant calling (VCF) and can be loaded into IGV as the reference for visualizing the Kivvi produced BAM and VCF files. 

## Demo and tutorials

Tutorials are available for understanding Kivvi output files:
- [JSON](docs/json.md)
- [Visualization (BAM/SVG)](docs/visualization.md)
- [VCF](docs/vcf.md)

In these tutorials, Sample HG03453 is used as an example. A bamlet of the KIV2 region for this sample is available [here](example/HG03453_kiv2_extract.bam). A bamlet of the D4Z4 region for this sample is available [here](example/HG03453_d4z4_extract.bam). Note that a warning of `Genome depth is unavailable` is expected when running Kivvi with these two bamlets, as these are not full WGS bams.

```bash
# Download human GRCh38 if you don't have one
wget https://downloads.pacbcloud.com/public/reference-genomes/human_GRCh38_no_alt_analysis_set.tar.2023-12-04.gz
tar -xpvf human_GRCh38_no_alt_analysis_set.tar.2023-12-04.gz
# Download the demo data for KIV2
wget https://raw.githubusercontent.com/PacificBiosciences/kivvi/main/example/HG03453_kiv2_extract.bam
wget https://raw.githubusercontent.com/PacificBiosciences/kivvi/main/example/HG03453_kiv2_extract.bam.bai
# Run Kivvi for KIV2
kivvi -b HG03453_kiv2_extract.bam -o ./kiv2_output/ -p HG03453 kiv2
# Download the demo data for D4Z4
wget https://raw.githubusercontent.com/PacificBiosciences/kivvi/main/example/HG03453_d4z4_extract.bam
wget https://raw.githubusercontent.com/PacificBiosciences/kivvi/main/example/HG03453_d4z4_extract.bam.bai
# Run Kivvi for D4Z4
kivvi -b HG03453_d4z4_extract.bam -o ./d4z4_output/ -p HG03453 d4z4
```

## DISCLAIMER

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
