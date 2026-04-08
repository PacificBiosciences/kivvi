# Kivvi JSON

Kivvi produces a single JSON file per sample, containing detailed information about alleles, variants, and other information such as methylation levels when available.

## Fields shared by KIV2 and D4Z4

- `allele_cn`: copy number of the repeat on assembled alleles. It’s empty when zero allele is assembled. When there are multiple assembled alleles, the numbers are `/`-separated for KIV2 (diploid) and comma-separated for D4Z4 (expects four alleles). 
- `depth_summary`: reports average depth of the genome and average depth of the repeat (after realignement to one repeat unit).
- `complete_alleles`: completely assembled alleles, represented as an ordered list of repeat unit IDs.
- `partial_alleles`: partially assembled alleles
- `supporting_reads`: reads supporting completely assembled alleles
- `complete_allele_variants`: variants on completely assembled alleles. Variants are reported per allele per repeat unit (`{unit index on allele};{unit ID};{variant list}`). 
- `other_unit_variants`: variants on any repeat units that are not assembled into complete alleles (`{unit ID};{variant list}`).

Part of the HG03453 KIV2 JSON file is shown below. HG03453 has 14 KIV2 units on one allele and 23 units on the other.
```json
{
  "allele_cn": "14/23",
  "complete_alleles": [
    "LeftFlank-8-11-9-6-1-6-10-1-1-4-7-2-3-5-RightFlank",
    "LeftFlank-14-13-31-30-12-27-19-28-23-29-18-17-15-16-15-25-24-20-24-23-18-26-22-RightFlank"
  ],
}
```

## D4Z4 specific fields

FSHD is caused by D4Z4 chromatin relaxation (hypomethylation) and ectopic expression of the DUX4 gene, encoded in D4Z4, in skeletal muscle, mediated by contraction of D4Z4 to 1-10 copies (FSHD1, 95% of FSHD cases) or mutations in chromatin factor genes (FSHD2, 5% of FSHD cases). Therefore, both the size and the methylation information matter for D4Z4. The D4Z4 repeat exists on both chr4 and chr10, while a contracted D4Z4 allele is only pathogenic on a specific haplotype (qA, with intact PolyA signal). Kivvi collects D4Z4 reads from both chromosomes and assembles them, and assigns each allele to one of the three allele types:
- `qAIntactPolyA`
- `qB`
- `qADisruptedPolyA`

Only `qAIntactPolyA` alleles are able to express DUX4 and thus cause FSHD. 

D4Z4 alleles can be very long (1-100 copies, 3.3kb each). Longer alleles are difficult to assemble and they are represented as partial alleles if not fully assembled - in this case, a lower bound of allele size is reported. 

The most important field to check for D4Z4 is the `allele_info` field:
- `allele_info`: reports all alleles, including partially assembled alleles. For each allele, it reports the size (a lower bound if partially assembled), originating chromosome, distal haplotype (status of polyA signal) and methylation level.
  - allele_name: alleles represented as an ordered list of repeat unit IDs.
  - chromosome: chromsome 4 or chromosome 10, followed by assignment of upstream haplotype groups, which are discussed in our preprint (link available soon). The upstream haplotype group provides extra information since `chr4:Group2.1` indicates the 4A166 haplotype, which is not pathogenic when contracted.
  - distal_haplotype: `qAIntactPolyA`, or `qADisruptedPolyA` or `qB`.
  - allele_type: `assembled`, `partial`, `merged` (two partial alleles merged into one) or `assembled_cis_duplication` (in-cis duplications).
  - allele_size: exact or lower bound (for partial or merged alleles)
  - methylation: methylation level

More information on methylation can be found in the `methylation` field.
- `methylation`: Kivvi reports the two fields below for complete alleles (`complete_alleles`), as well as all alleles ends (`all_allele_ends`), i.e. all D4Z4 alleles (including partial alleles) containing the last unit. 
  - `median_methylation_per_unit`: the median methylation level of each repeat unit on each allele (each repeat unit has a value). 
  - `methylation_per_site`: the methylation level of each CpG site (which is the median of values across supporting reads overlapping the CpG site) on each allele. Kivvi uses a set of 101 pre-selected CpG sites in the D4Z4 repeat unit for methylation analysis, so there are `repeat copy number * 101` values reported for each allele.

Other information on D4Z4 can be found in the `additional` field. 
- `qc_metrics`: two metrics indicating the HiFi data quality of the D4Z4 region.
  - `median_read_length`: median read length of the D4Z4 region
  - `per_allele_depth`: median depth per D4Z4 allele (can be different from the genome average depth as D4Z4 is highly GC rich)
- `median_methylation_all_sites`: median methylation level across all D4Z4 units. This value can be used to evaluate the overal methylation level of D4Z4 and is low in FSHD2 samples.
- `median_methylation_per_site`: median methylation of each CpG site
- `median_methylation_per_unit`: median methylation of each repeat unit
- `allele_background`: reports three sets of information:
  - The chromosome and distal haplotype for each complete allele (`complete_alleles`), in the format of `x-y` (`x` can be `qAIntactPolyA`, or `qADisruptedPolyA` or `qB`, and `y` can be `chr4` or `chr10`). 
  - The distal haplotype for all allele ends (`all_allele_ends`), i.e. all D4Z4 alleles (including partial alleles) containing the last unit. 
  - The chromosome for all allele starts (`all_allele_starts`), i.e. all D4Z4 alleles (including partial alleles) containing the first unit. 
- `phasing_of_flanking_regions`: Kivvi uses [Paraphase](https://github.com/PacificBiosciences/paraphase) to phase haplotypes for genomic regions upstream and downstream of D4Z4 (as those regions are also homologous between chr4 and chr10). This field reports the Paraphase phasing results. The upstream haplotypes are used to infer the originating chromosome of D4Z4 alleles. The downstream haplotypes are used to call in-cis duplications of D4Z4 arrays.
- `allele_starts_mapped_to_upstream_haplotypes`: for all D4Z4 alleles (including partial alleles) containing the first unit, report the upstream Paraphase haplotypes that they are connected to. 
- `d4z4_arrays_in_cis`: calls of in-cis duplications of D4Z4 arrays, where multiple D4Z4 alleles are reported as a group indicating that they are next to each other on the same chromosome.

Part of the HG03453 D4Z4 JSON file is shown below. Kivvi assembles three complete alleles and one partial allele. The two alleles (size 28 and 35) on chr4 are both `qAIntactPolyA`. The partial allele is a `qADisruptedPolyA` allele that is at least 14 units long. It's not fully assembled so the chromosome is unknown, but we can infer that it's on chr10 given that 1) `qADisruptedPolyA` are most commonly found on chr10 and 2) among the other three complete alleles in this sample, two are on chr4 and one is on chr10. 
```json
{
  "allele_cn": "28,56,35",
  "complete_alleles": [
    "LeftFlank-1-14-14-13-3-4-12-14-11-11-6-9-9-8-7-8-5-8-10-14-14-8-15-16-18-15-15-17-RightFlank",
    "LeftFlank-54-48-73-51-33-44-38-39-38-38-33-35-38-33-33-38-33-61-59-38-45-36-33-55-56-33-60-55-33-55-52-60-52-60-56-61-59-33-53-55-33-33-55-43-33-33-38-37-55-55-55-56-55-33-33-49-RightFlank",
    "LeftFlank-29-62-32-32-68-69-26-31-24-25-72-27-74-2-58-63-57-30-28-64-65-65-65-66-20-19-20-21-21-19-21-22-21-21-23-RightFlank"
  ],
  "partial_alleles": [
    "41-41-41-41-47-71-41-47-38-41-71-38-38-70-RightFlank"
  ]
  "allele_info": [
    {
      "allele_name": "LeftFlank-1-14-14-13-3-4-12-14-11-11-6-9-9-8-7-8-5-8-10-14-14-8-15-16-18-15-15-17-RightFlank",
      "chromosome": "chr4:Group1.1",
      "distal_haplotype": "qAIntactPolyA",
      "allele_type": "assembled",
      "allele_size": "28",
      "methylation": 0.871
    },
    {
      "allele_name": "LeftFlank-29-62-32-32-68-69-26-31-24-25-72-27-74-2-58-63-57-30-28-64-65-65-65-66-20-19-20-21-21-19-21-22-21-21-23-RightFlank",
      "chromosome": "chr4:Group2.1",
      "distal_haplotype": "qAIntactPolyA",
      "allele_type": "assembled",
      "allele_size": "35",
      "methylation": 0.906
    },
    {
      "allele_name": "LeftFlank-54-48-73-51-33-44-38-39-38-38-33-35-38-33-33-38-33-61-59-38-45-36-33-55-56-33-60-55-33-55-52-60-52-60-56-61-59-33-53-55-33-33-55-43-33-33-38-37-55-55-55-56-55-33-33-49-RightFlank",
      "chromosome": "chr10:Group2.1",
      "distal_haplotype": "qADisruptedPolyA",
      "allele_type": "assembled",
      "allele_size": "56",
      "methylation": 0.9
    },
    {
      "allele_name": "41-41-41-41-47-71-41-47-38-41-71-38-38-70-RightFlank",
      "chromosome": "unknown",
      "distal_haplotype": "qADisruptedPolyA",
      "allele_type": "partial",
      "allele_size": ">14",
      "methylation": 0.9
    }
  ],
}
```
