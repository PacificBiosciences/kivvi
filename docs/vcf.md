# Kivvi VCF

Kivvi produces a VCF file for each region for each sample. The Kivvi VCF follows the standard VCF format specifications. The QUAL field contains a placeholder value. We have repurposed the sample column so that each column represents a completely assembled allele. In addition, we report the `RU` INFO field to represent the repeat units where a variant is found (repeat units are separated by comma). The four values for each repeat unit are:
1. Repeat unit ID
2. Where to find this repeat unit, represented as `x.y`, where `x` is the index of the allele and `y` is the index of the repeat unit on the allele.
3. Read depth at this position for this repeat unit
4. Number of reads supporting the variant

## Example KIV2 VCF

| CHROM | POS       | ID | REF | ALT | QUAL | FILTER | INFO                         | FORMAT   | allele1 | allele2 |
|:------| :-------- | :- | :-- | :-- | :--- | :----- | :----------------------------| :--------| :-------| :-------|
| chr6  | 160613775 | .  | G   | C   | .    | PASS   | RU=31:2.3:31:31,30:2.4:30:30 | GT       | 0       | 1       |

This variant is found on two repeat units:
1. Repeat unit index #31, which is the third unit on the second allele (`2.3`). There are 31 reads covering this position for this repeat unit and 31 of those support the variant.
2. Repeat unit index #30, which is the fourth unit on the second allele (`2.4`). There are 30 reads covering this position for this repeat unit and 30 of those support the variant.

This variant is only found on the second allele, so its GT is `0` for allele1 and `1` for allele2.

## Example D4Z4 VCF

| CHROM    | POS  | ID | REF | ALT | QUAL | FILTER | INFO                         | FORMAT   | allele1 | allele2 | allele3 |
|:---------| :--- | :- | :-- | :-- | :--- | :----- | :----------------------------| :--------| :-------| :-------| :-------|
| d4z4_ref | 1611 | .  | C   | G   | .    | PASS   | RU=73:2.3:19:19,44:2.6:21:21 | GT       | 0       | 1       | 0       |

D4Z4 reference can be found [here](../data/d4z4/d4z4_ref.fa).
