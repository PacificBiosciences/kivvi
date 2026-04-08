use crate::detail::hapcmp::HapCompare;
use crate::detail::util::{DError, DResult};
use crate::io::json::GeneCall;
use crate::phaser::Phaser;
use crate::realign::reference_length;

use vstr::VStr;

use rand::{prelude::SliceRandom, SeedableRng};
use rust_htslib::bam::{self, record::Cigar, Read, Record};

use std::path::{Path, PathBuf};

pub mod colors {
    pub const READ: &str = "166,206,227";
    pub const READ_ALLELE1: &str = "178,223,138";
    pub const READ_ALLELE2: &str = "177,156,217";
}

pub struct BamWriter<'a> {
    phaser: &'a Phaser,
    call: &'a GeneCall,
}

pub struct IOTuple(pub PathBuf, pub PathBuf, pub String, pub bool);

impl IOTuple {
    #[must_use]
    pub fn source_bam(&self) -> &Path {
        &self.0
    }
    #[must_use]
    pub fn dest_bam(&self) -> &Path {
        &self.1
    }
    #[must_use]
    pub fn chromosome_name(&self) -> &str {
        &self.2
    }
    #[must_use]
    pub fn is_gene2(&self) -> bool {
        self.3
    }
}

impl std::fmt::Debug for IOTuple {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "IOTuple{{Source: {:?}. Dest: {:?}. Chrom: {:?}. Primary or secondary gene: {}",
            self.source_bam(),
            self.dest_bam(),
            self.chromosome_name(),
            if self.is_gene2() {
                "primary"
            } else {
                "secondary"
            }
        )
    }
}

impl<'a> BamWriter<'a> {
    #[must_use]
    pub fn new(phaser: &'a Phaser, call: &'a GeneCall) -> Self {
        Self { phaser, call }
    }

    fn add_tag_to_read(
        &self,
        record: &mut bam::Record,
        use_supp: bool,
        is_gene2: bool,
        rng: Option<&mut rand::rngs::SmallRng>,
    ) -> DResult {
        record.push_aux(b"RN", bam::record::Aux::String(self.phaser.gene_name()))?;
        let haps = &self.call.final_haplotypes;
        if haps.is_empty() {
            record.push_aux(b"HP", bam::record::Aux::String("Unassigned"))?;
            return Ok(());
        }
        let nonunique_reads = &self.call.nonunique_supporting_reads;
        //let alleles = &self.call.alleles_final;
        let mut alleles = Vec::new();
        let alleles_in_call = self.call.region_specific_info.get("raw_alleles");
        if !alleles_in_call.is_none() {
            for allele in alleles_in_call.and_then(|x| x.as_array()).unwrap() {
                let allele_parsed = allele.as_array().unwrap();
                if !allele_parsed.is_empty() {
                    let allele_string = allele_parsed
                        .iter()
                        .map(|x| x.as_str().unwrap().to_string())
                        .collect::<Vec<_>>();
                    alleles.push(allele_string);
                }
            }
        }
        let read_details = &self.call.read_details;
        /*
        hap_found = False
        read_name = read.qname
        if (
            read.is_supplementary
            and self.use_supplementary
            and gene2 is False
            # coordinates are all changed after realigning to gene2
        ):
            read_name = (
                read_name + f"_sup_{read.reference_start}_{read.reference_length}"
            )
        */
        let qname = VStr::from(record.qname()).to_string();
        let qname = if use_supp && !is_gene2 {
            // && record.is_supplementary()
            let ref_length = reference_length(record);
            let mut read_start_pos = 0;
            for x in record.cigar().iter() {
                match x {
                    Cigar::HardClip(_len) | Cigar::SoftClip(_len) => {
                        read_start_pos += i64::from(x.len())
                    }
                    _ => break,
                }
            }
            format!("{qname}_sup_{read_start_pos}_{ref_length}")
        } else {
            qname
        };
        for (hap, hap_name) in &self.call.final_haplotypes {
            if let Some(reads) = self.call.unique_supporting_reads.get(hap) {
                if reads.contains(&qname) {
                    record.push_aux(b"HP", bam::record::Aux::String(hap_name))?;
                    let color = bam::record::Aux::String(if alleles.is_empty() {
                        colors::READ
                    } else if alleles[0].contains(hap_name) {
                        colors::READ_ALLELE1
                    } else if alleles.len() > 1 && alleles[1].contains(hap_name) {
                        colors::READ_ALLELE2
                    } else {
                        colors::READ
                    });
                    record.push_aux(b"YC", color)?;
                    return Ok(());
                }
            }
        }
        if let Some(fingerprint) = read_details.get(&qname) {
            let mut mismatches = Vec::with_capacity(self.call.final_haplotypes.len());
            for (result, hap) in self
                .call
                .final_haplotypes
                .keys()
                .map(|hap| (HapCompare::from_haps(fingerprint, hap), hap))
            {
                let result = result?;
                mismatches.push((hap, result.mismatches));
            }
            mismatches.sort_by(|a, b| a.1.cmp(&b.1));
            if mismatches.len() > 1
                && (1..=2).contains(&mismatches[0].1)
                && mismatches[1].1 >= mismatches[0].1 + 2
            {
                let hp_match = self
                    .call
                    .final_haplotypes
                    .get(mismatches[0].0)
                    .ok_or("key not found")?;
                record.push_aux(b"HP", bam::record::Aux::String(hp_match))?;
                return Ok(());
            }
        }
        /*
        # find closest match
        if read_name in read_details:
            read_seq = read_details[read_name]
            keys = []
            mismatches = []
            for ass_hap in hp_keys.keys():
                match, mismatch, extend = VariantGraph.compare_two_haps(
                    read_seq, ass_hap
                )
                keys.append(ass_hap)
                mismatches.append(mismatch)
            mismatches_sorted = sorted(mismatches)
            if (
                len(mismatches_sorted) > 1
                and 0 < mismatches_sorted[0] <= 2
                and mismatches_sorted[1] >= mismatches_sorted[0] + 2
            ):
                best_match = keys[mismatches.index(mismatches_sorted[0])]
                read.set_tag("HP", hp_keys[best_match], "Z")
                hp_found = True
        */

        let assignment = rng
            .and_then(|rng| {
                nonunique_reads.get(&qname).map(|possible| {
                    &possible
                        .choose(rng)
                        .expect("Expected at least one possibility")[..]
                })
            })
            .unwrap_or("Unassigned");
        if self.call.final_haplotypes.contains_key(assignment) {
            let assignment1 = self
                .call
                .final_haplotypes
                .get(assignment)
                .expect("Check if haplotype exists");
            record.push_aux(b"HP", bam::record::Aux::String(assignment1))?;
        } else {
            record.push_aux(b"HP", bam::record::Aux::String(assignment))?;
        }
        /*
        if rng is False:
            read.set_tag("HP", "Unassigned", "Z")
        else:
            if read_name in nonunique:
                possible_haps = nonunique[read_name]
                random_hap = possible_haps[
                    random.randint(0, len(possible_haps) - 1)
                ]
                if random_hap in hp_keys:
                    read.set_tag("HP", hp_keys[random_hap], "Z")
            else:
                read.set_tag("HP", "Unassigned", "Z")
        */
        Ok(())
    }

    pub fn write_bams(&self) -> Result<Vec<Record>, DError> {
        let mut records_to_write = Vec::new();
        let gene1_input_bam = self.phaser.realigned_bam_path();
        let gene1_output_bam = self.phaser.realigned_tagged_bam_path();
        let chr = self.phaser.chr().unwrap();
        let gene1_inputs = IOTuple(gene1_input_bam, gene1_output_bam, chr.to_string(), false);
        // Assign non-uniquely supporting reads randomly to one haplotype
        let mut records_out: Vec<Record> = self.write_bam(gene1_inputs, Some(1))?;
        records_to_write.append(&mut records_out);
        bam::index::build(
            self.phaser.realigned_tagged_bam_path(),
            None,
            bam::index::Type::Bai,
            1,
        )?;

        Ok(records_to_write)
    }

    /// Attempts to write a bam
    /// # Inputs
    /// 1. `IOTuple` - inbam, outbam, isgene2
    /// 2. Seed - `Option<usize>`. If None, random assign is disabled. Otherwise, seeds a RNG.
    /// # Return
    /// 1. Path to output bam.
    /// 2. bool - whether or not written. If no haplotypes found, no bam is written.
    pub fn write_bam(&self, tuple: IOTuple, seed: Option<u64>) -> Result<Vec<Record>, DError> {
        let mut rng = seed.map(rand::rngs::SmallRng::seed_from_u64);
        let use_supp = self.phaser.use_supplementary();
        log::debug!("IOTuple: {tuple:?}");
        let mut reader = bam::IndexedReader::from_path(tuple.source_bam())?;
        let mut tmp_bam_writer = bam::Writer::from_path(
            tuple.dest_bam(),
            &bam::Header::from_template(reader.header()),
            bam::Format::Bam,
        )?;

        reader.fetch(tuple.chromosome_name())?;
        let mut records = Vec::new();
        let mut record = bam::Record::new();
        while let Some(rc) = reader.read(&mut record) {
            if rc.is_err() {
                continue;
            }
            if record.is_secondary() {
                continue;
            }
            self.add_tag_to_read(&mut record, use_supp, tuple.is_gene2(), rng.as_mut())?;
            records.push(record.clone());
        }
        records.sort_by_cached_key(|x| i64::from(x.tid()) << 32 | x.pos());
        for record in &records {
            tmp_bam_writer.write(record)?;
        }
        Ok(records)
    }
}
