use crate::util::{DError, RegionCoordinates};
use itertools::Itertools;
use log::{debug, info, warn};
use rust_htslib::bam::{self, IndexedReader, Read};
use std::path::PathBuf;

/// Summary of depth information
#[derive(Clone, Debug, Default, serde::Serialize, serde::Deserialize)]
pub struct DepthSummary {
    /// genome average depth
    pub genome_depth: Option<f32>,
    /// repeat average depth
    pub repeat_depth: Option<f32>,
    /*
    /// depth-based estimation of the copy number of the repeat
    pub depth_based_cn: Option<f32>,
    /// depth-based estimation of the copy number of the repeat, adjusted
    pub depth_based_cn_adjusted: Option<f32>,
    */
}

/// Calculated depth-bassed copy number
/// # Arguments
/// * `wgs_bam` - WGS bam, used for genome average depth
/// * `realigned_bam` - realigned bam, used for repeat depth
/// * `region_coordinates` - coordinated defined for this region
/// # Returns
/// * `DepthSummary` - summary of depth information
pub fn depth_based_cn(
    wgs_bam: PathBuf,
    realigned_bam: PathBuf,
    region_coordinates: RegionCoordinates,
) -> Result<DepthSummary, DError> {
    let mut gdepth: Vec<i32> = Vec::new();
    let mut rdepth: Vec<i32> = Vec::new();
    // genome depth
    let mut bam_reader = IndexedReader::from_path(wgs_bam)?;
    let depth_chrom = if bam_reader.header().tid(b"chr6").is_some() {
        "chr6"
    } else {
        "6"
    };
    for pos in region_coordinates.genome_depth_sites {
        bam_reader
            .fetch(bam::FetchDefinition::RegionString(
                depth_chrom.as_bytes(),
                pos + 500,
                pos + 1500,
            ))
            .unwrap_or_else(|e| panic!("Failed to fetch region {e}"));
        for p in bam_reader.pileup() {
            let pileup = p?;
            let this_depth = pileup.depth();
            let this_pos = pileup.pos() as i64;
            if this_pos < pos + 1500 && this_pos > pos + 500 {
                gdepth.push(this_depth as i32);
            }
        }
    }
    debug!("gdepth {:?} sites", gdepth.len());
    // repeat depth
    let mut bam_reader = IndexedReader::from_path(realigned_bam)?;
    let (nstart, nend) = region_coordinates.depth_regions;
    for p in bam_reader.pileup() {
        let pileup = p?;
        let this_depth = pileup.depth();
        let this_pos = pileup.pos() as i64;
        if this_pos <= nend && this_pos >= nstart {
            rdepth.push(this_depth as i32);
        }
    }
    debug!("rdepth {:?} sites", rdepth.len());
    let genome_median = median(&gdepth);
    let repeat_median = median(&rdepth);

    if let Some(genome_median_value) = genome_median {
        info!("Genome depth is {:?}", genome_median.unwrap());
        if genome_median_value < 15.0 {
            warn!("Genome depth is low. Recommend sequencing to a higher coverage (>20X).");
        } else if let Some(repeat_median_value) = repeat_median {
            // CN estimate is too noisy. Disabled.
            let _repeat_cn = (2.0 * repeat_median_value / genome_median_value).round();
            let _repeat_cn_adjusted = (_repeat_cn * 0.7 + 0.6).round();
        }
    } else {
        warn!("Genome depth is unavailable.");
    }
    Ok(DepthSummary {
        genome_depth: genome_median,
        repeat_depth: repeat_median,
    })
}

/// Compute the median of a slice of integers.
/// Assigns a float, averages midpoints in case of an even-sized slice.
/// ```
/// use kivvi::depth::median;
/// let x = &[0, 3i32, 4];
/// assert_eq!(3., median(x).unwrap());
/// assert!(median(&[] as &[i32]).is_none());
/// let x = &[4i32, 2];
/// assert_eq!(median(x).unwrap(), 3.);
/// let x = &[1i32, 2];
/// assert_eq!(median(x).unwrap(), 1.5f32);
/// ```
#[must_use]
pub fn median(x: &[i32]) -> Option<f32> {
    if x.is_empty() {
        return None;
    }
    let copied = x.iter().copied().sorted().collect::<Vec<_>>();
    let len = copied.len();
    let midpoint = len / 2;
    if len & 1 != 0 {
        Some(copied[midpoint] as f32)
    } else {
        Some((copied[midpoint] + copied[midpoint - 1]) as f32 * 0.5)
    }
}

/*
// From https://rust-lang-nursery.github.io/rust-cookbook/science/mathematics/statistics.html
// Functions below are used to calculate median
use std::cmp::Ordering;

fn partition(data: &[i32]) -> Option<(Vec<i32>, i32, Vec<i32>)> {
    match data.len() {
        0 => None,
        _ => {
            let (pivot_slice, tail) = data.split_at(1);
            let pivot = pivot_slice[0];
            let (left, right) = tail.iter().fold((vec![], vec![]), |mut splits, next| {
                {
                    let (ref mut left, ref mut right) = &mut splits;
                    if next < &pivot {
                        left.push(*next);
                    } else {
                        right.push(*next);
                    }
                }
                splits
            });

            Some((left, pivot, right))
        }
    }
}

fn select(data: &[i32], k: usize) -> Option<i32> {
    let part = partition(data);

    match part {
        None => None,
        Some((left, pivot, right)) => {
            let pivot_idx = left.len();

            match pivot_idx.cmp(&k) {
                Ordering::Equal => Some(pivot),
                Ordering::Greater => select(&left, k),
                Ordering::Less => select(&right, k - (pivot_idx + 1)),
            }
        }
    }
}

pub fn median(data: &[i32]) -> Option<f32> {
    let size = data.len();
    if data.is_empty() {
        return None;
    }

    match size {
        even if even % 2 == 0 => {
            let fst_med = select(data, (even / 2) - 1);
            let snd_med = select(data, even / 2);

            match (fst_med, snd_med) {
                (Some(fst), Some(snd)) => Some((fst + snd) as f32 / 2.0),
                _ => None,
            }
        }
        odd => select(data, odd / 2).map(|x| x as f32),
    }
}
*/
