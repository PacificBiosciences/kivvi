use crate::plot::{
    pipe_plot::{Beta, Color, Legend, Pipe, PipePlot, PipeSeg, Shape},
    svg,
};
use crate::util::{invalid_data_error, missing_data_error, DResult};
use crate::variant::{AlleleInfoForPlotting, ReadInfoForPlotting};
use log::trace;
use std::{collections::BTreeMap, path::PathBuf};

pub const FLANK_WIDTH: u32 = 30;

/// Represent a read in the data structure for plotting
#[derive(Clone, Debug, PartialEq)]
struct ReadForPlotting {
    /// start position on the allele (first position without missing info)
    start_position: i64,
    /// bases at variant sites following start_position
    bases: Vec<usize>,
    /// whether each base of read is nonunique
    is_nonuniq: Vec<usize>,
}

/// Arrange reads into the same row
/// # Arguments
/// * `reads` - reads represented as vectors
/// * `spacer` - spacer length
/// # Returns
/// * `Vec<ReadForPlotting>` - reads arranged into the same row
fn partition_reads(reads: Vec<ReadInfoForPlotting>, spacer: usize) -> Vec<ReadForPlotting> {
    let mut new_reads: Vec<ReadForPlotting> = Vec::new();
    for each_read in reads {
        let read_start = each_read.start_position;
        let read = each_read.bases;
        let mut is_nonuniq = Vec::new();
        if each_read.is_nonuniq {
            for _ in 0..read.len() {
                is_nonuniq.push(1);
            }
        } else {
            for _ in 0..read.len() {
                is_nonuniq.push(0);
            }
        }
        if new_reads.is_empty() {
            new_reads.push(ReadForPlotting {
                start_position: read_start,
                bases: read.clone(),
                is_nonuniq,
            });
        } else {
            let mut found_read_to_append = false;
            let mut j_index: usize = 0;
            for (j, existing_read) in new_reads.iter().enumerate() {
                j_index = j;
                let existing_read_end =
                    existing_read.bases.len() as i64 + existing_read.start_position;
                if existing_read_end < read_start {
                    found_read_to_append = true;
                    break;
                }
            }
            if found_read_to_append {
                // add this read to existing read
                let existing_read = &new_reads[j_index];
                let mut existing_read_copy = existing_read.bases.clone();
                let existing_read_end =
                    existing_read.bases.len() as i64 + existing_read.start_position;
                for _ in 0..(read_start - existing_read_end) {
                    existing_read_copy.push(spacer);
                }
                for a in &read {
                    existing_read_copy.push(*a);
                }

                // update is_nonuniq status
                let mut existing_read_nonuniq_copy = existing_read.is_nonuniq.clone();
                for _ in 0..(read_start - existing_read_end) {
                    existing_read_nonuniq_copy.push(0);
                }
                if each_read.is_nonuniq {
                    for _ in 0..read.len() {
                        existing_read_nonuniq_copy.push(1);
                    }
                } else {
                    for _ in 0..read.len() {
                        existing_read_nonuniq_copy.push(0);
                    }
                }

                new_reads[j_index] = ReadForPlotting {
                    start_position: existing_read.start_position,
                    bases: existing_read_copy,
                    is_nonuniq: existing_read_nonuniq_copy,
                };
            } else {
                let mut is_nonuniq = Vec::new();
                if each_read.is_nonuniq {
                    for _ in 0..read.len() {
                        is_nonuniq.push(1);
                    }
                } else {
                    for _ in 0..read.len() {
                        is_nonuniq.push(0);
                    }
                }
                new_reads.push(ReadForPlotting {
                    start_position: read_start,
                    bases: read.clone(),
                    is_nonuniq,
                });
            }
        }
    }
    new_reads
}

/// Plot out each allele and its supporting reads aligned to it
/// # Arguments
/// * `out_file_name` - plot file name
/// * `alleles_for_plot` - allele information
/// # Colors
/// * 0 -> reference, yellow
/// * 1 -> variant, black
/// * 2 -> missing info, pink
/// * 3 -> no data (reads not overlapping), white
/// * 4 -> flank, teal
pub fn plot_alleles_and_reads(
    out_file_name: &PathBuf,
    alleles_for_plot: AlleleInfoForPlotting,
) -> DResult {
    let nvar = alleles_for_plot.variant_count_per_copy;
    let reads = alleles_for_plot.reads;
    // 2 alleles -> 2 vec<Pipe>
    // each pipe is a line/read
    let legend = Legend {
        labels: vec![
            (String::from("Reference"), Color::Yellow),
            (String::from("Variant"), Color::Black),
            (String::from("OtherBase/MissingInfo"), Color::Pink),
            (String::from("Flank"), Color::Teal),
        ],
        height: 4,
    };

    let mut panels: Vec<Vec<Pipe>> = Vec::new();
    for allele_reads in reads {
        let mut panel: Vec<Pipe> = Vec::new();
        // get the scale for the allele
        let mut scale = Vec::new();
        let total_len = allele_reads
            .first()
            .ok_or_else(|| missing_data_error("first allele plot row", "empty allele read set"))?
            .bases
            .len() as i64;
        let cn_max: i64 = total_len / nvar;
        for j in 0..cn_max {
            let this_scale = (FLANK_WIDTH + (j * nvar) as u32, Some((j + 1) as u32));
            scale.push(this_scale);
        }
        let new_reads = partition_reads(allele_reads, 3);
        // update pipes while iterating through reads
        for (i, this_read) in new_reads.iter().enumerate() {
            let mut segs: Vec<PipeSeg> = Vec::new();
            let mut outlines: Vec<PipeSeg> = Vec::new();
            let read_start = this_read.start_position;
            let read = &this_read.bases;
            let read_isnonuniq = &this_read.is_nonuniq;
            // preprare segs
            if *read
                .first()
                .ok_or_else(|| missing_data_error("first plotted read base", format!("{read:?}")))?
                != 4
            {
                segs.push(PipeSeg {
                    width: FLANK_WIDTH - 1 + read_start as u32,
                    color: Color::White,
                    shape: Shape::Rect,
                });
            }
            let mut this_width = 0;
            let mut prev_color = match read
                .first()
                .ok_or_else(|| missing_data_error("first plotted read base", format!("{read:?}")))?
            {
                0 => Color::Yellow,
                1 => Color::Black,
                2 => Color::Pink,
                3 => Color::White,
                4 => Color::Teal,
                _ => Color::LightGray,
            };
            for base in read {
                let this_color = match base {
                    0 => Color::Yellow,
                    1 => Color::Black,
                    2 => Color::Pink,
                    3 => Color::White,
                    4 => Color::Teal,
                    _ => Color::LightGray,
                };
                if this_color == prev_color {
                    this_width += 1;
                } else {
                    let seg_width = if prev_color == Color::Teal {
                        FLANK_WIDTH
                    } else {
                        this_width
                    };
                    segs.push(PipeSeg {
                        width: seg_width,
                        color: prev_color,
                        shape: Shape::Rect,
                    });
                    this_width = 1;
                    prev_color = this_color;
                }
            }
            // last seg
            let seg_width = if prev_color == Color::Teal {
                FLANK_WIDTH
            } else {
                this_width
            };
            segs.push(PipeSeg {
                width: seg_width,
                color: prev_color,
                shape: Shape::Rect,
            });

            // prepare outlines
            let mut this_width = 0;
            let mut prev_color = match read_isnonuniq.first().ok_or_else(|| {
                missing_data_error(
                    "first plotted nonunique marker",
                    format!("{read_isnonuniq:?}"),
                )
            })? {
                0 => Color::White,
                1 => Color::Red,
                _ => Color::LightGray,
            };
            if *read
                .first()
                .ok_or_else(|| missing_data_error("first plotted read base", format!("{read:?}")))?
                != 4
            {
                outlines.push(PipeSeg {
                    width: FLANK_WIDTH - 1 + read_start as u32,
                    color: Color::White,
                    shape: Shape::Rect,
                });
            } else {
                outlines.push(PipeSeg {
                    width: FLANK_WIDTH,
                    color: prev_color.clone(),
                    shape: Shape::Rect,
                });
            }
            for base in read_isnonuniq {
                let this_color = match base {
                    0 => Color::White,
                    1 => Color::Red,
                    _ => Color::LightGray,
                };
                if this_color == prev_color {
                    this_width += 1;
                } else {
                    let seg_width = this_width;
                    outlines.push(PipeSeg {
                        width: seg_width,
                        color: prev_color,
                        shape: Shape::Rect,
                    });
                    this_width = 1;
                    prev_color = this_color;
                }
            }
            // last seg
            if *read
                .last()
                .ok_or_else(|| missing_data_error("last plotted read base", format!("{read:?}")))?
                == 4
            {
                this_width += FLANK_WIDTH - 1;
            }
            let seg_width = this_width;
            outlines.push(PipeSeg {
                width: seg_width,
                color: prev_color,
                shape: Shape::Rect,
            });

            // add pipes
            if i == 0 {
                // first pipe has scale
                let first_pipe = Pipe {
                    segs,
                    betas: vec![],
                    height: 1,
                    outline: vec![],
                    scale: scale.clone(),
                };
                panel.push(first_pipe);
            } else {
                panel.push(Pipe {
                    segs,
                    betas: vec![],
                    height: 1,
                    outline: outlines,
                    scale: vec![],
                });
            }
        }
        panels.push(panel);
    }
    let pipe_plot: PipePlot = PipePlot { panels, legend };
    svg::generate(
        &pipe_plot,
        out_file_name.to_str().ok_or_else(|| {
            invalid_data_error(format!(
                "Plot output path is not valid UTF-8: {}",
                out_file_name.display()
            ))
        })?,
    );
    Ok(())
}

/// Plot methylation levels
/// # Arguments
/// * `out_file_name` - plot file name
/// * `allele_methyl` - allele -> read -> methylation levels
pub fn plot_methyl(
    out_file_name: &PathBuf,
    allele_methyl: Vec<BTreeMap<String, Vec<usize>>>,
) -> DResult {
    let mut panels: Vec<Vec<Pipe>> = Vec::new();
    let mut allele_index = 0;
    for each_allele in allele_methyl {
        allele_index += 1;
        /*
        let allele_fp_methyl_values = allele_methyl
            .allele_fps_methyl_value
            .get(&allele_index)
            .ok_or("err")?;
        */
        let mut panel: Vec<Pipe> = Vec::new();
        let mut each_allele_for_plot = Vec::new();
        let mut read_index = 0;
        for (read, read_methyl) in &each_allele {
            //assert_eq!(allele_fp_methyl_values.len(), read_methyl.len());
            read_index += 1;
            // plot allele
            if read_index == 1 {
                let mut segs: Vec<PipeSeg> = Vec::new();
                let mut scale = Vec::new();
                let nvar: i64 = 101;
                let total_len: i64 = read_methyl.len() as i64;
                let cn_max: i64 = total_len / nvar;
                for j in 0..cn_max {
                    let this_scale = ((j * nvar) as u32, Some((j + 1) as u32));
                    scale.push(this_scale);
                }
                /*
                let mut betas: Vec<Beta> = Vec::new();
                for i in 0..total_len {
                    let beta = Beta {
                        pos: i as usize,
                        value: allele_fp_methyl_values[i],
                    };
                    betas.push(beta);
                }
                */
                segs.push(PipeSeg {
                    width: total_len as u32,
                    color: Color::Gray,
                    shape: Shape::Rect,
                });
                panel.push(Pipe {
                    segs,
                    betas: vec![],
                    height: 1,
                    outline: vec![],
                    scale: scale.clone(),
                });
            }

            let mut beginning_unknown: usize = 0;
            for a in read_methyl {
                if *a == 400 {
                    beginning_unknown += 1;
                } else {
                    break;
                }
            }
            let mut bases = read_methyl[beginning_unknown..]
                .into_iter()
                .map(|x| *x as usize)
                .collect::<Vec<usize>>();
            if !bases.is_empty() {
                while *bases.last().ok_or_else(|| {
                    missing_data_error(
                        "last methylation base in plotted read",
                        format!("{bases:?}"),
                    )
                })? == 400
                {
                    bases.pop();
                }
                let this_read = ReadInfoForPlotting {
                    start_position: beginning_unknown as i64,
                    bases: bases.clone(),
                    is_nonuniq: false,
                };
                each_allele_for_plot.push(this_read.clone());
                trace!(
                    "each_allele_for_plot allele_index {allele_index} read {read} {:?}",
                    this_read
                );
            }
        }
        each_allele_for_plot.sort_by(|a, b| a.start_position.cmp(&b.start_position));
        let new_reads = partition_reads(each_allele_for_plot.clone(), 300);
        trace!("new_reads {:?}", new_reads);

        for this_read in new_reads.iter() {
            let this_read_start = this_read.start_position;
            let mut segs: Vec<PipeSeg> = Vec::new();
            let mut betas: Vec<Beta> = Vec::new();
            segs.push(PipeSeg {
                width: this_read_start as u32,
                color: Color::White,
                shape: Shape::Rect,
            });
            let bases = &this_read.bases;
            for pos in 0..bases.len() {
                let this_base = bases[pos as usize];
                if this_base != 300 {
                    segs.push(PipeSeg {
                        width: 1,
                        color: Color::Gray,
                        shape: Shape::Rect,
                    });
                    if this_base < 256 {
                        let beta = Beta {
                            pos: pos + this_read_start as usize,
                            value: (this_base as f64) / 255.0,
                        };
                        betas.push(beta);
                    }
                } else {
                    segs.push(PipeSeg {
                        width: 1,
                        color: Color::White,
                        shape: Shape::Rect,
                    });
                }
            }
            panel.push(Pipe {
                segs,
                betas,
                height: 1,
                outline: vec![],
                scale: vec![],
            });
        }
        panels.push(panel);
    }
    let legend = Legend {
        labels: vec![
            ("Methylated".to_string(), Color::Grad(1.0)),
            ("Unmethylated".to_string(), Color::Grad(0.0)),
        ],
        height: 4,
    };
    let pipe_plot: PipePlot = PipePlot { panels, legend };
    svg::generate(
        &pipe_plot,
        out_file_name.to_str().ok_or_else(|| {
            invalid_data_error(format!(
                "Plot output path is not valid UTF-8: {}",
                out_file_name.display()
            ))
        })?,
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::env;

    #[test]
    fn test_partition_reads() {
        // overlapping
        let mut reads = Vec::new();
        reads.push(ReadInfoForPlotting {
            start_position: 1,
            bases: vec![1, 1, 1],
            is_nonuniq: false,
        });
        reads.push(ReadInfoForPlotting {
            start_position: 2,
            bases: vec![1, 1, 1],
            is_nonuniq: true,
        });
        let new_reads = partition_reads(reads.clone(), 3);
        assert_eq!(
            new_reads,
            vec![
                ReadForPlotting {
                    start_position: 1,
                    bases: vec![1, 1, 1],
                    is_nonuniq: vec![0, 0, 0],
                },
                ReadForPlotting {
                    start_position: 2,
                    bases: vec![1, 1, 1],
                    is_nonuniq: vec![1, 1, 1],
                }
            ]
        );

        // right next to each other
        let mut reads = Vec::new();
        reads.push(ReadInfoForPlotting {
            start_position: 1,
            bases: vec![1, 1, 1],
            is_nonuniq: true,
        });
        reads.push(ReadInfoForPlotting {
            start_position: 4,
            bases: vec![1, 1, 1],
            is_nonuniq: false,
        });
        let new_reads = partition_reads(reads.clone(), 3);
        assert_eq!(
            new_reads,
            vec![
                ReadForPlotting {
                    start_position: 1,
                    bases: vec![1, 1, 1],
                    is_nonuniq: vec![1, 1, 1],
                },
                ReadForPlotting {
                    start_position: 4,
                    bases: vec![1, 1, 1],
                    is_nonuniq: vec![0, 0, 0],
                }
            ]
        );

        // not touching
        let mut reads = Vec::new();
        reads.push(ReadInfoForPlotting {
            start_position: 1,
            bases: vec![1, 1, 1],
            is_nonuniq: false,
        });
        reads.push(ReadInfoForPlotting {
            start_position: 5,
            bases: vec![1, 1, 1],
            is_nonuniq: true,
        });
        let new_reads = partition_reads(reads.clone(), 3);
        assert_eq!(
            new_reads,
            vec![ReadForPlotting {
                start_position: 1,
                bases: vec![1, 1, 1, 3, 1, 1, 1],
                is_nonuniq: vec![0, 0, 0, 0, 1, 1, 1],
            }]
        );
    }

    #[test]
    fn test_plot_alleles_and_reads_errors_on_empty_allele_read_set() {
        let error = plot_alleles_and_reads(
            &env::temp_dir().join("kivvi-empty-allele-plot.svg"),
            AlleleInfoForPlotting {
                reads: vec![vec![]],
                variant_count_per_copy: 1,
            },
        )
        .expect_err("empty allele plot rows should error");

        assert!(
            error
                .to_string()
                .contains("missing first allele plot row: empty allele read set"),
            "unexpected error: {error}"
        );
    }
}
