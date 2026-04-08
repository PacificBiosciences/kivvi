///
use statrs::distribution::Discrete;
use statrs::distribution::DiscreteCDF;

///
/// Compute probability of genotypes given a set of depths.
///
/// # Panics
/// 1. Could not convert `haploid_depth` to `f64`.
/// 2. Expected `depth < 0`.
/// 3. Assert failure for sanity check.
pub fn depth_prob(nread: i32, haploid_depth: impl std::convert::TryInto<f64>) -> Option<[f32; 4]> {
    let haploid_depth = haploid_depth
        .try_into()
        .unwrap_or_else(|_| panic!("Failed to convert haploid_depth to f64"));
    if haploid_depth <= 0.0 {
        return None;
    }
    let probs = (1..5)
        .map(|x| {
            let expected_depth = f64::from(x) * haploid_depth;
            // scipy: pmf(k, mu, loc=0)
            // k = nread
            // mu = depthexpected
            let poisson =
                statrs::distribution::Poisson::new(expected_depth).expect("Depth should be > 0.");
            poisson.pmf(nread as u64)
        })
        .collect::<Vec<_>>();
    let prob_sum = probs.iter().sum::<f64>();
    if prob_sum == 0. {
        None
    } else {
        let mut res = [0.0f32; 4];
        debug_assert_eq!(res.len(), probs.len());
        probs.iter().zip(res.iter_mut()).for_each(|(pr, res)| {
            *res = (pr / prob_sum) as f32;
        });
        Some(res)
    }
}

/// Compute cdf for heterozygosity.
pub fn het_prob_cdf(total_count: u32, minor_count: u32) -> f64 {
    let distro = statrs::distribution::Binomial::new(0.5, total_count as u64).unwrap();
    distro.cdf(minor_count as u64)
}

// Both constants taken from pharmgoat.
const HET_CDF_CUTOFF: f32 = 0.001;
const HET_MAF_CUTOFF: f32 = 0.05;

/// Computes if a site is probably heterozygous based on on cdf.
pub fn site_probably_het_cdf(total_count: u32, minor_count: u32, cutoff: Option<f32>) -> bool {
    let cdf = het_prob_cdf(total_count, minor_count) as f32;
    cdf >= cutoff.unwrap_or(HET_CDF_CUTOFF)
}

/// Computes if a site is probably heterozygous based on on maf.
pub fn site_probably_het_maf(total_count: u32, minor_count: u32, cutoff: Option<f32>) -> bool {
    (minor_count as f32 / total_count as f32) > cutoff.unwrap_or(HET_MAF_CUTOFF)
}

/// Computes if a site is probably heterozygous based on on maf and cdf.
/// This is the recommended access point, as it is robust in low and high coverage cases.
pub fn site_probably_het(
    total_count: u32,
    minor_count: u32,
    cdf_cutoff: Option<f32>,
    maf_cutoff: Option<f32>,
) -> bool {
    site_probably_het_cdf(total_count, minor_count, cdf_cutoff)
        && site_probably_het_maf(total_count, minor_count, maf_cutoff)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn depth_prob_ok() {
        let prob = depth_prob(40, 20).expect("Failed to generate probs");
        assert_eq!(
            prob.iter()
                .enumerate()
                .max_by(|x, y| x.1.partial_cmp(y.1).unwrap_or(std::cmp::Ordering::Equal))
                .expect("empty iterator")
                .0,
            1
        );
    }
}
