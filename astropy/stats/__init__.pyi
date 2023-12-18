# Licensed under a 3-clause BSD style license - see LICENSE.rst
from .bayesian_blocks import (
    Events as Events,
    FitnessFunc as FitnessFunc,
    PointMeasures as PointMeasures,
    RegularEvents as RegularEvents,
    bayesian_blocks as bayesian_blocks,
)
from .biweight import (
    biweight_location as biweight_location,
    biweight_midcorrelation as biweight_midcorrelation,
    biweight_midcovariance as biweight_midcovariance,
    biweight_midvariance as biweight_midvariance,
    biweight_scale as biweight_scale,
)
from .circstats import (
    circcorrcoef as circcorrcoef,
    circmean as circmean,
    circmoment as circmoment,
    circstd as circstd,
    circvar as circvar,
    rayleightest as rayleightest,
    vonmisesmle as vonmisesmle,
    vtest as vtest,
)
from .funcs import (
    binned_binom_proportion as binned_binom_proportion,
    binom_conf_interval as binom_conf_interval,
    bootstrap as bootstrap,
    cdf_from_intervals as cdf_from_intervals,
    fold_intervals as fold_intervals,
    gaussian_fwhm_to_sigma as gaussian_fwhm_to_sigma,
    gaussian_sigma_to_fwhm as gaussian_sigma_to_fwhm,
    histogram_intervals as histogram_intervals,
    interval_overlap_length as interval_overlap_length,
    kuiper as kuiper,
    kuiper_false_positive_probability as kuiper_false_positive_probability,
    kuiper_two as kuiper_two,
    mad_std as mad_std,
    median_absolute_deviation as median_absolute_deviation,
    poisson_conf_interval as poisson_conf_interval,
    signal_to_noise_oir_ccd as signal_to_noise_oir_ccd,
)
from ._histogram import (
    calculate_bin_edges as calculate_bin_edges,
    freedman_bin_width as freedman_bin_width,
    histogram as histogram,
    knuth_bin_width as knuth_bin_width,
    scott_bin_width as scott_bin_width,
)
from .info_theory import (
    akaike_info_criterion as akaike_info_criterion,
    akaike_info_criterion_lsq as akaike_info_criterion_lsq,
    bayesian_info_criterion as bayesian_info_criterion,
    bayesian_info_criterion_lsq as bayesian_info_criterion_lsq,
)
from .jackknife import (
    jackknife_resampling as jackknife_resampling,
    jackknife_stats as jackknife_stats,
)
from .sigma_clipping import (
    SigmaClip as SigmaClip,
    sigma_clip as sigma_clip,
    sigma_clipped_stats as sigma_clipped_stats,
)
from .spatial import RipleysKEstimator as RipleysKEstimator
