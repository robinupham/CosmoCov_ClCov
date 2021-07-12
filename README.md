# CosmoCov_ClCov - shear harmonic space covariance

Fork of [CosmoLike/CosmoCov](https://github.com/CosmoLike/CosmoCov) with interface for shear harmonic space covariance.

**[CosmoLike/CosmoCov](https://github.com/CosmoLike/CosmoCov) is written by Xiao Fang, Elisabeth Krause & Tim Eifler. Please see [their list of papers to cite here](https://github.com/CosmoLike/CosmoCov#papers-to-cite).**

## Compilation
```shell
$ make shear_clcov
```
The GSL and FFTW3 libraries are required.

## Usage
The code is designed to produce one block-row of the full shear covariance matrix in a single command, where a single block corresponds to a pair of power spectra.

One block-row corresponds to holding one power spectrum of the pair fixed. This is the spectrum indexed by `spec1_idx` (see [Power spectrum indexing](#power-spectrum-indexing) section below). The code will automatically iterate over all `spec2_idx <= spec1_idx`.

For a given `spec1_idx`, the code can be run as follows:

```shell
$ ./get_shear_clcov {path_to_config.ini} {spec1_idx}
```
where `{path_to_config.ini}` is the path to the configuration file (see [Configuration](#configuration) section below).

## Power spectrum indexing

Since each power spectrum corresponds to a pair of shear fields, the set of power spectra may be laid out in the upper (or lower) triangle of a matrix.

The convention used here is to order the power spectra by diagonal of this matrix, then by row. Note that this corresponds to `new=True` ordering in healpy. An example is shown below for 5 shear fields, labelled `bin1` to `bin5`.

Zero-based indexing is used.

```
      bin1 bin2 bin3 bin4 bin5
bin1     0    5    9   12   14
bin2     -    1    6   10   13
bin3     -    -    2    7   11
bin4     -    -    -    3    8
bin5     -    -    -    -    4
```

## Configuration
See an example configuration file: [example_input/example.ini](example_input/example.ini).

The config file contains all the settings other than `{spec1_idx}`. Most are the same as the original CosmoCov settings. A complete list of settings is below.

Inherited from CosmoCov:

- `Omega_m`, `Omega_v`, `omb`, `sigma_8`, `n_spec`, `w0` `wa`, `h0`: cosmological parameters (Omega_v == Omega_lambda);

- `area`: survey area in square degrees;

- `c_footprint_file` : (optional) mask power spectrum used for super-sample covariance (is automatically normalised);

- `clustering_REDSHIFT_FILE`, `shear_REDSHIFT_FILE`, `lens_tomobins`, `source_tomobins`, `lens_n_gal`, `source_n_gal`: details of lens and source galaxy samples (file paths, the numbers of tomographic bins, the number densities in each bin); the redshift file has (number of tomo bin + 1) columns, in which the 1st column is the z_min of each z bin;

- `sigma_e`: total intrinsic shape dispersion;

- `lens_tomogbias`: linear galaxy bias parameter of each lens galaxy bin;

- `lens_tomo_bmag`: magnification bias parameter of each lens galaxy bin (with `b_mag` described in Section 5.1.3 of [Fang et al. (arXiv:1911.11947)](https://arxiv.org/abs/1911.11947));

- `IA`: 0 or 1, the switch of running the intrinsic alignment NLA model;

- `A_ia`, `eta_ia`:  parameters of the NLA model (see Eq. 4.9 of [Fang et al. (arXiv:1911.11947)](https://arxiv.org/abs/1911.11947), but with `A_ia` represented by `a_IA` in the equation).

ClCov-specific settings:

- `lmin`, `lmax`: continuous ell range required;

- `ell`: alternatively, specify discrete ells separated by commas, e.g. `2,5,10,100`;

- `do_g`: include Gaussian contribution (which includes shape noise);

- `do_ss`: include super-sample covariance;

- `do_cng`: include connected non-Gaussian covariance (slow and generally sub-dominant);

- `output_dir`: directory to output results to (must already exist).

## Output

For each (`spec1`, `spec2`) pair a text file is produced, with filename `cov_{contributions}_spec1_{spec1_idx}_spec2_{spec2_idx}.txt`, where `{contributions}` is the relevant combination of `g`, `ss` and `cng` separated by underscores, depending on the settings of `do_g`, `do_ss` and `do_cng`.

The file contains the covariance matrix of the two power spectra with indices `spec1_idx` and `spec2_idx`. It is symmetric if `spec2_idx == spec1_idx`.
