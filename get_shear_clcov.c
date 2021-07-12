/*************************************************************************************************************
Calculate shear Gaussian + non-Gaussian covariance using CosmoLike.
Gaussian, super-sample and connected non-Gaussian contributions can each be switched on and off independently.
This loops over all combinations of z bins, saving each to disk separately.

Provide two command line arguments:
  1. Path to ini file with all the configuration (see the example_input/example.ini).
  2. Index of the first power spectrum, spec1_idx. This will then iterate over all spec2_idx <= spec1_idx.
**************************************************************************************************************/

#include <string.h>

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_expint.h>

// CosmoLike includes (std & GSL includes must be above here)
#include <cosmolike_core/theory/structs.c>
#include <cosmolike_core/theory/basics.c>
#include <cosmolike_core/theory/recompute.c>
#include <cosmolike_core/theory/cosmo3D.c>
#include <cosmolike_core/theory/redshift_spline.c>
#include <cosmolike_core/theory/cosmo2D_fourier.c>
#include <cosmolike_core/theory/halo.c>
#include <cosmolike_core/theory/covariances_3D.c>
#include <cosmolike_core/theory/covariances_fourier.c>
#include <covs/init.c>


// Cl covariance params
typedef struct {
    int lmin;
    int lmax;
    char ell_str[200];
    char output_dir[200];
    int do_g;
    int do_ss;
    int do_cng;
} clcov_par;
clcov_par clcov_params = {.lmin = 2, .lmax = -1, .ell_str = "NULL", .do_g = 0, .do_ss = 0, .do_cng = 0};


// Set Cl covariance params from ini file
void set_clcov_parameters(char *paramfile, int output)
{
  char line[256];
  int iline=0;

  FILE* input = fopen(paramfile, "r");
  while(fgets(line, 256, input) != NULL) {
    char name[128], val[128];
    iline++;
    if(line[0] == '#') continue;
    sscanf(line, "%128s : %128s", name, val);

    if(strcmp(name, "lmin") == 0) {
      sscanf(val, "%d", &clcov_params.lmin);
      if(output == 1) printf("lmin %d \n", clcov_params.lmin);
    } else if(strcmp(name, "lmax") == 0) {
      sscanf(val, "%d", &clcov_params.lmax);
      if(output == 1) printf("lmax %d \n", clcov_params.lmax);
    } else if(strcmp(name, "output_dir") == 0) {
      sprintf(clcov_params.output_dir, "%s", val);
      if(output == 1) printf("output_dir %s \n", clcov_params.output_dir);
    } else if(strcmp(name, "c_footprint_file") == 0) {
      sprintf(covparams.C_FOOTPRINT_FILE, "%s", val);
      if(output == 1) printf("c_footprint_file %s \n", covparams.C_FOOTPRINT_FILE);
    } else if(strcmp(name, "ell") == 0) {
      sprintf(clcov_params.ell_str, "%s", val);
      if(output == 1) printf("ell_str %s \n", clcov_params.ell_str);
    } else if(strcmp(name, "do_g") == 0) {
      sscanf(val, "%d", &clcov_params.do_g);
      if(output == 1) printf("do_g %d \n", clcov_params.do_g);
    } else if(strcmp(name, "do_ss") == 0) {
      sscanf(val, "%d", &clcov_params.do_ss);
      if(output == 1) printf("do_ss %d \n", clcov_params.do_ss);
    } else if(strcmp(name, "do_cng") == 0) {
      sscanf(val, "%d", &clcov_params.do_cng);
      if(output == 1) printf("do_cng %d \n", clcov_params.do_cng);
    }
  }
}


// 2D version of numpy savetxt. Use header = "0" for no header
int savetxt_2d(char *path, double **array, int n_row, int n_col, char *header) {

    FILE *file = fopen(path, "w");

    if (strncmp(header, "0", 2) != 0) fprintf(file, "# %s\n", header);

    for (int i = 0; i < n_row; i++) {
      printf("Saving %d / %d ... \r", i + 1, n_row);
      for (int j = 0; j < n_col; j++) {
        fprintf(file, "%.10e", array[i][j]);
        if (j < n_col - 1) {
          fprintf(file, " ");
        }
      }
      fprintf(file, "\n");
    }

    fclose(file);
    return 0;
}


int main(int argc, char** argv)
{

  if (argc != 3){
    fprintf(stderr, "Syntax: %s config_file spec1_idx\n", argv[0]);
    exit(1);
  }

  setbuf(stdout, NULL);

  // Set this to 1 to output details about inputs for diagnostics
  int output = 1;

  // Set cosmo, survey and Cl covariance params from ini file
  char *inifile = argv[1];
  set_cosmological_parameters(inifile, output);
  set_survey_parameters(inifile, output);
  covparams.full_tomo = 1;
  set_clcov_parameters(inifile, output);
  init_source_sample(redshift.shear_REDSHIFT_FILE, tomo.shear_Nbin);
  init_lens_sample(redshift.clustering_REDSHIFT_FILE, tomo.clustering_Nbin);

  // Determine which contributions to include
  char cov_str[9] = "";
  printf("Including Gaussian contribution: ");
  if (clcov_params.do_g) {
    printf("yes \n");
    sprintf(cov_str, "%sg_", cov_str);
  } else printf("no \n");
  printf("Including super-sample contribution: ");
  if (clcov_params.do_ss) {
    printf("yes \n");
    sprintf(cov_str, "%sss_", cov_str);
  } else printf("no \n");
  printf("Including connected non-Gaussian contribution: ");
  if (clcov_params.do_cng) {
    printf("yes \n");
    sprintf(cov_str, "%scng_", cov_str);
  } else printf("no \n");
  covparams.ssc = clcov_params.do_ss;
  covparams.cng = clcov_params.do_cng;

  // Obtain the number of z bins from config
  int n_zbin;
  if (tomo.shear_Nbin != tomo.clustering_Nbin) error("tomo.shear_Nbin != tomo.clustering_Nbin");
  else n_zbin = tomo.shear_Nbin;

  // Determine ells, either from the list provided in the config or the provided lmin and lmax
  // Do this in 2 steps: first determine the number of ells, then the ells themselves
  int n_ell;
  if (strcmp(clcov_params.ell_str, "NULL") == 0) {
    if (clcov_params.lmax < clcov_params.lmin) error("lmax < lmin, make sure lmax is provided");
    n_ell = clcov_params.lmax - clcov_params.lmin + 1;
  }
  else { // Count the number of commas + 1
    n_ell = 1;
    const char *remaining_str = clcov_params.ell_str;
    while ((remaining_str = strstr(remaining_str, ",")))
    {
      n_ell++;
      remaining_str++;
    }
  }
  int ell[n_ell];
  if (strcmp(clcov_params.ell_str, "NULL") == 0) {
    for (int l = clcov_params.lmin; l <= clcov_params.lmax; l++) {
      ell[l - clcov_params.lmin] = l;
    }
  }
  else {
    char* rem_str = calloc(strlen(clcov_params.ell_str) + 1, sizeof(char));
    strcpy(rem_str, clcov_params.ell_str);
    int l_idx = 0;
    char *token = strtok(rem_str, ",");
    while (token != NULL) {
      ell[l_idx] = strtol(token, NULL, 10);
      token = strtok(NULL, ",");
      l_idx++;
    }
  }

  // Prepare a string about the ells to go in the header of the output text files
  char ell_header_str[200];
  if (strcmp(clcov_params.ell_str, "NULL") == 0) {
    sprintf(ell_header_str, "lmin %d, lmax %d", clcov_params.lmin, clcov_params.lmax);
  } else {
    sprintf(ell_header_str, "ell: %s", clcov_params.ell_str);
  }
  printf("%s \n", ell_header_str);

  // Generate list of power spectra within the shear data vector
  // Rows are the different spectra, cols are: z1, z2
  int n_spectra_sametype = n_zbin * (n_zbin + 1) / 2;
  int spectra_sametype[n_spectra_sametype][2];
  int spec_idx = 0;
  for (int diag = 0; diag < n_zbin; diag++) {
    for (int row = 0; row < n_zbin - diag; row++) {
      spectra_sametype[spec_idx][0] = row;
      spectra_sametype[spec_idx][1] = row + diag;
      spec_idx++;
    }
  }

  // Read in the provided spec1_idx and check it is allowed
  int spec1_idx = atoi(argv[2]);
  if (spec1_idx < 0 || spec1_idx >= n_spectra_sametype) {
    error("spec1_idx provided is not valid");
  } else {
    printf("Requested spec1_idx is %d \n", spec1_idx);
  }

  // Common variables for the covariance blocks
  int z1, z2, z3, z4, l1, l2;
  double cov;
  char save_path[250];
  char header[250];

  // Shear-shear covariance (symmetric so only do spec2_idx <= spec1_idx)
  for (int spec2_idx = 0; spec2_idx <= spec1_idx; spec2_idx++) {

    // Allocate memory for the block and fill with nans
    double** cov_shear_shear = (double**) malloc(n_ell * sizeof(double*));
    for (int i = 0; i < n_ell; i++) cov_shear_shear[i] = (double*) malloc(n_ell * sizeof(double));
    for (int l1_idx = 0; l1_idx < n_ell; l1_idx++) {
      for (int l2_idx = 0; l2_idx < n_ell; l2_idx++) {
        cov_shear_shear[l1_idx][l2_idx] = NAN;
      }
    }

    // Loop over all combinations of (l1, l2), unless the two spectra are the same in which case only need l2 <= l1
    for (int l1_idx = 0; l1_idx < n_ell; l1_idx++) {
      printf("spec2_idx %d / %d, l1_idx %d / %d \n", spec2_idx, spec1_idx, l1_idx, n_ell);
      int l2_idx_max;
      if (spec2_idx == spec1_idx) {
        l2_idx_max = l1_idx;
      } else {
        l2_idx_max = n_ell - 1;
      }
      for (int l2_idx = 0; l2_idx <= l2_idx_max; l2_idx++) {
        l1 = ell[l1_idx];
        l2 = ell[l2_idx];
        z1 = spectra_sametype[spec1_idx][0];
        z2 = spectra_sametype[spec1_idx][1];
        z3 = spectra_sametype[spec2_idx][0];
        z4 = spectra_sametype[spec2_idx][1];
        cov = 0.0;
        if (clcov_params.do_g) cov += cov_G_shear_shear_tomo(l1, l2, z1, z2, z3, z4);
        if (clcov_params.do_ss || clcov_params.do_cng) cov += cov_NG_shear_shear_tomo(l1, l2, z1, z2, z3, z4);
        cov_shear_shear[l1_idx][l2_idx] = cov;
        if (spec2_idx == spec1_idx) {
          cov_shear_shear[l2_idx][l1_idx] = cov;
        }
      }
    }

    // Save block to disk and free memory
    sprintf(save_path, "%s/cov_%sspec1_%d_spec2_%d.txt", clcov_params.output_dir, cov_str, spec1_idx, spec2_idx);
    sprintf(header, "cov_shear_shear for spec_idx %d and %d with do_g=%d, do_ss=%d, do_cng=%d, %s",
            spec1_idx, spec2_idx, clcov_params.do_g, clcov_params.do_ss, clcov_params.do_cng, ell_header_str);
    if (savetxt_2d(save_path, cov_shear_shear, n_ell, n_ell, header) == 0) printf("\nSaved %s\n", save_path);
    for (int i = 0; i < n_ell; i++) free(cov_shear_shear[i]);
    free(cov_shear_shear);

  } // next spec2_idx for this fixed spec1_idx

  printf("Done \n");
  return 0;
}
