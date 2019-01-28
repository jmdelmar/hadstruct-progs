#ifndef _PARSER_TYPES_H
#define _PARSER_TYPES_H 1
#include <qhg.h>
#define MG_MAX_LEVELS 6
#define NMOM 12

struct source_position {
  int coords[ND];
  int nmoms;
  int mom_vecs[NMOM][ND-1];
  int nsinks;
  qhg_thrp_nn_sink_params *sinks;
};

struct action_params {
  double mu, mu_l, mu_s, csw, kappa;
  int bc;
};

struct multigrid_params {
  int n_levels;
  int block[ND];
  int setup_iterations[MG_MAX_LEVELS];
  int n_basis_vectors[MG_MAX_LEVELS];
  double coarse_mu[MG_MAX_LEVELS];
  int verbosity;
};

struct smearing_params {
  int n_gauss;
  double alpha_gauss;
  int n_gauss_l;
  double alpha_gauss_l;
  int n_gauss_s;
  double alpha_gauss_s;
  int n_ape;
  double alpha_ape;
};

struct run_params {
  char config[256];
  char prop_dir[256];
  char corr_dir[256];  
  int dims[ND], procs[ND];
  int nsp;
  struct smearing_params smearing;
  struct source_position *spos;
  struct multigrid_params mg;
  struct action_params act;
};

#endif /* _PARSER_TYPES_H */
