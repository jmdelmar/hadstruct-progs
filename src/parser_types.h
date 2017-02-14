#ifndef _PARSER_TYPES_H
#define _PARSER_TYPES_H 1
#include <qhg.h>
#define MG_MAX_LEVELS 6

struct source_position {
  int coords[ND];
  int nsinks;
  qhg_thrp_nn_sink_params *sinks;
};

struct action_params {
  double mu, csw, kappa;
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
