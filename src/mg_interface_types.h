#ifndef _MG_INTERFACE_TYPES_H
#define _MG_INTERFACE_TYPES_H 1
#include <DDalphaAMG.h>
#define MG_MAX_LEVELS 6

enum mu_sign {
  plus,
  minus
};

typedef struct {
  DDalphaAMG_init init;
  DDalphaAMG_parameters params;
  DDalphaAMG_status status;
  enum mu_sign current_mu_sign;
} mg_state;

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

#endif /* _MG_INTERFACE_TYPES_H */
