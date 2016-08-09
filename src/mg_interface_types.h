#ifndef _MG_INTERFACE_TYPES_H
#define _MG_INTERFACE_TYPES_H 1
#include <DDalphaAMG.h>

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

#endif /* _MG_INTERFACE_TYPES_H */
