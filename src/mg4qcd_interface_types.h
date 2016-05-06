#ifndef _MG4QCD_INTERFACE_TYPES_H
#define _MG4QCD_INTERFACE_TYPES_H 1
#include <mg4qcd.h>

enum mu_sign {
  plus,
  minus
};

typedef struct {
  MG4QCD_Init mg_init;
  MG4QCD_Parameters mg_params;
  MG4QCD_Status mg_status;
  enum mu_sign current_mu_sign;
} mg4qcd_state;

#endif /* _MG4QCD_INTERFACE_TYPES_H */
