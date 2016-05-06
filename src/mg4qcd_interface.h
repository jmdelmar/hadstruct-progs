#ifndef _MG4QCD_INTERFACE_H
#define _MG4QCD_INTERFACE_H 1
#include <mg4qcd.h>
#include <parser_types.h>
#include <mg4qcd_interface_types.h>

mg4qcd_state mg4qcd_init(struct run_params, qhg_gauge_field);
void mg4qcd_invert(qhg_spinor_field, qhg_spinor_field, double, enum mu_sign, mg4qcd_state *);
void mg4qcd_finalize(void);

#endif /* _MG4QCD_INTERFACE_H */
