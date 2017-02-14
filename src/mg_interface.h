#ifndef _MG_INTERFACE_H
#define _MG_INTERFACE_H 1
#include <DDalphaAMG.h>
#include <parser_types.h>
#include <mg_interface_types.h>

mg_state mg_init(struct run_params, qhg_gauge_field);
void mg_invert(qhg_spinor_field, qhg_spinor_field, double, enum mu_sign, mg_state *);
void mg_finalize(void);

#endif /* _MG_INTERFACE_H */
