#ifndef _QQ_INTERFACE_H
#define _QQ_INTERFACE_H 1
#include <parser_types.h>
#include <qhg.h>
#include <qq_interface_types.h>

qq_state qq_init(struct run_params, qhg_gauge_field, enum qhg_fermion_bc_time);
void qq_invert(qhg_spinor_field, qhg_spinor_field, double, enum mu_sign, qq_state *);
void qq_finalize(qq_state);

#endif /* _QQ_INTERFACE_H */
