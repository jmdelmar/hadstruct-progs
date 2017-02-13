#ifndef _QQ_INTERFACE_TYPES_H
#define _QQ_INTERFACE_TYPES_H 1
#include <quda.h>

enum mu_sign {
  plus,
  minus
};

typedef struct {
  QudaGaugeParam gauge_param;
  QudaInvertParam invert_param;
  QudaMultigridParam mg_param;
} qq_state;

#endif /* _QQ_INTERFACE_TYPES_H */
