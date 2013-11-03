
#ifndef BIO_GSL_H_
#define BIO_GSL_H_

#include "bio/defs.h"
#include "bio/raii.h"

#include <gsl/gsl_integration.h>

BIO_NS_START

void
gsl_init();

template <>
RAII<gsl_integration_workspace *>::~RAII();

#define BIO_GSL_EXP(x) ((x) < GSL_LOG_DBL_MIN ? 0.0 : gsl_sf_exp(x))

BIO_NS_END


#endif //BIO_GSL_H_

