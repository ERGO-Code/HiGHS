#ifndef HIGHS_METIS_SORT_H_
#define HIGHS_METIS_SORT_H_

#include <stddef.h>

#include "metislib.h"
#include "metis_mangle.h"

#ifdef __cplusplus
extern "C" {
#endif

void isortd(size_t n, idx_t* base);
void ikvsorti(size_t n, ikv_t* base);

#ifdef __cplusplus
}
#endif

#endif