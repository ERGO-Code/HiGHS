#ifndef FACTORHIGHS_RETURN_VALUES_H
#define FACTORHIGHS_RETURN_VALUES_H

namespace hipo {

enum RetValue {
  kRetOk,
  kRetInvalidInput,
  kRetOutOfMemory,
  kRetInvalidPivot,
  kRetMetisError,
  kRetIntOverflow,
  kRetGeneric
};

}

#endif