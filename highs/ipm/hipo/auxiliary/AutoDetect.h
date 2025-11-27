#ifndef HIPO_AUTO_DETECT_H
#define HIPO_AUTO_DETECT_H

#include <string>

namespace hipo {
// Detect BLAS integer model
enum class IntegerModel { not_set, unknown, lp64, ilp64 };
IntegerModel getBlasIntegerModel();

// Detect Metis integer type
IntegerModel getMetisIntegerModel();

std::string getIntegerModelString(IntegerModel i);

}  // namespace hipo

#endif