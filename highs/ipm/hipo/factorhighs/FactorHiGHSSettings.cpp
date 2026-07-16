#include "ipm/hipo/factorhighs/FactorHiGHSSettings.h"

namespace hipo {

HipoTuning& hipoTuning() {
  static HipoTuning tuning{};
  return tuning;
}

}  // namespace hipo
