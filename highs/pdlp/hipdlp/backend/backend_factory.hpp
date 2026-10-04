#ifndef PDLP_HIPDLP_BACKEND_BACKEND_FACTORY_HPP_
#define PDLP_HIPDLP_BACKEND_BACKEND_FACTORY_HPP_

#include <cstddef>
#include <memory>

#include "pdlp/hipdlp/backend/backend.hpp"

namespace hipdlp {
namespace backend {

struct FactoryResult {
  api::Error error;
  std::unique_ptr<Backend> backend;

  FactoryResult() = default;
  FactoryResult(FactoryResult&&) = default;
  FactoryResult& operator=(FactoryResult&&) = default;

  FactoryResult(const FactoryResult&) = delete;
  FactoryResult& operator=(const FactoryResult&) = delete;

  bool ok() const;
};

const Descriptor* registrations(std::size_t& count);
const Descriptor* findRegistration(Kind kind);
FactoryResult createBackend(Kind kind);

}  // namespace backend
}  // namespace hipdlp

#endif
