#ifndef PDLP_HIPDLP_BACKEND_COMMAND_BUFFER_HPP_
#define PDLP_HIPDLP_BACKEND_COMMAND_BUFFER_HPP_

#include <cstddef>

#include "pdlp/hipdlp/backend/backend.hpp"

namespace hipdlp {
namespace backend {

struct CommandBufferStatistics {
  std::size_t replay_calls;
  std::size_t direct_operation_enqueue_calls;
  std::size_t allocation_calls;
  std::size_t descriptor_creation_calls;
  std::size_t host_to_device_transfer_calls;
  std::size_t device_to_host_transfer_calls;
  std::size_t synchronization_calls;

  CommandBufferStatistics();
};

class CommandBufferRecorder {
 public:
  virtual ~CommandBufferRecorder() {}

  virtual api::Error beginRecording() = 0;
  virtual api::Error endRecording(CommandBufferHandle& command_buffer) = 0;
  virtual api::Error cancelRecording() = 0;
  virtual api::Error replay(const CommandBufferHandle& command_buffer,
                            CommandBufferStatistics& statistics) = 0;
};

CommandBufferRecorder* commandBufferRecorder(Backend& backend);
void addCommandBufferStatistics(const CommandBufferStatistics& addition,
                                CommandBufferStatistics& total);

}  // namespace backend
}  // namespace hipdlp

#endif