#include "pdlp/hipdlp/backend/command_buffer.hpp"

namespace hipdlp {
namespace backend {

CommandBufferStatistics::CommandBufferStatistics()
    : replay_calls(0),
      direct_operation_enqueue_calls(0),
      allocation_calls(0),
      descriptor_creation_calls(0),
      host_to_device_transfer_calls(0),
      device_to_host_transfer_calls(0),
      synchronization_calls(0) {}

CommandBufferRecorder* commandBufferRecorder(Backend& backend) {
  return dynamic_cast<CommandBufferRecorder*>(&backend);
}

void addCommandBufferStatistics(const CommandBufferStatistics& addition,
                                CommandBufferStatistics& total) {
  total.replay_calls += addition.replay_calls;
  total.direct_operation_enqueue_calls +=
      addition.direct_operation_enqueue_calls;
  total.allocation_calls += addition.allocation_calls;
  total.descriptor_creation_calls += addition.descriptor_creation_calls;
  total.host_to_device_transfer_calls += addition.host_to_device_transfer_calls;
  total.device_to_host_transfer_calls += addition.device_to_host_transfer_calls;
  total.synchronization_calls += addition.synchronization_calls;
}

}  // namespace backend
}  // namespace hipdlp