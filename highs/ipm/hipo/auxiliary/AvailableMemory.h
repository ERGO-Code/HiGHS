#ifndef HIPO_AVAILABLE_MEMORY_H
#define HIPO_AVAILABLE_MEMORY_H

namespace hipo {

// Memory that can be used without swapping, in bytes, or negative if unknown.
// On Linux this respects the cgroup limits of a container or systemd unit. A
// process limit (ulimit) needs no check, since allocation then fails
double availableMemoryBytes();

}  // namespace hipo

#endif
