# Security Policy

## Reporting a vulnerability

Please **do not** report security vulnerabilities through public GitHub issues,
discussions or pull requests.

Instead, report them privately via GitHub's
[private vulnerability reporting](https://github.com/ERGO-Code/HiGHS/security/advisories/new)
(the **Security** tab → **Report a vulnerability**).

Please include:

- the affected HiGHS version, interface (C++, C, Python/`highspy`, etc.) and platform
- a description of the issue and its impact
- the attack vector, including what an attacker needs to control (see *Scope* below)
- a minimal reproducer (e.g. a model file or short script), if possible

The HiGHS maintainers are a small team whose expertise is in optimization rather
than security. We will acknowledge reports as soon as we can, work with you on a
fix, and publish a GitHub Security Advisory once a fix is released.

## Supported versions

Security fixes are made to the latest release only.

## Scope

HiGHS is a library for solving optimization problems. It is not designed to be
exposed directly to untrusted users: anyone who can freely call HiGHS functions
must be considered to have the same privileges as the calling process.

The following are generally **bugs, not security issues** (please report them as
normal issues):

- long run times or high memory use on difficult models (many problems HiGHS
  solves are NP-hard, so this is expected behaviour)
- crashes, hangs or wrong results caused by misuse of the API or invalid option
  values from a caller who already controls the process

The following **are** in scope:

- memory-safety issues (e.g. buffer overflows, out-of-bounds reads/writes)
  triggered by reading a crafted model file (MPS, LP, etc.)
- compromise of official release artifacts or the build/release pipeline
  (e.g. PyPI wheels, GitHub release binaries)

If you solve models from untrusted sources, run HiGHS in an isolated process and
set resource limits such as the `time_limit` option.
