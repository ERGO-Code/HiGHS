# Parallel Search Design

General design and summary of draft-implementation of parallel tree search in HiGHS

## HighsMipWorker

This is the only new class! The goal of HighsMipWorker is that nearly everything after presolve should run through it.
Consider it a surrogate for
HighsMipData. When you want to evaluate the root node, run heuristics, run separators, or dive, the code will be
accessing all "global" data from HighsMipWorker. It is therefore now commonly passed as a parameter.

The HighsMipWorker contains the following relevant information:

- HighsDomain (A representation of the global domain that is globally valid but local to the worker)
- HighsCutPool (A pool of globally valid conflicts only accessed by the worker excluding sync points)
- HighsConflictPool (A pool of globally valid cuts only accessed by the worker excluding sync points)
- HighsLpRelaxation (A copy of the LP. Note that this will have cuts both from the global pool and from the worker pool)
- HighsPseudoCost (A copy of the global HighsPseudoCost object with local changes made)
- HighsSearch (A unique pointer)
- HighsSepa (A unique pointer)
- HighsNodeQueue (A local queue that is used when a worker performs backtrack plunge)
- HeurStatistics / SepaStatistics (A way to store statistics and buffer them without changing the global state)
- upper_bound / upper_limit / optimality_limit (Benefit from found solutions without changing the global state)
- Const references to HighsMipSolver and HighsMipSolverData

HighsDomain / HighsCutPool / HighsConflictPool / HighsLpRelaxation / HighsPseudoCost are all currently pointers and
stored in std::deque objects in HighsMipSolverData. This is done for two reasons:

- When starting the parallel search, we need to reassign the master HighsMipWorker to point at new not-the-true-global
  objects. References can't be reassigned, so we'd have to destroy and recreate the worker. On the opposite side, if we
  never reach the tree search, then there's no need to create any new objects.
- Cuts need to have an associated index of where they come from when they're in the LP. Therefore we need to have a
  unique identifier from cut -> cutpool. If the pools are stored in a std:deque, then the index of the pool in the deque
  is that identifier. If the CutPool was only stored to the worker, then there'd need to be a more confusing mapping
  that first goes through the workers.

## General Design Principles

- Parallelism only starts after processing the root node. Before that the master worker has pointers to all the "true"
  global data structures
- No information is shared between workers excluding the dedicated sync points
- There is a central parallel lock, called `parallel_lock` in `HighsMipSolverData`. It is accessed by a function
  `parallelLockActive`, which also consider whether there are additional workers.
- There is a central spawn task function, called `runTask`. It sets the lock, then depending on the amount of tasks and
  user parameters, e.g., `simulate_concurrency`, it either runs them sequentially or in parallel.
- Currently only cuts from a worker's pool that are added to the LP are copied into the global pool.
- Currently, all conflicts from a worker's pool are flushed to the global pool.
- There is a parameter `mip_search_concurrency`, which will create `x` workers per core that HiGHS is using, so if your
  machine has 16 cores, where HiGHS by default will use 8 (half of what's available), and the parameter is set to 2 (
  default currently without any testing), then 16 workers will be spawned, and at most 8 workers will be running at
  once.
- We want more workers than cores (threads that HiGHS will use). That way the chance of waiting on a single task for a
  long time is minimised, and we hope to get more reliable average case performance.
- The general pseudocode (subject to change) of the entire parallel search would be:
    - while (true)
        - Run heuristics (parallel)
        - Sync solutions + sync statistics + prune suboptimal nodes + break if infeasible (serial)
        - Dive (parallel)
        - Sync solutions + sync statistics + prune suboptimal nodes + break if infeasible or some node limit for
          dive round reached (serial)
        - Backtrack plunge (parallel)
        - Break if cant backtrack
    - Push open nodes from workers to global queue (serial)
    - Flush statistics (serial)
    - Sync pools + sync global domain changes found by workers (serial)
    - Propagate global domain with new information (serial)
    - Consider restart
    - While node queue is not empty
        - Sync pseudo costs (serial)
        - Install nodes (serial)
        - Evaluate nodes (parallel)
        - Handle pruned nodes (parallel)
        - If all nodes pruned, sync domain changes + continue
        - Separate (parallel)
        - Stop if infeasible
        - Store basis (parallel)
        - break

## Changes to existing classes

### HighsCutPool

- A std::atomic has been added representing the number of LPs a cut is in. This seems to be a necessary evil (
  alternative would be to create a boolean per worker per cut). Consider the global cut pool, which all workers will
  separate. If worker A and worker B both add the cut into their respective LP, then we need a way to track exactly how
  many LPs each cut is currently in.
- There's also now a buffer called `ageResetWhileLocked_`, which is a flag that mentions that some worker used
  information related to this cut, and it should therefore be having its age reset, but we haven't because doing this
  would produce non-determinism. This has affected the aging code.
- General sync function introduced.

### HighsConflictPool

- Similar to HighsCutPool, there's also now a buffer called `ageResetWhileLocked_`.
- General sync function introduced.

### HighsDomain

- Minor changes to propagation logic to accommodate multiple cut and conflict pools.

### HighsLpRelaxation

- An LpRow now has an index associated to which CutPool the row is from

### HighsPseudoCost

- Has two additional vectors to keep track of which columns have been actively changed. This was done to make the
  syncing phase quicker
- General sync function introduced

### HighsSeparation

- Clique table is now only cleaned up during root node separation
- `separateImpliedBounds` does not produce any new implications during the tree search

### HighsTransformedLp

- Is now passed the global domain (or what the worker believes the global domain currently is)

### HighsPostsolveStack

- Added a thread safe option (copy some data structures) when undoing transformation for a given solution.

### HighsHash

- No clue at all. Some C++ wizardry.

## Expected Questions:

- Why is everything a lambda function `HighsMipSolver`? Answer: Because it was easy to prototype and I didn't have to
  worry about what parameters to pass. They could be changed to standard functions.
- Why have some many files been touched? Answer: Many files were accessing something like `mipsolver.mipdata_->domain`,
  and now `worker.globaldom` has had to be passed through all the functions along that path.
- Is the current code deterministic? Answer: It should be, but there's likely still some non-determinism bugs.
- If I run in serial is the code identical to v1.12? Answer: No. I have tried to make it so, but tracking down these
  small variations is time-consuming and not necessarily beneficial.
- What still needs to be done?
    - Answers:
        - Many hard-coded parameters will need to now consider `num_workers`. The only one I've currently changed is
          `numPlungeNodes`, which had an incredible impact on performance. Which ones and to what values are an open
          question, but this will be necessary for performance. Current example observation: Some problems will
          instantly restart five times in a row because so many nodes are now dumped quickly to the global queue.
        - Extensive determinism and bug checks. I'd like to believe there's not many bugs, but with this much code that
          cannot be.
        - General design review. Examples: (1) Should `HighsMipWorker` be changed (2) Is the pointer to `HighsMipWorker`
          in `HighsLpRelaxation` acceptable (3) Does the way cuts and conflicts are synced make sense w.r.t.
          performance? (4) Are there any potential sync opportunities being missed? (5) Do we have a new performance
          bottleneck?
        - Timers. I have not hacked them in to HighsMipWorker, and have therefore commented many ouy.
        - Merging latest into the branch is going to be a few hours of annoyance. I've held off on doing it so we have
          v1.12 to compare to.
- Is the code in a reviewable state? Answer: I've cleaned up everything, although there's still some lingering WARNINGS
  and TODOS that I've left in case the design changes, e.g., for cut management. So it should be readable and relatively
  clean.