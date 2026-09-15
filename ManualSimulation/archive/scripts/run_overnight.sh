#!/bin/bash
# run_overnight.sh
# Top-level orchestrator for the whole overnight run. Runs everything as
# a SEQUENCE of separate MATLAB processes (never concurrent ones -- this
# machine has deadlocked before on concurrent MATLAB launches, and has no
# Parallel Computing Toolbox anyway) so it's safe to leave unattended.
#
# Steps:
#   1. Run_GridSearch.m -- main sweep: best pressure per sea state for
#      PassivePump/DHD(2/3/4), plus EHA mech/elec-opt run directly.
#      Produces results/gridSearch/summary.csv and pressure-sweep plots.
#   2. prepareAstarTasks.m -- builds the A*-horizon-length study's task
#      list (DHD2, 3 sea states x 7 horizon lengths = 21 tasks), using
#      step 1's best-pressure results.
#   3. One MATLAB process PER task (runSingleAstarTask.m), each wrapped
#      in `timeout` -- a horizon length that's taking too long is simply
#      killed and skipped, per your instruction. This is a real OS-level
#      kill, which is the only way to enforce that without Parallel
#      Computing Toolbox.
#   4. Aggregate + plot the A*-horizon results.
#
# Usage: ./run_overnight.sh   (run from the ManualSimulation directory,
# or anywhere -- it cd's to its own location first)

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

LOG=overnight_log.txt
ASTAR_TIMEOUT=600  # seconds; per your instruction, skip a horizon length rather than block

echo "=== $(date): Starting overnight run ===" | tee -a "$LOG"

echo "--- Step 1: Run_GridSearch.m (auto-retry on failure/OOM) ---" | tee -a "$LOG"
# This machine is already tight on memory from other running apps (not
# this job) -- an OOM kill mid-sweep is a real possibility, not just a
# theoretical one (already observed once during testing tonight). Each
# task saves its own file immediately and Run_GridSearch.m skips
# already-finished tasks on a fresh run, so simply re-launching after a
# crash resumes safely rather than redoing completed work.
MAX_ATTEMPTS=6
attempt=1
RC=1
while [ $attempt -le $MAX_ATTEMPTS ]; do
    echo "$(date): Run_GridSearch.m attempt $attempt/$MAX_ATTEMPTS" | tee -a "$LOG"
    matlab -nodisplay -nosplash -batch "Run_GridSearch" >> "$LOG" 2>&1
    RC=$?
    if [ $RC -eq 0 ]; then
        break
    fi
    echo "$(date): Run_GridSearch.m exited with code $RC -- retrying in 30s (already-completed tasks are preserved)." | tee -a "$LOG"
    attempt=$((attempt+1))
    sleep 30
done
if [ $RC -ne 0 ]; then
    echo "$(date): Run_GridSearch.m still failing after $MAX_ATTEMPTS attempts -- stopping. Check $LOG." | tee -a "$LOG"
    exit 1
fi

echo "--- Step 2: prepareAstarTasks.m ---" | tee -a "$LOG"
matlab -nodisplay -nosplash -batch "prepareAstarTasks" >> "$LOG" 2>&1

TASKLIST=results/astarHorizon/taskList.csv
if [ ! -f "$TASKLIST" ]; then
    echo "$(date): $TASKLIST not found -- prepareAstarTasks.m must have failed. Skipping A* study." | tee -a "$LOG"
    exit 1
fi

echo "--- Step 3: A* horizon-length tasks (one process each, ${ASTAR_TIMEOUT}s timeout) ---" | tee -a "$LOG"
NTASKS=$(($(wc -l < "$TASKLIST") - 1))  # minus header row
for ((i=1; i<=NTASKS; i++)); do
    echo "$(date): task $i/$NTASKS" | tee -a "$LOG"
    timeout "$ASTAR_TIMEOUT" matlab -nodisplay -nosplash -batch "runSingleAstarTask($i)" >> "$LOG" 2>&1
    if [ $? -eq 124 ]; then
        echo "$(date): task $i TIMED OUT after ${ASTAR_TIMEOUT}s -- skipped." | tee -a "$LOG"
    fi
done

echo "--- Step 4: aggregate + plot A*-horizon results ---" | tee -a "$LOG"
matlab -nodisplay -nosplash -batch "aggregateAstarHorizonResults(); plotAstarHorizonResults();" >> "$LOG" 2>&1

echo "=== $(date): Overnight run complete ===" | tee -a "$LOG"
