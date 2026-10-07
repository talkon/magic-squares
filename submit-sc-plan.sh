#!/bin/bash
# usage: PLAN=[plan file] LLsub ./submit-sc-plan.sh [LLsub options]
#
# Runs a static plan made by the scheduler, e.g.
#   python3 scripts/scheduler.py emit --units 2000 --unit-time 600 > plan.txt
# Each line of the plan is a list of msearch arguments (one unit of work,
# writing its own output file under data/sched/units/). Task LLSUB_RANK runs
# lines LLSUB_RANK, LLSUB_RANK + LLSUB_SIZE, ... (LLSUB_SIZE defaults to 1
# task per line). Afterwards, copy data/sched/units/ back and run
#   python3 scripts/scheduler.py report    (or fit, plan, emit, run)

source /etc/profile

: ${LLSUB_RANK:=0}
: ${LLSUB_SIZE:=$(wc -l < "$PLAN")}

mkdir -p data/sched/units
line_no=0
while IFS= read -r line; do
  if [ $((line_no % LLSUB_SIZE)) -eq "$LLSUB_RANK" ]; then
    echo bin/msearch $line
    bin/msearch $line
  fi
  line_no=$((line_no + 1))
done < "$PLAN"
