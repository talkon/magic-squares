#!/bin/bash
# usage: INPUT=[input file] COUNT=[count] LLsub ./submit-sc.sh [LLsub options]
# input file format: each line contains a P tuple

source /etc/profile
module add anaconda/2023a

PYPY3=python
P_TUPLE=$(head "-$LLSUB_RANK" $INPUT | tail -1)

echo ./run.sh "$P_TUPLE" -c --count $COUNT
./run.sh "$P_TUPLE" -c --count $COUNT