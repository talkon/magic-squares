#!/bin/bash
# usage: INPUT=[input file] COUNT=[count] LLsub ./submit-sc.sh [LLsub options]
# input file format: each line contains a P tuple

source /etc/profile
module add anaconda/2023a

export PYPY3=python
P_TUPLE=$(head "-$((1 + $LLSUB_RANK))" $INPUT | tail -1)

echo ./run.sh "$P_TUPLE" -c -e -af -pf --count $COUNT
./run.sh "$P_TUPLE" -c -e -af -pf --count $COUNT
