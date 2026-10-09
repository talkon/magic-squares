#!/bin/sh
# The gate of the class support in the V_d searches (msearch
# --class-support-min-labels L, default 137: only the V_d with at least L
# labels; dfirst_set_class_min_labels).
# usage: class_gate_test.sh MSEARCH   (ctest fast_msearch_class_gate)
m="$1 --diag-first --diag-first-min-n 0 --dfirst-star 0"
fail() {
  echo "FAIL: $*"
  exit 1
}
field() { # field NAME: the value of NAME in the dsum record on stdin
  grep '"type":"dsum"' | grep -o "\"$1\":[0-9]*" | cut -d: -f2
}
pairs() { # the (d, hash) of the dsquare records on stdin, sorted
  grep '"type":"dsquare"' | grep -o '"d":[0-9]*\|"hash":"[0-9a-f]*"' | paste -d' ' - - | sort
}

# 16 5 4 2 / 849: 102 labels, every V_d below the gate: no class support,
# the nodes of --no-class-support; with L = 0 every V_d has it, fewer
# nodes, the same pair
a="849 16 5 4 2"
o=$($m --no-class-support --sums $a)
g=$($m --sums $a)
z=$($m --class-support-min-labels 0 --sums $a)
[ "$(echo "$g" | field nd_class)" = 0 ] || fail "849: nd_class $(echo "$g" | field nd_class)"
[ "$(echo "$g" | field nodes)" = "$(echo "$o" | field nodes)" ] || fail "849: gated nodes"
[ "$(echo "$z" | field nd_class)" = "$(echo "$z" | field nd)" ] || fail "849: L 0 nd_class"
[ "$(echo "$z" | field nodes)" -lt "$(echo "$o" | field nodes)" ] || fail "849: L 0 nodes"
[ "$(echo "$o" | pairs)" = "$(echo "$g" | pairs)" ] && [ "$(echo "$o" | pairs)" = "$(echo "$z" | pairs)" ] ||
  fail "849: pairs"
[ -n "$(echo "$o" | pairs)" ] || fail "849: no pair"
echo "$o" | grep '"type":"dsum"' | grep -q nd_class && fail "849: nd_class without class support"

# 13 7 4 3 1 1 / 1900, d 0..100: 137 labels, V_d of 129-137 labels (d 23
# and 91 of 137): the gate splits them; the same pairs in all three runs,
# the gated nodes between the two others
a="1900 13 7 4 3 1 1"
o=$($m --no-class-support --d-range 0:100 --sums $a)
g=$($m --d-range 0:100 --sums $a)
z=$($m --class-support-min-labels 0 --d-range 0:100 --sums $a)
nd=$(echo "$g" | field nd)
nc=$(echo "$g" | field nd_class)
[ "$nd" = 100 ] && [ "$nc" -gt 0 ] && [ "$nc" -lt "$nd" ] || fail "1900: nd $nd nd_class $nc"
[ "$(echo "$z" | field nd_class)" = 100 ] || fail "1900: L 0 nd_class"
no=$(echo "$o" | field nodes)
ng=$(echo "$g" | field nodes)
nz=$(echo "$z" | field nodes)
[ "$nz" -lt "$ng" ] && [ "$ng" -lt "$no" ] || fail "1900: nodes $no $ng $nz"
[ "$(echo "$o" | pairs)" = "$(echo "$g" | pairs)" ] && [ "$(echo "$o" | pairs)" = "$(echo "$z" | pairs)" ] ||
  fail "1900: pairs"
# a gate above every V_d: as without the class support
h=$($m --class-support-min-labels 100000 --d-range 0:100 --sums $a)
[ "$(echo "$h" | field nodes)" = "$no" ] && [ "$(echo "$h" | field nd_class)" = 0 ] || fail "1900: L max"
# the plain root (no top-label root): the class support cannot run, so no
# d counts in nd_class, and the nodes are those without it
r=$($m --class-support-min-labels 0 --d-plain-root --d-range 0:100 --sums $a)
rn=$($m --no-class-support --d-plain-root --d-range 0:100 --sums $a)
[ "$(echo "$r" | field nd_class)" = 0 ] && [ "$(echo "$r" | field nodes)" = "$(echo "$rn" | field nodes)" ] ||
  fail "1900: plain root nd_class $(echo "$r" | field nd_class)"
echo "class support gate ok"
