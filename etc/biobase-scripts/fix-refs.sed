#!/bin/ksh

input=${1:?Did not specify input file}
print "Fixing references in place in file ${input}"
sed -i '/^BF / s/ Species:.*Species:/ Species:/' ${input}
