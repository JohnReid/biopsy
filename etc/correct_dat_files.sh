#!/bin/sh

# Fix species lines in matrix and site tables
#sed -r 's/(BF  .*Species:.*) Species:.*/\1./g' ~/Data/serialised/BF.txt
for dat in transfac/matrix.dat transfac/site.dat
do
  sed -i.bak -r 's/(BF  .*Species:.*) Species:.*/\1./g' $dat
done

# remove spare semi-colon in transcompel
sed -i.bak -r 's/^DR  EMBL:;/DR  EMBL:/' transcompel/compel.dat
