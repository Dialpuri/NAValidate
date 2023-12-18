#!/bin/sh
source "/Applications/ccp4-8.0/bin/ccp4.setup-sh"
export PATH=/usr/bin:$PATH
make -f Makefile.arch -j
./fix_library.sh

./navalidate -pdbin data/1hr2.pdb