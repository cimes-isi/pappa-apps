# simple driver for tce.py
# Copyright (C) 2013, 2014 ET International, Inc.
# written by Mark Glines

import sys
import tce

if(len(sys.argv) < 3):
    print 'Usage: '+sys.argv[0]+' <functionname> <inputfile>'
    sys.exit(1)

basename = sys.argv[1]
fn = sys.argv[2]

print 'Reading file:', fn
ltc = tce.readfromfile(fn)
# ltc is an object of type tce.ListTensorContractions
# uncomment this to dump latex for the input math:
#tce.writetofile(ltc.tex(),basename+'.in.tex')
ot = ltc.breakdown()
# ot is an object of type tce.OperationTree
ot.fullyfactorize()
# dump latex for the output (force-reduction-optimized) math
tce.writetofile(ot.tex(),basename+'.tex')
#ot.fortran77(basename).writetofile(basename)
#ot.fortran90(basename)
#ot.pythongen(basename)
ot.isoc99(basename,[],[],[],[],0).writetofile(basename)
