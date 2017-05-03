#!/bin/bash

toymodel="/afs/cern.ch/user/s/swretbor/test/cern/synchrotron"
resources="$toymodel/resources"

# copy them over locally to the batch so we don't overload AFS
cp "$resources/LHC_ramp.dat" .
cp "$resources/motor_tcp_ir3_f5433b1.txt" .
cp "$resources/motor_tcp_ir7_f5433b1.txt" .
cp "$toymodel/2dsynch" .

toymodelBin="$PWD/2dsynch"

jobs=250

LSFerrFile="errfile.txt"
LSFoutFile="outfile.txt"

md=$PWD

for ((j = 1; j <= jobs; j++)) ; do
    rm -rf job$j
    mkdir job$j
    cd job$j
    echo job$j
    cat > 2dsynch_$j.job << EOF
#!/bin/bash

#inputs
cp "$md/LHC_ramp.dat" .
cp "$md/motor_tcp_ir3_f5433b1.txt" .
cp "$md/motor_tcp_ir7_f5433b1.txt" .

$toymodelBin lossmap > stdout.txt

#copy back
cp -p stdout.txt $PWD
cp -p startdist.dat $PWD
cp -p enddist.dat $PWD
cp -p *.ch $PWD # collfiles
exit
EOF
    spamString="-o $PWD/${LSFoutFile} -e $PWD/${LSFerrFile}"
    chmod 777 2dsynch_$j.job
    # -sp : set priority, between 1 and 100
    # -R "rusage[pool=30000]" : allocate 30 GB hard drive space
    bsub -q 2nd -sp 100 -R "rusage[pool=500]" ${spamString} 2dsynch_$j.job
    cd ..
done
