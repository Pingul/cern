#!/bin/bash

LocalPWD=`pwd`
LocalMachine="lxplus.cern.ch"
LocalUser="swretbor"

toymodel="/afs/cern.ch/user/s/swretbor/test/cern/2d-synchrotron"
resources="$toymodel/resources"
toymodelBin="$toymodel/main"

# No spaces inbetween file names
input_files="$resources/LHC_ramp.dat,$resources/motor_tcp.txt"
copy_back="stdout.txt,startdist.dat,coll.dat"

jobs=500

LSFerrFile="errfile.txt"
LSFoutFile="outfile.txt"

for ((j = 1; j <= jobs; j++)) ; do
    rm -rf job$j
    mkdir job$j
    cd job$j
    echo job$j
    cat > 2dsynch_$j.job << EOF
#!/bin/bash

cp -p {$input_files} .

$toymodelBin lossmap > stdout.txt

cp -p {$copy_back} $PWD
exit
EOF
    spamString="-o $PWD/${LSFoutFile} -e $PWD/${LSFerrFile}"
    chmod 777 2dsynch_$j.job
    # -sp : set priority, between 1 and 100
    # -R "rusage[pool=30000]" : allocate 30 GB hard drive space
    bsub -q 2nd -sp 100 -R "rusage[pool=500]" ${spamString} 2dsynch_$j.job
    sleep 0.3
    cd ..
done
