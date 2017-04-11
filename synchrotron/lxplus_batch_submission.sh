#!/bin/bash

toymodel="/afs/cern.ch/user/s/swretbor/test/cern/synchrotron"
resources="$toymodel/resources"

declare -a precopy_for_batch=(
    "$resources/LHC_ramp.dat"
    "$resources/motor_tcp.txt"
    "$toymodel/2dsynch")
# copy them over locally to the batch so we don't overload AFS
for file in "${precopy_for_batch[@]}"
do
    cp "$file" .
done

toymodelBin="$PWD/2dsynch"
# No spaces inbetween file names
input_files="$PWD/LHC_ramp.dat,$PWD/motor_tcp.txt"
copy_back="stdout.txt,startdist.dat,enddist.dat,coll.dat"

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
    cd ..
done
