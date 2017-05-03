#!/bin/bash

toymodel="/afs/cern.ch/user/s/swretbor/test/cern/synchrotron"
resources="$toymodel/resources"

declare -a input_files=(
        "LHC_ramp.dat"
        "motor_tcp_ir3_f5433b1.txt"
        "motor_tcp_ir7_f5433b1.txt")
# copy them over locally to the batch so we don't overload AFS
for file in "${input_files[@]}"
do
    cp "$resources/$file" .
done
cp "$toymodel/2dsynch" .

toymodelBin="$PWD/2dsynch"
copy_back="stdout.txt,startdist.dat,enddist.dat,TCP*.dat"

jobs=250

LSFerrFile="errfile.txt"
LSFoutFile="outfile.txt"

for ((j = 1; j <= jobs; j++)) ; do
    rm -rf job$j
    mkdir job$j
    cd job$j
    echo job$j
    cat > 2dsynch_$j.job << EOF
#!/bin/bash

for file in "${input_files[@]}"
do
    cp -p $PWD/$file .
end

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
