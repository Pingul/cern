# Purpose to compare the Sixtrack results with the toy model. 

#dr=$1
denergy=$1
dr="/Users/swretbor/Workspace/work_afs/sixtrack/simulations/2dsynch_comp/energy-$denergy"
comp_dir=`pwd`

echo "--- Current directory: '$comp_dir'"

# Sixtrack
echo "--- SixTrack: reading from '$dr'"
dynkset=$dr/dynksets.dat
dump=$dr/DUMP.txt
six=$dr/six

# Getting coordinate offsets
echo "--- Parsing '$six'"
grep -A 10 "CLOSED ORBIT AND" $six | grep "/CLS" | grep -o "[0-9-]*[.][0-9]*" | paste -d " " - - > 1p.txt

# Getting coordinates
echo "--- Parsing '$dump'"
grep "^\s*1\s" $dump > __1p.txt

echo "--- Parsing '$dynkset'"
grep "energy" $dynkset > __dynkset.dat
paste __dynkset.dat __1p.txt | awk 'BEGIN {print "#turn Eref Epart phase"} {printf("%d %.17g %.17g %.17g\n", $8, $6, ($6 + $6*$15), $14)}' >> 1p.txt

rm -f __1p.txt __dynkset.dat

# Toy model
tm_dir=../2d-synchrotron
echo "--- ToyModel: path '$tm_dir'"
echo "--- Running simulation"
cd $tm_dir
make
./main sixtrack-comp $denergy
echo "--- Copy back data"
mv calc/toymodel_track.dat $comp_dir
cd $comp_dir

echo "--- Plotting"
python3.5 comp.py $dr
