# Purpose to compare the Sixtrack results with the toy model. 
# Needs a Sixtrack simulation with DEBUG in the DYNK module, and with DUMP sampling every turn.

dr=$1
comp_dir=`pwd`

echo "Curent directory: '$comp_dir'"

# Sixtrack
echo "Fetching data from Sixtrack..."
echo "	reading from '$dr'"
dynkset=$dr/dynksets.dat
dump=$dr/DUMP.txt

cat $dump | grep "^\s*1\s" > __1p.txt
tail -n +2 $dynkset > __dynkset.dat
paste __dynkset.dat __1p.txt | awk 'BEGIN {print "#turn Eref Epart phase"} {printf("%d %f %f %f\n", $8, $6, ($6 + $6*$15), $14)}' > 1p.txt

rm -f __1p.txt __dynkset.dat

# Toy model
tm_dir=../2d-synchrotron
echo "Running toy model simulation..."
echo "	path '$tm_dir'"

echo "--- START"
cd $tm_dir
make
./main sixtrack-comp
mv calc/toymodel_track.dat $comp_dir
cd $comp_dir
echo "--- FINISH"

echo "Plotting..."
python3.5 comp.py $dr