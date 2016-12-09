dr=$1
stdout=$dr/six

echo "Extracting energies from '$stdout'..."

grep -A 1 "DDD>ejfv\s*1\s" $stdout | grep -o '[0-9-]*[.][0-9]*' | paste -d " " - - > part_1_energy.txt

echo "Plotting..."

python3.5 comp_ejv_ejfv.py $dr