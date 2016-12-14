# Compares the internal Sixtrack total energy with the total 
# momentum of a particle (both given by printouts in sixtrack.s).
#
# Assuming we find 
#	'dynksets.dat' 	: result from using DEBUG in dynk module
#	'six'  			: output from sixtrack containing lines like
#					: 'DDD>ejfv		id		value'
#					  'DDD>ejv 		id 		value' to parse

dr=$1
stdout=$dr/six

echo "Extracting energies from '$stdout'..."

grep -A 1 "DDD>ejfv\s*1\s" $stdout | grep -o '[0-9-]*[.][0-9]*' | paste -d " " - - > part_1_energy.txt

echo "Plotting..."

python3.5 comp_ejv_ejfv.py $dr