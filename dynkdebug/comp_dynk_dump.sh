dr=$1

# Meant to compare if the real particle energy (given by printouts in
# sixtrack.s) can be reconstructed using the DUMP + ramp

# Can only handle 2 particle simulations

# Assuming we find 
#	'dynksets.dat' 	: result from using DEBUG in dynk module
#	'six'  			: output from sixtrack containing lines like
#					: 'DDD>		id		energy' to parse
#	'DUMP.txt'		: normal DUMP file. Sampling once every lap.

# These needs to exists for the python script
# dynkset=$dr/dynksets.dat
# dump=$dr/DUMP.txt

stdout=$dr/six

echo "Extracting energies from '$stdout'..."

cat $stdout | grep "DDD>\s*1\s" | sed 's/DDD>[ \t]*//g' > "p1.txt"
cat $stdout | grep "DDD>\s*2\s" | sed 's/DDD>[ \t]*//g' > "p2.txt"

echo "Plotting..."

python3.5 comp_dynk_dump.py $dr