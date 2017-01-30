
possibleActions="animate lossmap energy lossmap-analysis startdist"
action="$1"
save_file="$2"

workingdir=`pwd`
resultdir="$workingdir/calc"
cache="$workingdir/simulations/$(date +%Y.%m%.%d.%H.%M.%S)"


actionOK=0
for a in $possibleActions
do
	if [ "$a" == "$action" ]; then
		actionOK=1;
	fi
done

if [ "$actionOK" == 0 ]; then
	echo " --- $action is not a valid action, choose from "
	echo " --- \t$possibleActions"
else
	echo " --- Running '$action'..."
	if [ ${#save_file} -gt 0 ]; then
		echo " --- Will save potential output to '$save_file'."
	else
		echo " --- Will not save data."
	fi

	echo " --- Compiling c++... "
	make 
	echo " --- Running c++... "
	echo " --- Caching data in '$cache'."
	mkdir $cache
	./main $action | tee $cache/stdout.txt
	cp $resultdir/* $cache
	echo " --- Plotting...  "
	python3.5 plot.py $action $save_file
fi
