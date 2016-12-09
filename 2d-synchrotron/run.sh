
possibleActions="animate lossmap energy lossmap-analysis"
action="$1"
save_file="$2"

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
	./main $action
	echo " --- Plotting...  "
	python3.5 plot.py $action $save_file
fi
