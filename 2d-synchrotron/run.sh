
possibleActions="animate lossmap energy lossmap-analysis startdist phasespace"
action="$1"
arg="$2"
save_file="$2"

workingdir=`pwd`
resultdir="$workingdir/calc"
moviedir="$workingdir/movies"
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
    if [ "$action" == "phasespace" ]; then
        for (( i=500; i <= 5000000; i+=11245 ))
        do
            echo " --- Turn $i "
            ./main $action $i
            echo " --- Saving image to '$moviedir/ramp_turn$i.png"
            python3.5 plot.py $action $moviedir/ramp_turn$i.png
            #echo " --- Copying result to '$resultdir/lines$i.dat'"
            #mv $resultdir/lines.dat $resultdir/lines$i.dat
        done
    else
        echo " --- Caching data in '$cache'."
        mkdir $cache
        ./main $action $arg | tee $cache/stdout.txt
        cp $resultdir/* $cache
        echo " --- Plotting...  "
        python3.5 plot.py $action $save_file
    fi
fi
