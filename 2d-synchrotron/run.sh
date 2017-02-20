
possibleActions="animate animate-background lossmap energy lossmap-analysis startdist phasespace phasespace-mov lost"
action="$1"
arg="$2"

workingdir=`pwd`
resultdir="$workingdir/calc"
moviedir="$workingdir/movies"
cache="$workingdir/simulations/cache/$(date +%Y.%m%.%d.%H.%M.%S)"


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
    echo " --- Running '$action $arg'..."
    echo " --- Compiling c++... "
    make 
    echo " --- Running c++... "
    #if [ "$action" == "phasespace" ]; then
        #for (( i=500; i <= 5000000; i+=11245 ))
        #do
            #echo " --- Turn $i "
            #./main $action $i
            #echo " --- Saving image to '$moviedir/ramp_turn$i.png"
            #python3.5 plot.py $action $moviedir/ramp_turn$i.png
            ##echo " --- Copying result to '$resultdir/lines$i.dat'"
            ##mv $resultdir/lines.dat $resultdir/lines$i.dat
        #done
    #else
        echo " --- Caching data in '$cache'."
        mkdir $cache
        ./main $action $arg | tee $cache/stdout.txt
        cp $resultdir/* $cache
        echo " --- Plotting...  "
        python3.5 plot.py $action $arg
    #fi
fi
