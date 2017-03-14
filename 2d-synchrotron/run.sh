
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
    echo " --- Caching data in '$cache'."

    mkdir $cache
    ./main $action $arg | tee $resultdir/stdout.txt
    cp $resultdir/* $cache
    echo " --- Plotting...  "
    python3.5 plot.py $action $arg
fi
