
possibleActions="animate animate-background lossmap energy lossmap-analysis startdist phasespace phasespace-mov lost"
action="$1"
arg="$2"

workingdir=`pwd`
resultdir="$workingdir/calc"
moviedir="$workingdir/movies"
cache="$workingdir/simulations/cache/$action-$(date +%Y.%m%.%d.%H.%M.%S)"
stdout=$resultdir/stdout.txt

echo "2D-SYNCHROTRON $(date +%Y.%m%.%d.%H.%M.%S)\n\tAction: $action\n\tArg: $arg\n***************" > $stdout

actionOK=0
for a in $possibleActions
do
    if [ "$a" == "$action" ]; then
        actionOK=1;
    fi
done

if [ "$actionOK" == 0 ]; then
    echo " --- $action is not a valid action, choose from " | tee -a $stdout
    echo " --- \t$possibleActions" | tee -a $stdout
else
    echo " --- Running '$action $arg'..." | tee -a $stdout
    echo " --- Compiling c++... " | tee -a $stdout
    make 
    echo " --- Running c++... " | tee -a $stdout
    echo " --- Caching data in '$cache'." | tee -a $stdout

    mkdir $cache
    ./main $action $arg | tee -a $stdout
    cp $resultdir/* $cache
    echo " --- Plotting...  " | tee -a $stdout
    python3.5 -u plot.py $action $arg | tee -a $stdout
fi
