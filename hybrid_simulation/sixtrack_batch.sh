#!/bin/bash

# -------------------------------------------------------------------
# user settings
# NB: ALWAYS check:
# - batch queue
# - copy of simulation files (and scripts, manually) to local machine
# - post-processing programs to run
# - which output files should be copied back
# - should previous dirs run* be deleted
# -------------------------------------------------------------------

# - copy simulation files to local machine (leave empty in case you don't want to copy back!):
#   NB: 
#   * if you require to avoid LSF emails, scp-ing won't actually involve the LSF files;
#     thus, each numerical directory will be recreated at the very end of the job,
#     containing only the LSF files;
#   * if you don't scp, the output on terminal of sixtrack is zipped, to save disk space;
#   . FULLPATH (it will be automatically created by the script)
LocalPWD=
#   . (FULL) machine name 
LocalMachine=
#   . user
LocalUser=`whoami`

# - queue type:
#   choose between LSF and HTCONDOR
batchSys=HTCONDOR
#   pick up a value amongst:
#   . LSF: test 8nm 1nh 8nh 1nd 2nd 1nw 2nw (see http://lsf-rrd.cern.ch/lrf-lsf/)
#   . HTCONDOR: espresso (20min) microcentury (1h) longlunch (2h) workday (8h) tomorrow (1d) testmatch (3d) nextweek (1w) (see http://batchdocs.web.cern.ch/batchdocs/local/lsfmigratepractical.html)
queue=tomorrow

# - list of files to be kept:
fileList="DUMP.txt dist0.dat screenout FirstImpacts.dat first_imp_average.dat LPI_BLP_out.s Coll_Scatter_real.dat collgaps.dat coll_summary.dat impacts_real.dat"

# - random seed is randomly chosen (false) or set equal to run index (true)
lseed=false

# - range of numerical dirs:
LIMITLOW=1
LIMITHIGH=2000

# - executables:
SixTOOLS="/afs/cern.ch/user/s/swretbor/collsoft_utilities"
SixExeCollimat="/afs/cern.ch/user/s/swretbor/sixtrack/bin/SixTrackCollimat"
SixExeNonCollimat="/afs/cern.ch/user/s/swretbor/sixtrack/bin/SixTrackNonCollimat"
beamLossPatt="${SixTOOLS}/Post_SixTrack/trunk/BeamLossPattern"
cleanInel="${SixTOOLS}/Post_SixTrack/trunk/CleanInelastic"
cleanCollScatter="${SixTOOLS}/Post_SixTrack/trunk/CleanCollScatter"
cleanCollSum="${SixTOOLS}/Post_SixTrack/trunk/correct_coll_summary.sh"
FirstImpactsAve="${SixTOOLS}/Post_SixTrack/trunk/FirstImpacts--Average.sh"
fort3lib="${SixTOOLS}/Submission_script/fort_3_lib.sh"
HTCONDORtemplate="${SixTOOLS}/Submission_script/htcondor.sub"
extract_chits="/afs/cern.ch/work/s/swretbor/sixtrack/simulations/oml_study/extract_chits.sh"

# Toy model
TMParticleGenerator="/afs/cern.ch/user/s/swretbor/test/cern/synchrotron/exportTool"

# - avoid LSF mails (leave empty in case you want to be spammed):
LSFerrFile=stderr.txt
LSFoutFile=stdout.txt

# -------------------------------------------------------------------
# checks & settings
# -------------------------------------------------------------------

# some sanity checks:
echo "sanity checks about provided infos..."
labort=false
for tmpExe in ${SixExe} ${beamLossPatt} ${cleanInel} ${cleanCollScatter} ${cleanCollSum} ${FirstImpactsAve} ; do
    if [ ! -e ${tmpExe} ] ; then
	echo "fatal: ${tmpExe} does not exists!"
	labort=true
    fi
done
if ${labort} ; then
    echo "aborting..."
    exit 
else
    echo "...all fine!"
fi

# path on lxplus where simulation case is located:
PWD=`pwd`
echo "current dir: $PWD"

# scp? remind the user of the settings:
if [ -n "${LocalPWD}" -a -n "${LocalMachine}" -a -n "${LocalUser}" ] ; then
    echo "warning: copy simulation files to local machine:"
    echo "machine: ${LocalMachine}"
    echo "path: ${LocalPWD}"
    echo "user: ${LocalUser}"
    lcopy=true
else
    echo "warning: no copy of simulation files to local machine required"
    lcopy=false
fi

# batch system
batchSys=`echo "${batchSys}" | awk '{print (toupper($1))}'`
if [ "${batchSys}" != "LSF" ] && [ "${batchSys}" != "HTCONDOR" ] ; then
    echo "unknown batch system ${batchSys}!"
    echo "only LSF and HTCONDOR are presently supported!"
    echo "aborting..."
    exit 1
fi

# LSF emails:
doNotSpamMe=false
if [ "${batchSys}" == "LSF" ] ; then
    if [ -n "${LSFerrFile}" -a -n "${LSFoutFile}" ] ; then
	doNotSpamMe=true
	echo "warning: you won't receive emails from the LSF system;"
	echo "         thus, the same info is contained in ${LSFerrFile} and ${LSFoutFile} files,"
	echo "         locally saved in the numerical directory;"
    fi
fi

# -------------------------------------------------------------------
# actual script
# -------------------------------------------------------------------

# in case, set up directory on local machine
if ${lcopy} ; then
    ssh ${LocalUser}@${LocalMachine} "mkdir -p ${LocalPWD}"
    scp -r $PWD/clean_input ${LocalUser}@${LocalMachine}:"${LocalPWD}"

    # scp executables
#     for tmpExe in ${SixExe} ${beamLossPatt} ${cleanInel} ${cleanCollScatter} ${cleanCollSum} ${FirstImpactsAve} ; do
# 	scp $tmpExe ${LocalUser}@${LocalMachine}:"${LocalPWD}clean_input"
#     done

    # copy this script
    thisscript=`basename $0`
    thisdir=`dirname $0`
    scp ${thisdir}/${thisscript} ${LocalUser}@${LocalMachine}:"${LocalPWD}"/clean_input
fi

# some cleanup (-rf to avoid useless echos...)
if [ "${batchSys}" == "LSF" ] ; then
    rm -rf LSFJOB_*
fi

print_date="echo \$(date)"

# loop over numerical dirs
for ((a = LIMITLOW; a <= LIMITHIGH ; a++)) ; do
    index=`printf "%04i" "$a"`
    if [ -d run$index ] ; then
	echo "...removing exiting run$index !"
	rm -rf run$index
    fi
    mkdir run$index
    echo run$index
    # do not cd run$index, otherwise $PWD changes value

    cat > run$index/SixTr--$index.job << EOF
#!/bin/bash
lcopy=${lcopy}
lseed=${lseed}

# ---------------------------------------------------------------------------------------------------
# load bash library for manipulating fort.3
source ${fort3lib}

# ---------------------------------------------------------------------------------------------------
# copy all needed files in tmp dir (just for the time of running)
cp $PWD/clean_input/*.* .
surveyFile=\`\ls -1 . | \grep -i survey | \head -1\`
collPosFile=\`\ls -1 . | \grep -i collpos | \head -1\`
apeFile=\`\ls -1 . | \grep -i allapert | \head -1\`
if [ "\${surveyFile}" != "SurveyWithCrossing_XP_lowb.dat" ] ; then
   mv \${surveyFile} SurveyWithCrossing_XP_lowb.dat
fi

# ---------------------------------------------------------------------------------------------------
# set seed
if \${lseed} ; then
    # seed is based on index of numerical dir
    set_seed $a
else
    # seed is randomly chosen by SixTrack
    set_seed 0
fi


# ---------------------------------------------------------------------------------------------------
# Part that does the work
# 1. Generate a particle distribution
# 2. Run distribution for n turns with SixTrack
# 3. Convert the output to collimation input
# 4. Run collimation version of SixTrack

echo $print_date
# 1
echo "Generating particles..."
$TMParticleGenerator nocoll . 1 2000 lin+exp ver 0
mv 1.txt fort.13

echo $print_date
# 2
echo "Running non-collimation version..."
mv fort.3_nocoll fort.3
$SixExeNonCollimat > screen_nocoll

echo $print_date
# 3
echo "Converting output..."
python convert_dump_to_input.py ip1_dump.txt coll_input.txt 2000 451209.0195450667524710 # Last number is expected energy after 11 seconds of ramping

echo $print_date
# 4
echo "Running collimation version..."
mv fort.3_coll fort.3
$SixExeCollimat > screen_coll

echo $print_date

# ---------------------------------------------------------------------------------------------------
# post-processing

$beamLossPatt lowb tracks2.dat BLP_out \${apeFile}
# remove binary characters in LPI file:
perl -pi -e 's/\0/ /g' LPI_BLP_out.s 

# clean lists of events in collimators from protons being lost in machine aperture beforehand
$cleanInel FLUKA_impacts.dat LPI_BLP_out.s \${collPosFile} coll_summary.dat
$cleanCollScatter Coll_Scatter.dat LPI_BLP_out.s \${collPosFile} coll_summary.dat

# create clean coll_summary with awk script
$cleanCollSum

# have a first feeling about the average beam impact parameter
${FirstImpactsAve} > first_imp_average.dat

# Extract coll hits from the screenout
#$extract_chits screenout

# keep this for debugging in case of crashes
ls -lh

# copy files back to lxplus dir (-p: preserve date/time)
cp -p ${fileList} $PWD/run$index/
cp -p core* $PWD/run$index/

if \${lcopy} ; then
    gzip $PWD/run${index}/*
    scp -r $PWD/run${index} ${LocalUser}@${LocalMachine}:"${LocalPWD}"
    rm -r $PWD/run${index}/*
else
    # save some disk space
    gzip $PWD/run${index}/screenout
fi

exit
EOF
    if [ -d "run$index" ]; then
	chmod 777 run$index/SixTr--$index.job
	if [ "${batchSys}" == "LSF" ] ; then
	    if ${doNotSpamMe} ; then
		spamString="-o $PWD/run${index}/${LSFoutFile} -e $PWD/run${index}/${LSFerrFile}"
	    else
		spamString=""
	    fi
	    cd run$index
            # -sp : set priority, between 1 and 100
            # -R "rusage[pool=30000]" : allocate 30 GB hard drive space
	    bsub -q $queue -sp 100 -R "rusage[pool=50000]" ${spamString} SixTr--$index.job
            # echo the command line on a txt file, in case you need to submit the job again
	    echo "bsub -q $queue -sp 100 -R \"rusage[pool=50000]\" ${spamString} SixTr--$index.job" > .command.txt
	    cd ../
	elif [ "${batchSys}" == "HTCONDOR" ] ; then
	    echo run$index/SixTr--$index.job >> jobs.txt
	fi
    fi
  
done

if [ "${batchSys}" == "HTCONDOR" ] ; then
    echo " Submitting jobs to htcondor..."
    cp ${HTCONDORtemplate} .
    sed -i "s/^+JobFlavour.*/+JobFlavour = \"${queue}\"/" htcondor.sub
    condor_submit htcondor.sub
    if [ $? -eq 0 ] ; then
	rm -f jobs.txt
    fi
fi
