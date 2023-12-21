#!/bin/sh


#  Path to Canu.

syst=`uname -s`
arch=`uname -m | sed s/x86_64/amd64/`

bin="/home/u0jana01/anaconda3/envs/IGv2/bin/$syst-$arch/bin"

if [ ! -d "$bin" ] ; then
  bin="/home/u0jana01/anaconda3/envs/IGv2/bin"
fi

#  Report paths.

echo ""
echo "Found perl:"
echo "  " `which perl`
echo "  " `perl --version | grep version`
echo ""
echo "Found java:"
echo "  " `which /home/u0jana01/anaconda3/envs/IGv2/bin/java`
echo "  " `/home/u0jana01/anaconda3/envs/IGv2/bin/java -showversion 2>&1 | head -n 1`
echo ""
echo "Found canu:"
echo "  " $bin/canu
echo "  " `$bin/canu -version`
echo ""


#  Environment for any object storage.

export CANU_OBJECT_STORE_CLIENT=
export CANU_OBJECT_STORE_CLIENT_UA=
export CANU_OBJECT_STORE_CLIENT_DA=
export CANU_OBJECT_STORE_NAMESPACE=
export CANU_OBJECT_STORE_PROJECT=




#  Discover the job ID to run, from either a grid environment variable and a
#  command line offset, or directly from the command line.
#
if [ x$CANU_LOCAL_JOB_ID = x -o x$CANU_LOCAL_JOB_ID = xundefined -o x$CANU_LOCAL_JOB_ID = x0 ]; then
  baseid=$1
  offset=0
else
  baseid=$CANU_LOCAL_JOB_ID
  offset=$1
fi
if [ x$offset = x ]; then
  offset=0
fi
if [ x$baseid = x ]; then
  echo Error: I need CANU_LOCAL_JOB_ID set, or a job index on the command line.
  exit
fi
jobid=`expr -- $baseid + $offset`
if [ x$CANU_LOCAL_JOB_ID = x ]; then
  echo Running job $jobid based on command line options.
else
  echo Running job $jobid based on CANU_LOCAL_JOB_ID=$CANU_LOCAL_JOB_ID and offset=$offset.
fi

if [ $jobid -eq 1 ] ; then
  bat="001"
  job="001/000001"
  opt="-h 1-34 -r 1-34 --hashdatalen 168635"
fi


if [ ! -d ./$bat ]; then
  mkdir ./$bat
fi


if [ -e $job.ovb ]; then
  exists=true
else
  exists=false
fi
if [ $exists = true ] ; then
  echo Job previously completed successfully.
  exit
fi

#  Fetch the frequent kmers, if needed.
if [ ! -e ../0-mercounts/canu.ms22.dump ] ; then
  mkdir -p ../0-mercounts
  cd ../0-mercounts
  cd -
fi


$bin/overlapInCore \
  -t 8 \
  -k 22 \
  -k ../0-mercounts/canu.ms22.dump \
  --hashbits 22 \
  --hashload 0.8 \
  --maxerate  0.01 \
  --minlength 500 \
  $opt \
  -o ./$job.ovb.WORKING \
  -s ./$job.stats \
  ../../canu.seqStore \
&& \
mv ./$job.ovb.WORKING ./$job.ovb


exit 0
