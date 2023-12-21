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

if [ $jobid -gt 4 ]; then
  echo Error: Only 4 partitions, you asked for $jobid.
  exit 1
fi

if [ $jobid -le 3 ] ; then
  tag="ctg"
else
  tag="utg"
  jobid=`expr $jobid - 3`
fi

jobid=`printf %04d $jobid`

if [ ! -d ./${tag}cns ] ; then
  mkdir -p ./${tag}cns
fi

if [ -e ./${tag}cns/$jobid.cns ] ; then
  exit 0
fi



$bin/utgcns \
  -R ../canu.${tag}Store/partition.$jobid \
  -T ../canu.${tag}Store 1 \
  -P $jobid \
  -O ./${tag}cns/$jobid.cns.WORKING \
  -maxcoverage 40 \
  -e 0.05 \
  -pbdagcon \
  -edlib    \
  -threads 4 \
&& \
mv ./${tag}cns/$jobid.cns.WORKING ./${tag}cns/$jobid.cns \


exit 0
