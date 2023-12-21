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

if [ -e ../canu.ctgStore/seqDB.v001.tig -a -e ../canu.utgStore/seqDB.v001.tig ] ; then
  exit 0
fi

#
#  Check if the outputs exist.
#
#  The boilerplate function for doing this fails if the file isn't
#  strictly below the current directory, so some gymnastics is needed.
#

cd ..

if [ -e canu.ctgStore/seqDB.v001.tig ]; then
  exists=true
else
  exists=false
fi

if [ -e canu.utgStore/seqDB.v001.tig ]; then
  exists=true
else
  exists=false
fi

cd 4-unitigger

#
#  Run if needed.
#

if [ $exists = false ] ; then
  $bin/bogart \
    -S ../../canu.seqStore \
    -O    ../canu.ovlStore \
    -o     ./canu \
    -gs 18790 \
    -eg 0.01 \
    -eM 0.01 \
    -mo 500 \
    -covgapolap 500 \
    -lopsided nobest 50  \
    -minolappercent   0.0  \
    -dg 12 \
    -db 12 \
    -dr 6 \
    -ca 2100 \
    -cp 200 \
    -threads 4 \
    -M 16 \
    -unassembled 2 0 1.0 0.5 3 \
    -eg 0.0003 -sb 0.01 -dg 0 -db 3 -dr 0 -ca 50 -cp 5 \
    > ./unitigger.err 2>&1 \
  && \
  mv ./canu.ctgStore ../canu.ctgStore \
  && \
  mv ./canu.utgStore ../canu.utgStore
fi

if [ ! -e ../canu.ctgStore -o \
     ! -e ../canu.utgStore ] ; then
  echo bogart appears to have failed.  No canu.ctgStore or canu.utgStore found.
  exit 1
fi

if [ ! -e ../canu.ctgStore/seqDB.v001.sizes.txt ] ; then
  $bin/tgStoreDump \
    -S ../../canu.seqStore \
    -T ../canu.ctgStore 1 \
    -sizes -s 18790 \
   > ../canu.ctgStore/seqDB.v001.sizes.txt
fi



cd ../canu.ctgStore
cd -

cd ../canu.utgStore
cd -


exit 0