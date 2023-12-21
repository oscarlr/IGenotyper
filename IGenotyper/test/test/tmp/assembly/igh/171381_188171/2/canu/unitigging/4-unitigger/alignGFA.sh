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








if [ ! -e ./canu.unitigs.aligned.gfa ] ; then

  $bin/alignGFA \
    -T ../canu.utgStore 2 \
    -i ./canu.unitigs.gfa \
    -e 0.05 \
    -o ./canu.unitigs.aligned.gfa \
    -t 4 \
  > ./canu.unitigs.aligned.gfa.err 2>&1
fi


if [ ! -e ./canu.unitigs.aligned.bed ] ; then

  $bin/alignGFA -bed \
    -T ../canu.utgStore 2 \
    -C ../canu.ctgStore 2 \
    -i ./canu.unitigs.bed \
    -o ./canu.unitigs.aligned.bed \
    -t 4 \
  > ./canu.unitigs.aligned.bed.err 2>&1
fi


if [ -e ./canu.unitigs.aligned.gfa -a \
     -e ./canu.unitigs.aligned.bed ] ; then
  echo GFA alignments updated.
  exit 0
else
  echo GFA alignments failed.
  exit 1
fi
