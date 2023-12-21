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





/home/u0jana01/anaconda3/envs/IGv2/bin/sqStoreCreate \
  -o ./canu.seqStore.BUILDING \
  -minlength 1000 \
  -genomesize 3000 \
  -coverage   50 \
  -bias       0 \
  -raw -trimmed -pacbio-hifi reads /home/u0jana01/IGenotyper/test/test/tmp/assembly/chr15/35927334_35928334/0/reads.fasta \
&& \
mv ./canu.seqStore.BUILDING ./canu.seqStore \
&& \
exit 0

exit 1
