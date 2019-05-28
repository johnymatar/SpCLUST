n=$(nproc --all)
BASEDIR=$(dirname "$0")
if [ "$1" = "" ]
then
    arg1=$BASEDIR/ND3.fasta
	arg2=$BASEDIR/clusters.txt
else
    if [ "$2" = "" ]
	then
		arg1=$1
		arg2=$BASEDIR/clusters.txt
	else
		arg1=$1
		arg2=$2
	fi
fi
mpirun -np $n --allow-run-as-root $BASEDIR/gclust -alignMode fast -mdist PAM250 -in $arg1 -out $arg2
