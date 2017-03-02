
#!/bin/bash

pushd `dirname $0` > /dev/null
TPATH=$( cd .. ; pwd -P )
popd > /dev/null
pushd `dirname $0` > /dev/null
SCRIPTPATH=$(pwd -P )
popd > /dev/null


RESULT="$TPATH/results_btc/"
mkdir -p $RESULT
rm -rf $RESULT*
conffile="est_2000_5.conf"
datafile="btc_2000_5"



declare -a arr_n=("2000_5")
for nn in ${arr_n[@]}
do
    opt.py $conffile Input.parameterfile -v  ./parameters/param_$nn.dat
    opt.py $conffile Input.datafile -v ./data/$datafile.dat
    opt.py $conffile Domain.nx -v 128
    opt.py $conffile intervalsc1.1.end -v 7980
    opt.py $conffile intervalsc1.2.start -v 7980
    opt.py $conffile Init.maxiterations -v 15
    ./estimation $conffile > $RESULT$nn.log


    gnuplot -e "results='$TPATH/results/'" -e "datapath='$TPATH/data/$datafile'" -e "outputpath = '$RESULT/'"  -e "oname = '$nn'" $SCRIPTPATH/output.gpl

done

cd $RESULT
for f in *.eps;
do
    new_file=$(echo "$f" | sed "s/.eps//g" | sed 's/ *$//g');
    convert -density 300x300 -colorspace RGB -flatten -quality 100% $f $new_file.png;
done
cd $TPATH
