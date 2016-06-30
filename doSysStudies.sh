for var in `seq 1 10`; do
    for cut in `seq 1 3`; do
	echo Topol $var $cut;
	./runMe Topol $var $cut;
    done
done

