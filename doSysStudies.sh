# run topological cut analysis

for var in `seq 1 10`; do
    for cut in `seq 1 3`; do
	echo Topol $var $cut;
	./runMe Topol $var $cut;
    done
done

# run avg sep cut analysis
for var in `seq 0 5`; do
    for cut in `seq 1 3`; do
	echo AvgSep $var $cut;
	./runMe AvgSep $var $cut;
    done
done
