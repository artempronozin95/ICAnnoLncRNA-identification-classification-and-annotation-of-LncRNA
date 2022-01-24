test=$(ls -lrt $1 | cut -f 5 -d ' ')
minimumsize=4000000000
if  [ $test -gt $minimumsize ]; then
	mkdir ./data/input/tmhmm
	awk -v size=10000 -v pre=prefix -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("./data/input/tmhmm/%s%0"pad"d.fa",pre,n)}}{print>>f}' $1
	parallel -j 40 './tmhmm/bin/tmhmm {} > {.}.csv' ::: ./data/input/tmhmm/*.fa
	cat ./data/input/tmhmm/*.csv > $2/tmhmm.csv
#	rm -R ./data/input/tmhmm
else
	./tmhmm/bin/tmhmm $1 > $2/tmhmm.csv
fi
