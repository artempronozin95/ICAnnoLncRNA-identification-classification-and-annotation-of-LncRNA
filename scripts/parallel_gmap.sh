test=$(ls -lrt $2 | cut -f 5 -d ' ')
minimumsize=4000000
if  [ $test -gt $minimumsize ]; then
	mkdir ./data/input/gmap
	awk -v size=10000 -v pre=prefix -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("./data/input/gmap/%s%0"pad"d.fa",pre,n)}}{print>>f}' $2
	parallel -j 40 'python scripts/GMAP.py '$1' {} '$3' '$4' {.}.gff' ::: ./data/input/gmap/*.fa
	sed -i '1,2d' ./data/input/gmap/*.gff
	cat ./data/input/gmap/*.gff > $5
	rm -R ./data/input/gmap
else
	python scripts/GMAP.py $1 $2 $3 $4 $5
fi
