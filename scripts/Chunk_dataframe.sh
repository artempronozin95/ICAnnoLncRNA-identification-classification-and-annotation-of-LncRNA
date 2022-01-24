test=$(ls -lrt $2 | cut -f 5 -d ' ')
minimumsize=4000000000
if  [ $test -gt $minimumsize ]; then
	mkdir ./data/input/chunk
	awk -v size=10000 -v pre=prefix -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("./data/input/chunk/%s%0"pad"d.fa",pre,n)}}{print>>f}' $2
	Rscript scripts/lncFind.r $1 ./data/input/chunk $3 $4
else
	Rscript scripts/lncFind.r $1 $2 $3 $4
fi
