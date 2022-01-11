cat IDlist2orthofinder | while read i ;
do
echo "Rscript gggenes_orthofinder.R -i $i -o /jelly_data/yichun/myriapod/HGT/seqids/figure/"
Rscript 5gggenes_orthofinder.R -i $i -o /jelly_data/yichun/myriapod/HGT/seqids/figure/
done
