cat IDlist | while read i;
do
echo "$i"
ls $i | grep "\.fa" | wc -l
done

#sh count.sh > count.log
#cat count.log | sed ':t;N;s/\n/\t/;b t' | sed 's/\/jelly_data/\n\/jelly_data/g' > fastacount
#awk -F$'\t' '($2 > 1){print}' fastacount | cut -f1 > IDlist2orthofinder