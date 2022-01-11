find /jelly_data/yichun/myriapod/HGT/seqids/ | xargs ls -lrt | grep "1:$" | sed 's/\://' > IDlist
cat IDlist | while read i;
do
echo "$i"
ls $i

ls $i | grep -v "\." | while read j;
do
grep -v "HGT" $i/$j | cut -f5 | sort | uniq > $i/$j.id
seqtk subseq /jelly_data/yichun/myriapod/longest_protH/$j.fa $i/$j.id > $i/$j.fa
find $i -name "*" -type f -size 0c | xargs -n 1 rm
done

done
