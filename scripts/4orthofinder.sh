cat IDlist2orthofinder | while read i ;
do
echo "$i"
orthofinder -f $i -t 38 -a 38
cp $i/OrthoFinder/Results_Oct05/Orthogroups/Orthogroups.txt $i/Orthogroups.txt
done
