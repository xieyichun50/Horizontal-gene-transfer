# Run scripts 
Run in following order and be aware of working directory and file location.

#### Formating the Alienness result table
```
Rscript 1extract_up_down_stream_genes.R 
```
#### Extracting upstream and downstream geneID (by default, +/- 5 genes of the HGT)
```
Rscript 2extract_up_down_stream_genes_jelly.R 
```
#### Extract protein sequence, one HGT event site per folder, one species one fasta file
```
sh 3getfasta.sh
```
#### Exclude HGTs only found in one species
```
sh 4count.sh > count.log
cat count.log | sed ':t;N;s/\n/\t/;b t' | sed 's/\/jelly_data/\n\/jelly_data/g' > fastacount
awk -F$'\t' '($2 > 1){print}' fastacount | cut -f1 > IDlist2orthofinder
```
#### Indentify orthologs on HGT neighbour genes
```
sh 4orthofinder.sh
```
#### Visualise the gene order and orthologues relationship
```
sh 5gggenes_orthofinder.sh #(batch run 5gggenes_orthofinder.R)
```
