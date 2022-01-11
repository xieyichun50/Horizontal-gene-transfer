# Horizontal gene transfer visualisation
Visualising horizontal gene transfer results from Alienness [Rancurel et al. 2017](https://www.mdpi.com/2073-4425/8/10/248).
### Require the following files:

1. Alienness summary table : txt
```
Ato_Hho_Nno_Tco WP_052600908.1 / none
Ato_Hho_Nno_Tco AQQ75061.1 / none
Ato_Hho_Nno     WP_095524423.1 / none
Ato_Hho_Nno     WP_163310079.1 / none
Ato_Hho_Nno     WP_052955415.1 / none
Ato_Hho_Nno     RPH48231.1 / none
Ato_Hho_Nno     WP_034862301.1 / none
Ato_Hho_Nno     WP_129730820.1 / none
Ato_Hho_Nno     WP_143274636.1 / none
Ato_Hho_Nno     OYY43986.1 / none
```
2. Alienness species result table : \*_stat_queries_1_likely_hgt.xls
```
ai      hgtscore        query   best_donor_acc / best_toi_acc   best_donor_pident       best_donor_orgname      best_donor_taxonomy     nb_hits_supporting_taxo nb_hits_between_donor_and_possible_toi nb_total_hits   nb_unknown_acc (dr)     nb_excluded_acc (dr)
460.52  881.70  Hho_000036-T1   WP_158239238.1 / none   41.6    Uliginosibacterium sp. TH139    Bacteria        25      25      25      0       0
460.52  953.70  Hho_000093-T1   WP_079420767.1 / none   67.3    Thiomonas intermedia    Bacteria        22      25      25      3       0
460.52  2400.90 Hho_000154-T1   AMP10133.1 / none       68.8    Collimonas arenae       Bacteria        22      25      25      3       0
460.52  727.20  Hho_000425-T1   WP_185821637.1 / none   69.7    Xanthomonas sp. GW      Bacteria        25      25      25      0       0
460.52  836.60  Hho_000432-T1   WP_095524367.1 / none   66.5    Candidatus Dactylopiibacterium carminicum       Bacteria        24      25      25      1       0
460.52  837.80  Hho_000433-T1   WP_170020812.1 / none   38.8    Uliginosibacterium sp. IMCC34675        Bacteria        25      25      25      0       0
460.52  807.40  Hho_000462-T1   WP_170021018.1 / none   68.4    Uliginosibacterium sp. IMCC34675        Bacteria        23      25      25      2       0
460.52  1805.00 Hho_000465-T1   OYY43986.1 / none       48.1    Gallionellales bacterium 35-53-114      Bacteria        21      25      25      4       0
```
3. Gene annotation file : gff3, longest protein filtered, for determination of up-/down-stream genes
4. Protein sequence file : fasta, longest protein filtered

### Require the following environment/software
1. #!/bin/bash : to run \*.sh 
2. [orthofinder](https://github.com/davidemms/OrthoFinder)
3. [seqtk](https://github.com/lh3/seqtk)
4. Rscript : to run \*.R  
5. In R, [gggenes](https://github.com/wilkox/gggenes)
6. R packages: 'tidyr', 'dplyr', 'stringr', 'optparse', 'ggplot2', 'ggforce', 'RColorBrewer'

