# IgnoroMeNot
IgnoroMeNot outputs a list of ignorome genes highly associated with other well-annotated genes. Ignorome genes are genes have little to no validated experimental Gene Ontology annotations. Strong associations between ignorome genes and well-annotated genes can help fulfill the protein function space and minimize the gap in knowledge between "rich" and "poor" genes.

### Required files:
IgnoroMeNot requires a tab-separated file of a list of genes with any metric (GO annotation count, information content, article count). Most and least annotated genes are defined based on this metric.

Additionally, a STRING interaction network and a protein alias file of the organism to be examined can be downloaded from string-db.org/cgi/download. Example: E. coli alias file __https://stringdb-static.org/download/protein.aliases.v11.5/511145.protein.aliases.v11.5.txt.gz__ and interaction network __https://stringdb-static.org/download/protein.links.full.v11.5/511145.protein.links.full.v11.5.txt.gz__

## IgnoroMeNot usage:
`` usage: ignoromenot.py [-h] --ranktable RANKTABLE --idtable IDTABLE --stringppi STRINGPPI
                      [--threshold_top THRESHOLD_TOP | --percentile_top PERCENTILE_TOP]
                      [--threshold_bot THRESHOLD_BOT | --percentile_bot PERCENTILE_BOT]
                      [--threshold_ppi THRESHOLD_PPI | --percentile_ppi PERCENTILE_PPI] ``

### Example usage:
Demo data of E.coli are included in the GitHub repository.
``python3 ignoromenot.py --ranktable WyattClarkIC-perprotein.tsv --idtable 511145.protein.aliases.v11.5.txt --stringppi 511145.protein.links.full.v11.5.txt -ptop 90 -pbot 1 -tppi 850``
This command reads 3 input files, where the genes are ranked based on their Wyatt Clark 

