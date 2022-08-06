# IgnoroMeNot
IgnoroMeNot outputs a list of ignorome genes highly associated with other well-annotated genes. Ignorome genes are genes that have little to no validated experimental Gene Ontology annotations. Strong associations between ignorome genes and well-annotated genes can help fulfill the protein function space and minimize the gap in knowledge between "annotation-rich" and "annotation-poor" genes.

### Required files:
IgnoroMeNot requires a tab-separated file of a list of genes with any single annotation metric (e.g. GO annotation count, information content, article count). Highest and lowest annotated genes are defined based on this metric. Below is an example of how the rank table should look like:

| Genes  | Metric            |
|--------|-------------------|
| Q46822 | 55.34890734396416 |
| P18843 | 99.78875455312165 |

Additionally, a STRING interaction network and a protein alias file of the organism to be examined can be downloaded from [string-db.org/cgi/download]. Example: E. coli alias file [https://stringdb-static.org/download/protein.aliases.v11.5/511145.protein.aliases.v11.5.txt.gz] and interaction network [https://stringdb-static.org/download/protein.links.full.v11.5/511145.protein.links.full.v11.5.txt.gz]

## IgnoroMeNot usage:
```
usage: ignoromenot.py [-h] --ranktable RANKTABLE --idtable IDTABLE --stringppi STRINGPPI
                      [--threshold_top THRESHOLD_TOP | --percentile_top PERCENTILE_TOP]
                      [--threshold_bot THRESHOLD_BOT | --percentile_bot PERCENTILE_BOT]
                      [--threshold_ppi THRESHOLD_PPI | --percentile_ppi PERCENTILE_PPI] 

optional arguments:
-h, --help 
```
### Example usage:
Demo data of E.coli are included in the GitHub repository.
``$ python3 ignoromenot.py --ranktable WyattClarkIC-perprotein.tsv --idtable 511145.protein.aliases.v11.5.txt --stringppi 511145.protein.links.full.v11.5.txt --percentile_top 90 --percentile_bot 1 --threshold_ppi 850``

This command reads 3 input files, where the genes coming from E.coli (511145) are ranked based on their Wyatt Clark information content (``WyattClak-per-protein.tsv``). ``--percentile_top 90`` indicates that the genes at the top 10% with respct to WC infromation content are taken,
``--percentile_bot 1`` takes the bottom 1 percent annotated genes based on WC information content (those are the ignorome genes), and ``--threshold_ppi 850`` chooses STRING coexpression scores (``511145.protein.links.full.v11.5.txt``) of 850 and above. The protein alias file (``511145.protein.aliases.v11.5.txt``) makes sure that protein names from different databases have their IDs mapped properly to STRING.
IgnoroMeNot outputs a list of ignorome genes based on these parameters.
