# MtMYBs
scripts associated with the genome-wide investigation of MYBs in Medicago

## Co-expression analysis ##
This script performs a co-expression analysis to identify genes with a similar expression pattern based on a given bait gene. Expression data of various forms can be subjected to this analysis e.g. raw counts, TPMs, RPKMs, or FPKMs.


```
python coexp_medicago.py --in <FILE> --out <DIR> --exp <FILE>

Mandatory:
  --in      STR         Sequence ID input file
  --out     STR         Output folder
  --exp     STR         Expression data file

Optional:
  --mapping STR         Gene name mapping file
  --ann     STR         Functional annotation file
```

`--in` specifies a text file that contains gene IDs of interest. One ID per line is expected.

`--out` specifies the output folder. This folder will be created if it does not exist already.

`--exp` specifies a expression data file. The first row contains the sample names. The first column contains the gene IDs. These gene IDs need to match the IDs specified in the input file.

`--mapping` specifies a two column text file with a gene ID in the first column and the display name in the second column. These pairs of IDs and names will te used to replace the gene ID with the name when producing outputs.

`--ann` specifies a test file with functional annotation for each gene. The gene ID is expected in the first column with functional annotation in the following. This annotation text will be included in the output file.


### Reference ###


