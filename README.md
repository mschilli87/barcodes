
[//]: # (=====================)

[//]: # ( general information )

[//]: # (=====================)

[//]: # ( file:         README.md                                                                 )

[//]: # ( created:      2016-01-28                                                                )

[//]: # ( last update:  2016-01-28                                                                )

[//]: # ( author(s):    Nikolaos Karaiskos <nikolaos.karaiskos@mdc-berlin.de> (NK),               )

[//]: # (               Marcel Schilling <marcel.schilling@mdc-berlin.de> (MS)                    )

[//]: # ( purpose:      automate building of C program for barcode collapsing for drop-seq data   )


[//]: # (====================================)

[//]: # ( change log (reverse chronological) )

[//]: # (====================================)

[//]: # ( 2016-01-28: adjusted main title to reflect repository name change ('barcodes' -->        )

[//]: # (             'dropseq_barcodes') [MS]                                                     )

[//]: # (             fixed change log comment block (was missing author indication for previous   )

[//]: # (             four changes) [MS]                                                           )

[//]: # (             standardized To-do list style [MS]                                           )

[//]: # (             fixed general information comment block (was missing last update comment)    )

[//]: # (             [MS]                                                                         )

[//]: # (             fixed To-do (still contained part of Pythom To-Do) [MS]                      )

[//]: # (             fixed purpose (still referred to Pythom implementation) [MS]                 )

[//]: # (             replaced Python scripts by C program files / added Markdown 'comments' as    )

[//]: # (             per http://stackoverflow.com/a/32190021 [MS]                                 )

[//]: # (             initial version (purpose, (Python) scripts & To-do list) [NK]                )


[//]: # (==========)

[//]: # ( document )

[//]: # (==========)

# dropseq_barcodes
Collapsing barcodes generated by faulty dropseq beads.
Implemented in C.

## Files
* `collapse_barcodes.c`: collapses barcodes from bulk onto top ones (only mutations & single
deletions are currently considered)
* `Makefile`: builds `collapse_barcodes` from `collapse_barcodes.c`

## To-do list

### High priority
* add UMI correction
* collapse at least 2 deletions before checking mutations
* collapse a set of top barcodes among themselves

### Medium priority
* introduce dynamically thresholds for number of reads and/or percentages of Ts
* introduce check for percentage of Gs in the last positions of barcodes when Ts are found in the
UMIs
* introduce computeKnee function to have an estimate of the number of tops

### Low priority
* allow for more deletions
* standardize/document command line argument parsing
