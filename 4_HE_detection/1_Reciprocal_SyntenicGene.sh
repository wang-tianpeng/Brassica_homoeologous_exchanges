#!/bin/sh 
set -euo pipefail

### 1. files preparation, including:

SeqBra=Chiifu_v3.1.pep.fa
BedBra=bra31.bed.syn

SeqBol=JZS_0626.pep.fa
BedBol=bol20.bed.syn


### 2. Use SynOrths software to detect syntenic genes

SynOrths -a $SeqBra -b $SeqBol -p $BedBra -q $BedBol -o Bnapus_rapa2ole_synortho_a2c

SynOrths -b $SeqBra -a $SeqBol -q $BedBra -p $BedBol -o Bnapus_ole2rapa_synortho_c2a


### 3. Get the reciprocal Homoeologous gene pairs
perl get_reciprocalMatch_synorth.pl Bnapus_rapa2ole_synortho_a2c Bnapus_ole2rapa_synortho_c2a >Bnapus_ole2rapa_synortho.reciprocalmatch.txt


### Output the syntenic gene pair list: 
### Bnapus_ole2rapa_synortho.reciprocalmatch.txt