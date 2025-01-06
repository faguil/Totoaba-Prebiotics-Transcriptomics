#!/bin/bash

cat Totoaba_transcriptome_Sizes.fa | awk '{print $1"\tTotoaba\texon\t1\t"$2"\t.\t+\t.\tgene_id \""$1"\""}' > Totoaba_transcriptome_Fake.gtf