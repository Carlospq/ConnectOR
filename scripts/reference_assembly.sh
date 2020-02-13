#!/bin/bash
awk '$3=="transcript" || $3=="exon"{ if ($3=="transcript" &&  $0 !~ "RPKM") print $0" RPKM \"0.1\";"; else print $0 }' $1 | sort -k1,1 -k12,12 -k3,3r -k4,4n > assembly/reference.gtf
