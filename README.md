# ConnectOR (v2.0 - 28/07/2021)
Multiple Species Orthology Finder

![Image not available](https://raw.githubusercontent.com/Carlospq/ConnectOR/master/raw/ConnectOR.png "ConnectOR summary")

Changes in V2:
* Depuration of code
* Simplification and optimization of gene clustering step
* Enhanced summary plot visualization
* Enhanced summary gene sheet

## Prerequisites:
- Python3:
  * pickle
  * pandas
  * numpy
  * networkx
  * matplotlib
  * json
  * tqdm
 
- R:
  * ggplot2
  * gridExtra
  * ggplotify
  * dplyr
  * scales
  * ggrepel
  
- BedTools

## Set Up:
To execute connector clone this repository locally.
Following folders/files must be in the same path from command is executed:
  - ConnectOR.py
  - config
  - dictionaries
  - ./scripts/
 
## config example (Tab delimited file):
| species | assembly_version |annotation | chainmap |
| --- | --- | --- | --- |
| Human | hg38 | /path/to/human.gtf |  |
| Mouse	| mm10	| /path/to/mouse.gtf |  |
| Zebrafish | danRer10 | /path/to/zebrafish.gtf |  |

   **species**: name given to each species (must be unique)  
   **assembly_version**: ENSEMBL version of the species  
   **annotation [optional]**: path to the gtf (recommended to use full path)  
   **chainmap [optional]**: comma separated path to each chainmap needed for the analysis (eg for human: "path/to/human_to_mouse.chainmap,path/to/human_to_zebrafish.chainmap"  

## Execution command:
usage: ConnectOR.v2.py [-h] [-mM MINMATCH] [-g]

ConnectOR (v.2.0)
By: Carlos Pulido (carlos.pulido@dbmr.unibe.ch)

optional arguments:
  -h, --help            show this help message and exit
  -mM MINMATCH, --minMatch MINMATCH
                        0.N Minimum ratio of bases that must remap in liftOver
                        step.Default: 30 (0.30 minimum ratio)
  -g, --gene            Generate results at gene level along with exon level
                        results (default: False)
  
## Gene Clustering
![Image not available](https://raw.githubusercontent.com/Carlospq/ConnectOR/master/raw/GENCODE_hg38_mm10.png "ConnectOR summary")
