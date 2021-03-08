# ConnectOR
Multiple Species Orthology Finder


![alt text](https://raw.githubusercontent.com/Carlospq/ConnectOR/master/ConnectOR.png "ConnectOR summary")

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
You can execute ConnectOR from the command line as follows:
   **$ python ConnectOR.py**  
  

