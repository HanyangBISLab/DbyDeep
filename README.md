# DbyDeep
DbyDeep: Exploration of MS detectable peptides via deep learning  


## Hardware
DbyDeep requires  
* a GPU with CUDA support


## Installation
DbyDeep requires GPU setting and conda environment.  
1. GPU setting for using tensorflow  
DbyDeep was tested on Ubuntu 18.04, CUDA 11.1, CUDNN 8.0.5 with Nvidia RTX 8000 and RTX A6000 graphic cards with the dependencies above.
> nvidia GPU driver (https://www.nvidia.co.kr/Download/Find.aspx)  
CUDA >= 11.0  (https://developer.nvidia.com/cuda-toolkit-archive)
cudnn >= 8.0  (https://developer.nvidia.com/rdp/cudnn-download)
https://www.tensorflow.org/install/source#linux

2. conda environment.  
tensorflow = 2.4.0  
python >=3.6, <=3.8  
> conda env create -f environment.yml  


## Model
DbyDeep assumes your models are in directories that look like this:
> DbyDeep.h5 - a saved keras model and weights  


## Usage
1. dataset  
use ./scripts/dataset.sh file or python script.  
Currently two output formats are supported: a COMET style db_result.tsv and a MSGF+ style db_result.tsv file.
> bash dataset.sh  

> python dbydeep_data.py \
    --save-path /path/to/save/ \
    --protein-fasta /path/to/proteinDB/ \
    --peptide-tsv /path/to/SearchResult/ \
    --tool-name msgfplus  # [msgfplus, comet]

2. prediction  
use ./scripts/model.sh file or python script.  
> bash model.sh

> python dbydeep_model.py \
    --retrain-flag False \
    --data-path ./data/data.csv \
    --model-path ./data/DbyDeep.h5 \
    --save-path ./data/ \
    --job-name data_result

3. Using DbyDeep on your data  
You can retrain DbyDeep model to your own needs.  
> python dbydeep_model.py \
    --retrain-flag True \
    --data-path ./data/data.csv \
    --model-path ./data/DbyDeep.h5 \
    --save-path ./data/ \
    --job-name data_result

Please note: Sequences except 20 amino acids are not supported. Modifications are not supported.  


## Example
Please find an example input file at ./data/test.csv.

> python dbydeep_model.py \
    --retrain-flag False \
    --data-path ./data/test.csv \
    --model-path ./data/DbyDeep.h5 \
    --save-path ./data/ \
    --job-name test_result