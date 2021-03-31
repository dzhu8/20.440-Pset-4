# Integrative Analysis of inCITE-seq Data to Reveal Transcription Factor-Gene Relationships

This repository contains files needed to replicate preliminary analyses for the term project for the course 20.440 at MIT. Thus far, these analyses have consisted of pre-processing of snRNA-seq data (scaling as well as selection 
of highly-variable genes) as well as batch correction and subsequent visualization by uniform manifold approximation and projection (UMAP). Preprocessing and visualization are performed with Scanpy (Wolf, F. A., Angerer, P., & 
Theis, F. J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome biology*, 19(1), 1-5), and batch correction is performed using the harmony-pytorch package (Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., ... & Raychaudhuri, S. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature methods*, 16(12), 1289-1296.). 

## Data

The data is a .h5ad format file generated within the Regev lab at the Broad Institute of MIT and Harvard, and is the output of recent experiments characterizing the Intranuclear Cellular Indexing of Transcriptomes and Epitopes (inCITE-seq) platform (Chung, H., Parkhurst, C., Magee, E. M., Phillips, D., Habibi, E., Chen, F., ... & Regev, A. (2021). Simultaneous single cell measurements of intranuclear proteins and gene expression. *bioRxiv*.) for the simultaneous measurement of intranuclear protein levels and gene expression via a combined hashed-antibody and droplet-based scRNA-seq approach. It can be found within this repository in `data/inCITE.h5ad`. This AnnData object (for a primer on AnnData objects: https://anndata.readthedocs.io/en/latest/anndata.AnnData.html) contains several components, including .X (the raw data matrix, which contains scRNA-seq measurements), .obs (annotations of the observations/cells in this case) and .var (annotations of the variables). The .obs object contains information on the levels of the different transcription factors measured (as inCITE-seq involves simultaneous measurement of nuclear protein levels and mRNA counts), the fraction of mitochondrial counts, the number of total counts per cell, among other useful pieces of information. 

## Code

The outputs of the file can be found in `figures`, and the .py file used to generate each of these figures is `src/inCITE_analysis.py`. 

# Reproducing analyses

1. Make a virtual environment with conda (assuming conda has been installed). If not, it can be downloaded with Miniconda (https://docs.conda.io/en/latest/miniconda.html) or Anaconda (https://docs.anaconda.com/anaconda/install/) by following the directions provided within the appropriate URL. 

```
conda create --name <environment_name>
conda activate <environment_name>
```

where 'environment_name' is replaced by the name you choose to give the environment. 

2. The required packages can also be downloaded using the following series of commands. Python/pip should already be installed on the system if conda is working, but if not, the latest version of Python can be found and downloaded here: https://www.python.org/downloads/:

```
conda install numpy
conda install matplotlib
pip3 install scanpy
pip3 install harmony-pytorch
pip3 install leidenalg
```

3. Download everything in the repository. 

4. Navigate to the location of the downloaded file in command prompt (make sure you have navigated to within the folder housing the data/figures/src folders) and run the Python file:

```
python src\\inCITE_analysis.py
```
**NOTE**: there is a chance this error may result on running this command:
OSError: Unable to open file (file signature not found)

This means that the .h5ad file was corrupted on download (this can also be seen if the size of the downloaded file is much smaller than its actual size of ~1.5 GB). In this case, separately download `inCITE.h5ad` from `data` in this repository and manually move it from the Downloads folder to the location of the other downloaded files. After that, run the above command again. 

This Python file will generate three figures, each of them UMAP representations of the scRNA-seq data: one is generated before pre-processing steps, one is generated after pre-processing but before batch correction, and the final is generated after batch correction. 


5. To leave the virtual environment, run:

```
conda deactivate
```

# Directory structure

The structure of this repo loosely follows [Cookie Cutter Data Science](https://drivendata.github.io/cookiecutter-data-science/)'s recommended structure.

```  
|-- README.md
|
|-- data: README with link to data files
|
|-- figures: Folder in which to save generated figures.
|    
|-- src
     |-- inCITE_analysis: script to perform all analyses outlined in the introductory section.
```
