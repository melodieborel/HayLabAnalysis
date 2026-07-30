
# README

## General info
This repository contains jupyter notebooks to analyze ePhys and calcium imaging data. It uses the [Miniscope pipeline Minian](https://github.com/miniscope/minian)


## Requirements
For easiest use of the notebooks, make sure to have installed on your computer:
- [git](https://git-scm.com/downloads)
- vscode
- python 3.13 (with or without anaconda)


## Installation
To use the notebooks you can clone the repo github either manually from terminal:

```
git clone git@github.com:melodieborel/HayLabAnalysis.git
```

or [using VSCode](https://www.youtube.com/watch?v=bz1KauFlbQI): melodieborel/HayLabAnalysis


This will download all codes into the local folder of your choice. Best now is to create you own branch to not risk screwing up others work. To do so, click at the bottom left on the git button, + create a new branch..., you can give it your name so that everyone know its yours.

> [!Note]
> If you will work closely with someone who already has a branch, it could make sens to create yours from their...


In order to use the notebooks located in the python subfolder you will need to install a virtual python environment containing all the required dependencies. You can do that using venv or conda. If you don't know what it is or don't know what to choose, please use **venv** as it will make a few things easier...


### With venv (recommended)
On the top right, with a notebook open (for instance python/0_RawDATANumpyViewer.ipynb) you have to select a kernel > Python environments > create python environment > venv > (recreate >) python 3.13.*
then select python/requirements.txt

VSCode will create a subfolder .venv, download and install all packages that are needed to use the notebooks.


> [!Note]
> If you have an issue with creating the venv on a mac, please follow the procedure described on the last comment of the following [link](https://github.com/pyFFTW/pyFFTW/issues/314)


### With conda

Create a fresh conda environment that won't screw up any other environment you might use. Make sure to use this command on the Anaconda terminal if conda isn't added to your path:

```
conda create -n minian313 python=3.13
conda activate minian313
conda install -c conda-forge minian
pip install -r python/requirements.txt
python -m ipykernel install --user --name=minian313
```



## update

### with conda
We are not there yet, but it will probably involve commands such as :

```
conda install --file requirements.txt
conda env update --file local.yml --prune
```


