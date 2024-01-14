
# README

## General info
This repository contains jupyter notebooks to analyze ePhys and calcium imaging data. It uses the [Miniscope pipeline Minian](https://github.com/melodieborel/minian)


## Installation

```
git clone git@github.com:melodieborel/interfaceJupyter.git
```


> [!CAUTION]
> Hopefully, the folder downloaded will contain a directory name "minian". If this isn't the case, I have to figure out how subtrees work...
And you will have to manually add the subtree with the command:
>
> ```
> git subtree add --prefix minian git@github.com:melodieborel/minian.git python311 --squash
> ```


Once done, we create a fresh conda environment that won't screw up any other environment you might use. Make sure to use this command on the Anaconda terminal if conda isn't added to your path:

```
conda env create -n minian311 -f minian/environment.yml
```


After the environment is created, you can activate it, install an extra package that wasn't included in the minian requirements, and export the environment to be used in jupyter:

```
conda activate minian311
conda install conda-forge::ipyfilechooser
python -m ipykernel install --user --name=minian311
```



## update

We are not there yet, but it will probably involve commands such as :

```
conda install --file requirements.txt
conda env update --file local.yml --prune
```
