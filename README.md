
python -m ipykernel install --user --name=venv311


conda install --file requirements.txt


conda create --name minianMB python=3.11


conda create --name minianMB --clone minion


conda env update --file local.yml --prune





 conda env create -n minian38 -f environment.yml
 conda activate minian38
 conda install conda-forge::ipyfilechooser
 python -m ipykernel install --user --name=minian38
conda install conda-forge::panel   --> to update to panel 0.8.1 which is necessary


==

conda create --name minianDev --clone minian38
python -m ipykernel install --user --name=minianDev


