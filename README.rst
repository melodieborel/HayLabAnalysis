HayLabAnalysis
==============

HayLabAnalysis contains everything someone in Hay Lab should need to easily analyze ePhys and calcium imaging data. It defines an experiment class
regrouping all types of acquired data and metadata, deals with all aspects of the analysis workflow, and provides a series of jupyter notebooks.

* Calcium imaging analysis uses the `Miniscope pipeline Minian`_.
* Spike sorting uses the `SpikeInterface`_.

.. _Miniscope pipeline Minian: https://github.com/melodieborel/minian
.. _SpikeInterface: https://spikeinterface.readthedocs.io/en/latest/

Requirements
------------

For easiest use of the notebooks, make sure to have installed on your computer:

* `vscode`_.
* `python`_. 3.11 (with or without anaconda)
* Optionnally (if you wish to use version control or have collaborators rights) `git`_.

.. _vscode: https://code.visualstudio.com/
.. _python: hhttps://realpython.com/installing-python/
.. _git: https://git-scm.com/downloads


Installation
------------
For simplicity, the installation guidelines assume that you are using vscode to open and run the notebooks. If you are using another editor, please adapt the instructions accordingly.

With pip (recommended for non-coders; no version control)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For non-coders, the easiest way to use the notebooks is to install the package directly from github using pip. This way you will always have the latest version of the package and all dependencies will be installed automatically.

#. Create a folder where you want to work and store your analysis.
#. Open vscode and install the Python extension if you haven't already.
#. Make sure to open the folder you created in vscode.
#. Create a virtual environment:
    #. On the command palette (``Cmd+Shift+P`` or ``Ctrl+Shift+P``), type "Python: Create Environment" and select it.
    #. Select "venv" as the environment type.
    #. (Optional) Select "Recreate" if you want to recreate the environment.
    #. (Optional) Select "Python 3.11" as the Python version.
#. Open a terminal in vscode (``Ctrl+Shift+<``)
#. Install the package using the command: 
    .. code-block:: bash
        pip install git+https://github.com/melodieborel/HayLabAnalysis.git


With github
^^^^^^^^^^^
For coders or people who want to contribute to the development of the package, the best way to use the notebooks is to clone the repository from github. This way you will be able to push your changes and create pull requests.

Clone the repository
""""""""""""""""""""

#. Make sure you have a github account and that you have access to the repository.
#. Make sure you have git installed on your computer,  and that you have set up your ssh keys with github. If you haven't done that yet, please follow the instructions on the following link: https://docs.github.com/en/authentication/connecting-to-github-with-ssh
#. Open vscode and install the Python extension if you haven't already.
#. From the explorer view (``Cmd+Shift+E`` or ``Ctrl+Shift+E``), click on "Clone Repository" (or from the command palette ``Cmd+Shift+P`` or ``Ctrl+Shift+P`` and type "Git: Clone").
#. Select "Clone from GitHub" (you might be asked to sign in to github).
#. Enter the repository URL: 
    .. code-block:: bash
       git@github.com:melodieborel/HayLabAnalysis.git

This will download all codes into the local folder of your choice.


Create your own branch
""""""""""""""""""""""
Best now is to create your own branch to not risk screwing up other's work.

.. note::
    If you will work closely with someone who already has a branch, it could make sens to create yours from their... I

#. Click on the branch name at the bottom left of the window (it probably says "main" or "master").
#. Optional yet recommended: in the dropdown menu, select the branch that is likely to be the closest to your work. Click again on the branch name at the bottom left of the window (now it should say the name of the branch you just selected).
#. In the dropdown menu, select "Create new branch".
#. Give your branch a name (e.g. your username) and click "Create".

Create a virtual environment
"""""""""""""""""""""""""""""
#. On the command palette (``Cmd+Shift+P`` or ``Ctrl+Shift+P``), type "Python: Create Environment" and select it.
#. You can use venv or conda as the environment type.
#. Make sure to select a python version ~= 3.11.
#. Install the required packages by selecting the requirements.txt file provided in the repository.

VSCode will create a subfolder .venv, download and install all packages that are needed to use the notebooks. When you open a notebook, it should automatically use the interpreter from the virtual environment.
If not, you can manually select it by clicking on the interpreter name at the top right of the notebook window.

.. note::
    If you have an issue with creating the venv on a mac, please follow the procedure described on the last comment of the following [link](https://github.com/pyFFTW/pyFFTW/issues/314)

Regularly push your modifications to the remote repository
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
.. important::
    To take advantage of version control and to be able to contribute to the development of the package, you need to regularly push your modifications to the remote repository.

On the source control view (``Cmd+Shift+G`` or ``Ctrl+Shift+G``), you can see all the changes you made to the code since your last commit. From there, you can:
#. Stage your changes: select the files you want to include in the commit
#. Commit your changes: provide a commit message and confirm the commit
#. Push your changes: synchronize your branch with the remote repository

Keep up to date
---------------

With pip
^^^^^^^^^
To ensure you have the latest version of the package, you can run the following command in the terminal:
    .. code-block:: bash
       pip install --upgrade git+https://github.com/melodieborel/HayLabAnalysis.git

With github
^^^^^^^^^^^
Make sure to regularly pull the latest changes from the main branch to your branch.
