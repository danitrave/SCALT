Conda Installation
=====
 
SCARLET is a command-line running tool written in python, R and bash programming language that, currently, runs on UNIX systems only. 
Given that it requires detailed packages with specific versions, we decided to implement SCARLET on a **Conda** environment.

To install Conda, visit the `Conda Installation Page <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ and follow the instructions step by step for the regular installation on your operating system.


SCARLET download
==============

Subsequent to the installation of Conda, make sure that Conda is activated with the following command:

::

  conda activate

Then, download SCARLET either from the following `Github repository <https://github.com/danitrave/SCARLET>`_ as a zipped file or directly using the following command on the command line:

::

  git clone https://github.com/danitrave/SCARLET.git


Move inside the SCARLET directory making use of **cd** command as follows:

::

  cd path/to/directory/SCARLET

Inside the directory, the configuration file **SCARLET_conda_envSetup.yml** enables the installation of the proper SCARLET environment in which all the packages and programs required for the tool are provided. Install the environment running the following command:

::

  conda env create -f SCARLET_conda_envSetup.yml

Finally, make sure to activate the environment if not already activated using the following commad:

::

  conda activate SCARLET

Now, you are ready to play with SCARLET!

