# Developing with Conda/Mamba

To make packaging and publishing easier, we are using [Poetry](https://python-poetry.org/) as our dependency and package manager. You can use Poetry alongside Conda/Mamba, by using a Conda environment as the virtual environment used by Poetry.

When Poetry installs your package's dependencies, it checks whether you are already in a virtual environment (including a Conda environment) and, if so, uses that as the virtual environment for your project. This means that you can use Poetry as a package manager from within a Conda environment.

For example, for this project, you could create a Conda env with just Python 3.9 installed:

```shell
$ conda create -n fmnp python=3.9 --no-default-packages
$ conda activate fmnp
```

You can then install the rest of the dependencies using Poetry:

```shell
(fmnp) $ poetry install
```

Running `poetry env info` will verify that you are using the Conda environment as your virtual environment:

```shell
(fmnp) $ poetry env info

Virtualenv
Python:         3.9.12
Implementation: CPython
Path:           /home/username/miniconda3/envs/fmnp
Valid:          True
```

Note that it is highly recommended that you do not use Conda to install any packages within this environment (except for Python itself), as this is likely to lead to dependency issues.
