# Developers quickstart

We are using [Poetry](https://python-poetry.org/) for packaging and dependency management and this is the easiest way to get started. You can install Poetry following [these instructions](https://python-poetry.org/docs/#installation). This presumes you already have Python installed.

After you've installed Poetry, grab a copy of the code from this repo:

```shell
$ git clone https://github.com/microplastics-cluster/fragment-mnp.git
$ cd fragment-mnp
```

Then install the relevant dependencies using Poetry:

```shell
$ poetry install 
```

If you are not already in a virtual environment (i.e. you are not using Conda or virtualenv), this will create a virtual environment in which the package and its dependencies are installed. [These instructions](https://python-poetry.org/docs/basic-usage/#using-your-virtual-environment) detail how you can use this environment. If you are already in a virtual environment, the package and its dependencies will be installed into this.

You can test the installation has worked correctly by running `pytest` in the project root.

```shell
$ pytest
```

## Talking of Conda...

Check out [](conda.md) for specifics of using Conda to manage your environment.