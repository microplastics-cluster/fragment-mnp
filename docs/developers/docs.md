# Writing documentation

We aim to capture all functionality of the model in this documentation. As such, any new feature or change to the code should include associated changes to the documentation. The documentation is stored in the `docs` directory.

The documentation is built using [Jupyter Book](https://jupyterbook.org/en/stable/intro.html), which is based on [Sphinx](https://www.sphinx-doc.org/en/master/). For more information on how to write Sphinx/Jupyter Book documentation, visit those links.

The online documentation is automatically built each time a new tag is released. It is recommend that you test the build routinely whilst pushing updates to the code. If you are developing the project using Poetry, and have run `poetry install` (see [](./quickstart.md)), then you will already have Jupyter Book installed. To build the docs, you can run the following from the root directory of the model:

```bash
$ jb build docs/
```

This will build the docs into the `docs/_build/html` directory. You can either open the `index.html` file directly in your browser, or start up a local server to serve the docs:

```bash
$ python -m http.server -d docs/_build/html
```

By default, this command will start a server at port 8000, and so you will be able to load the docs by visiting [localhost:8000](http://localhost:8000).