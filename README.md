# The neuroendocrine stress response at single-cell resolution. (Lopez et al., 2020)

> scRNA-seq analyses

***INSERT ANOTHER GRAPHIC HERE***

[![]()]()

## Structure
```
HPA_reproducibility
|-- README.md
|-- code
|   \-- mast_scripts.py
|-- docker
|   |-- Dockerfile
|   \-- python-packages.txt
\-- notebooks
    |-- adrenal
    |   |-- 01-preprocessing.ipynb
    |   |-- 02-visualization.ipynb
    |   |-- 03-DE.ipynb
    |   \-- 04-ambient.ipynb
    |-- pituitary
    |   |-- 01-preprocessing.ipynb
    |   |-- 02-visualization.ipynb
    |   |-- 03-DE.ipynb
    |   \`-- 04-ambient.ipynb
    |-- pvn
    |   |-- 01-preprocessing.ipynb
    |   |-- 02-visualization.ipynb
    |   |-- 03-DE.ipynb
    |   \`-- 04-ambient.ipynb
    |-- rep_plots
    |   |-- main.ipynb
    |   \-- supplementary.ipynb
    \-- rep_tables
        \-- tables.ipynb
```

## Installation

We provide a Dockerfile with all the packages required to run the analyses. To build the Docker image use the command:
```shell
docker build docker -t hpa:latest
```

After the image is compiled, you can run it interactively using:

```shell
docker run --interactive --tty --name hpa --publish 8888:8888 --volume $HOME:/root/host_home --workdir /root hpa:latest  /bin/bash
```

This will start a container, in order to start a JupyterLab session use the alias:

```shell
jl
```

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**
