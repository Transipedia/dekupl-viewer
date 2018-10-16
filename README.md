![dekupl-viewer-logo](dekupl-viewer-logo.png)

# DE-kupl Annotation Viewer

DE-kupl annotation Viewer is part of the DE-kupl package, and performs to make interpretation of contigs annoted by DE-kupl annotation, using [R Shiny](https://shiny.rstudio.com/) framework.

![dekupl-viewer-screebbc](dekupl-viewer-logo.png)

## Installation

#### Required dependencies

* R (version >= version 3.3.2) with libraries defined in `install_r_packages.R`

#### Install from the sources
You have to clone this project and run `Rscript install_r_packages.R`.

## Usage

```
cd src

Rscript app.R -c /path/to/DiffContigsInfos.tsv -s /path/to/sample_conditions_full.tsv
```

Open a browser at this address `http://0.0.0.0:8080` and enjoy it.

## Run dekupl-viewer with conda
#### Install conda (miniconda or anaconda)

First you need to install conda, miniconda is harder to use because it comes with nothing installed

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
#### Install dekupl-viewer

```
conda install -n dekupl -y -m --override-channels -c transipedia -c bioconda -c conda-forge -c https://repo.anaconda.com/pkgs/main -c https://repo.anaconda.com/pkgs/free -c https://repo.anaconda.com/pkgs/pro dekupl-viewer
```
This will create a conda environment dekupl (if missing) and install dekupl-viewer inside, the order of the parameters is important.

### Run dekupl-viewer
```
source activate dekupl
dekupl-viewer -c ${PWD}/toy/DiffContigsInfos.tsv -s ${PWD}/toy/sample_conditions_full.tsv
```


## Run in virtual env (docker)
You could also run the app in a docker virtual environnement.
```
docker pull transipedia/dekupl-viewer:latest

docker run --rm -i -v ${PWD}:/${PWD} -p 8080:8080 transipedia/dekupl-viewer:latest -c ${PWD}/toy/DiffContigsInfos.tsv -s ${PWD}/toy/sample_conditions_full.tsv
```

## Run in virtual env (singularity)
You could also run the app in a docker virtual environnement.
```
singularity pull docker://transipedia/dekupl-viewer

./dekupl-viewer.simg -c ${PWD}/toy/DiffContigsInfos.tsv -s ${PWD}/toy/sample_conditions_full.tsv
```

## Toys dataset

Toy files are available with this repository to test the app.

```
cd src
Rscript app.R -t
```

## Profil of Filters

It is possible to define a preset of filters in a TSV file and load it into the interface using option `-f`. Default file is `src/preset-filters/transipedia.tsv`.

`Rscript app.R -h` for help.
