![dekupl-viewer-logo](dekupl-viewer-logo.png)

[![pipeline status](https://gitlab.com/transipedia/dekupl-viewer/badges/master/pipeline.svg)](https://gitlab.com/transipedia/dekupl-viewer/commits/master) [![docker pull](https://img.shields.io/docker/pulls/transipedia/dekupl-viewer.svg)](https://hub.docker.com/r/transipedia/dekupl-viewer/) [![conda install](https://anaconda.org/transipedia/dekupl-viewer/badges/downloads.svg)](https://anaconda.org/Transipedia/dekupl-viewer)

DE-kupl annotation Viewer is part of the DE-kupl package, and performs to make interpretation of contigs annoted by DE-kupl annotation, using [R Shiny](https://shiny.rstudio.com/) framework.

![dekupl-viewer-screencast](screencast.gif)

- [Usage](#usage)
- [Installation](#installation)
  - [Option 1: Install with conda](#option-1-install-with-conda)
  - [Option 2: Use dekupl-viewer with Docker](#option-2-use-dekupl-viewer-with-docker)
  - [Option 3: Use dekupl-viewer with singularity](#option-3-use-dekupl-viewer-with-singularity)
  - [Option 4: Install from the sources (not recommended)](#option-4-install-from-the-sources-not-recommended)
- [Toy dataset](#toy-dataset)
- [Profile of filters](#profile-of-filters)

## Usage

If dekupl-viewer was installed with the recommended conda package, all you need to do is to activate the dekupl environement and run dekupl-viewer with DiffContigsInfos.tsv file from [dekupl-annot](https://github.com/Transipedia/dekupl-annotation) output and sample_conditions_full.tsv from [dekupl-run](https://github.com/Transipedia/dekupl-run) output.

```
source activate dekupl
dekupl-viewer -c DiffContigsInfos.tsv -s sample_conditions_full.tsv
```

Open a browser at this address `http://0.0.0.0:8080` and enjoy it.

## Installation

### Option 1: Install with conda

- **Step 1: Install conda.** If you do not have a conda distribution installed, we recommend to install miniconda as follows. See [Miniconda website](https://conda.io/miniconda.html) for other installation instructions (ex. for OSX).
    ```
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    ``` 
- **Step 2: Install dekupl-viewer**. This will create a dekupl conda environment (if missing) and install dekupl-run inside. The order of parameters is important.
    ```
    conda install -n dekupl -y -m --override-channels -c transipedia \
     -c bioconda -c conda-forge -c https://repo.anaconda.com/pkgs/main \
     -c https://repo.anaconda.com/pkgs/free \
     -c https://repo.anaconda.com/pkgs/pro dekupl-viewer
    ```
- **Step 3: Run dekupl-viewer**. We first activate the conda environement where dekupl-run was installed, then we run the software.
    ```
    source activate dekupl
    dekupl-viewer -c DiffContigsInfos.tsv -s sample_conditions_full.tsv
    ```

### Option 2: Use dekupl-viewer with Docker

- **Step 1: Retrieve the docker image.**
    ```
    docker pull transipedia/dekupl-viewer:latest
    ```
- **Step 2: Run dekupl-viewer**.
  DiffContigsInfo.tsv and samples_conditions_full.tsv must be located in the directory
  that is mounted in the docker.
    ```
    docker run --rm -i -v ${PWD}:/${PWD} -p 8080:8080 \
    transipedia/dekupl-viewer:latest -c ${PWD}/DiffContigsInfos.tsv \
    -s ${PWD}sample_conditions_full.tsv
    ```

### Option 3: Use dekupl-viewer with singularity

One can create a singularity container from the docker image. Two methods are available, they should both work. 

A difference with docker image is that with Singularity, you don't need to mount any volume, but you must have your config.json and your inputs file in the directory where you are running dekupl-run.

```
singularity pull docker://transipedia/dekupl-viewer
./dekupl-viewer.simg -c DiffContigsInfos.tsv -s sample_conditions_full.tsv
```

### Option 4: Install from the sources (not recommended)

- **Step 1: Install dependancies**. Before using Dekupl-run, install these dependencies:
    - R (version >= version 3.3.2)
    - R packages; `Rscript install_r_packages.R`
- **Step 2: Clone this repository including submodules.**
  `git clone --recursive https://github.com/Transipedia/dekupl-viewer.git`
- **Step 3: Launch dekupl-viewer application**
  ```
  cd src
  Rscript app.R -c /path/to/DiffContigsInfos.tsv -s /path/to/sample_conditions_full.tsv
  ```

## Toy dataset

Toy files are available with this repository to test the app.

```
cd src
Rscript app.R -t
```

## Profile of filters

It is possible to define a preset of filters in a TSV file and load it into the interface using option `-f`. Default file is `src/preset-filters/transipedia.tsv`.

`Rscript app.R -h` for help.
