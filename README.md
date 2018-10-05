# DE-kupl Annotation Viewer

DE-kupl annotation Viewer is part of the DE-kupl package, and performs to make interpretation of contigs annoted by DE-kupl annotation, using [R Shiny](https://shiny.rstudio.com/) framework.

### Installation

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

#### Run in virtual env (docker)
You could also run the app in a docker virtual environnement.
```
docker pull transipedia/dekupl-viewer:latest

docker run --rm -i -v ${PWD}:/${PWD} -p 8080:8080 transipedia/dekupl-viewer:latest -c ${PWD}/toy/DiffContigsInfos.tsv -s ${PWD}/toy/sample_conditions_full.tsv
``` 
### Toys dataset

Toy files are available with this repository to test the app.

```
cd src
Rscript app.R -t
```

### Profil of Filters

It is possible to define a preset of filters in a TSV file and load it into the interface using option `-f`.

`Rscript app.R -h` for help.