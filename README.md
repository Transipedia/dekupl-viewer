# DE-kupl Annotation Viewer

DE-kupl annotation Viewer is part of the DE-kupl package, and performs to make interpretation of contigs annoted by DE-kupl annotation, using [R Shiny](https://shiny.rstudio.com/) framework.

### Installation

#### Required dependencies

* R (version >= version 3.2.3) with libraries `shiny`, `shinydashboard`, `ggplot2`, `optparse` and `DT`

#### Install from the sources
You just need to clone this project.

## Usage

```
cd src

Rscript app.R -c /path/to/DiffContigsInfos.tsv -s /path/to/sample_conditions_full.tsv
```

Open a browser at this address `http://0.0.0.0:8080` and enjoy it. 

#### Run in virtual env (docker)
You also could run the app in docker virtual environnement.
```
docker build -t dekupl-viewer:0.1 .

docker run --rm -i -v /[partitionA]:/partitionB]/:z -p 8080:8080 dekupl-viewer:0.1 -c /path/to/DiffContigsInfos.tsv -s /path/to/sample_conditions_full.tsv
``` 
### Tutorial & toys

Toy files are available with this repository to test the app.

```
cd src
Rscript app.R -t
```
