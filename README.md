# msqrob2PTMGUI
R package for launching shiny app to perform differential PTM abundance and usage analysis.

# The package
This R package contains all code to launch a shiny app that contains all functionalities of msqrob2PTM. 
The user can upload a peptide intensity file and a metadata file containing information about the experimental design. 
With the app, the user can perform preprocessing (normalisation, log transformation), visualisation, model building and hypothesis testing.

For more information about the package functionalities, please see our preprint: https://doi.org/10.1101/2023.07.05.547780 
For more information about the included example data, please see https://www.nature.com/articles/s41597-022-01736-1 

# How to install
In RStudio (or another IDE), use 

```
if (!require("devtools")) {install.packages("devtools")}
devtools::install_github("ndmeulem/msqrob2PTMGUI") 
```

This will install the package and its dependencies.

# How to use
To load package: 

```
library(msqrob2PTMGUI)
```

To launch shiny app, use: 

```
launchmsqrob2PTMGUI()
```

# msqrob2 resources
On the following webpage, some examples about different experimental designs are included: https://statomics.github.io/msqrob2Examples/index.html. 
Here, users can see how to formulate the model formula and contrasts for these designs.

