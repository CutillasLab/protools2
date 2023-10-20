# protools2
**A Set Of Tools For Proteomics And Phosphoproteomics Data Analysis**

Reads and normalises proteomic and phosphoproteomic data. Does statistical analysis using `limma` or `t-test`. It also incoroporates functions for KSEA, pathway analysis and network visualisation tools.

## Installation
### RStudio
1. Download the Package Archive File (**protools2_x.x.x.tar.gz**) of the latest [release](https://github.com/iibadshah/protools2/releases/latest).
   - *Not the Source code*
3. In RStudio, click: **Tools** menu
4. Select: **Install Packages...**
5. In the **Install from** list box, select: **Package Archive File (.zip; .tar.gz)**
6. Click: **Browse** to select the downloaded `protools2` Package Archive File
7. Select: **Install**

### R console
#### One-liner
Run: `devtools::install_github("CutillasLab/protools2@*release")`
#### Manual
1. Download the Package Archive File (**protools2_x.x.x.tar.gz**) of the latest [release](https://github.com/iibadshah/protools2/releases/latest).
   - *Not the Source code*
2. Run: `devtools::install_local(path = "C:/path/to/protools2_x.x.x.tar.gz")`
   - *Replace the string argument to* `path` *with the actual location*
   - *You may first need to install `devtools`: `install.packages("devtools")`*

### Errors
If you encounter errors stating that certain required packages are missing, try manually installing the missing packages listed prior to reattempting to install `protools2`:
- `install.packages("package_name")`
- *Replace* `"package_name"` *with each of the packages listed in the R console in turn*
