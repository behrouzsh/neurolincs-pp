# neurolincs-pp
This repository is for reformating MS analysis output files for data release of Neurolincs proteomics data.
Please see LINCS_VanEyk_reformat_dataoutput_protocol(1).pdf for the details of the procedure.

---
### How to run
So this package needs three different input files, 
1. first is the actual input .xlsx file, for the provided example the input file is "outputmatrix_OpenSWATH.xlsx", it is assumed    that the input data sheet is in "IPSC_DDA_6600_QE_combined_Canon" tab. This input file is in the input folder of the package.

2. A unimod file that maps unimod id of PTMs to mass difference, for this example: "srep07896-s9.xls" in the resource folder.

3. A meta data file that maps sample ids to another id name, for this example this is "Renaming_neuroLINCS.xlsx" in the resource folder.

The output of the package is stored in the output folder and for this example this is "outputmatrix_OpenSWATH.csv".

---
### Using alternative input, resource and output address to run the package
The package can be run with the following command: 
```
python main.py 
```
Alternatively, this command does the same thing and specifies the relative addresses for the input files
```
python main.py -i input/outputmatrix_OpenSWATH.xlsx -o output/outputmatrix_OpenSWATH.csv -m resource/Renaming_neuroLINCS.xlsx -u resource/srep07896-s9.xls -t IPSC_DDA_6600_QE_combined_Canon
```
The definition of flags are: -i for input, -o for output, -m for renaming file(meta data), -u for the unimod spread sheet and -t for the data tab name in the input file.
This package can be run with alternative input, output file addresses or tab names, for example it can be run with:
```
python main.py -i /Users/shamsabz/Dropbox/neurolincs-stuff/outputmatrix_OpenSWATH.xlsx -o /Users/shamsabz/Dropbox/neurolincs-stuff/outputmatrix.csv -m /Users/shamsabz/Dropbox/neurolincs-stuff/Renaming_neuroLINCS.xlsx -u /Users/shamsabz/Dropbox/neurolincs-stuff/srep07896-s9.xls -t IPSC_DDA_6600_QE_combined_Canon 
```
