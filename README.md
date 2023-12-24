# NakamuraVisiumR-Analysis
An analytical tool to be used with the 10X Visium platform to compare spatial gene expression values.

This tool was designed with the following considerations:
1. A 10X Visium spatial gene expression assay was completed to generate .CLOUPE file(s).
2. Using the Loupe Browser application from 10X Genomics, a user exported .CSV files of all the gene expression levels in their tissue regions of interest.
3. A user wants to know which genes are differentially expressed between 2 groups within defined regions of interest.

Getting Started:
1. Make sure that your computer has the following installed: R (this was tested on version 4.3.0) and RStudio (this was tested on version 2023.06.0+421)
2. Download the whole repository as a .ZIP file, and extract the contents of the file into a working directory of your choice.
3. Copy the CSV files exported from the Loupe Browser application into the included "CSV" directory.
4. Fill out the "SampleInfo.xlsx" file with the requested information 
- CSV file addresses
- Group names (for example Experimental and Control)
- Tissue region names or cluster names
- Note: it is important that whatever names/labels used in the SampleInfo.xlsx file are also included in the headers of the CSV files. 
5. Install any necessary libraries through CRAN before executing the library commands to load them. 
6. Execute the commands in the Pre-Processing section followed by the commands in the Analysis section. 
7. Example commands are included in the Analysis section, with the expected outputs given in the "Outputs" directory.

