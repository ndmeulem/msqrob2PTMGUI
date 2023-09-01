Metadata input file
======================

What is expected?
------------------------

We expect a file that contains the names of the intensity columns in the first column and in the other columns identifying information. These columns should contain at least enough information so that each row represents one unique ID (without taking the intensity column into account).


Example 
-----------

| Intensity   columns 	|   species 	  |   treatment     	|   sample 	  |
|:-------------------:	|:-------:	|:---------:	|:------:	|
| 210112_hESCs_DDA_P2 	|  human  	|     S     	|    1   	|
| 210114_hESCs_DDA_P3 	|  human  	|     S     	|    2   	|
| 210114_hESCs_DDA_P4 	|  human  	|     S     	|    3   	|
| 210115_hESCs_DDA_P6 	|  human  	|     S     	|    4   	|
| 210115_hESCs_DDA_P8 	|  human  	|     NS     	|    5   	|
| 210114_hESCs_DDA_P5 	|  human  	|     NS     	|    6   	|
| 210115_hESCs_DDA_P7 	|  human  	|     NS     	|    7   	|
| 210118_hESCs_DDA_P1 	|  human  	|     NS     	|    8   	|
| 210118_hESCs_DDA_P9 	|  human  	|     NS     	|    9   	|
| 210112_hESCs_DDA_N1 	|  mouse  	|     S     	|   10   	|
| 210112_hESCs_DDA_N2 	|  mouse  	|     S     	|   11   	|
| 210114_hESCs_DDA_N3 	|  mouse  	|     S     	|   12   	|
| 210114_hESCs_DDA_N4 	|  mouse  	|     S     	|   13   	|
| 210115_hESCs_DDA_N5 	|  mouse  	|     NS     	|   14   	|
| 210115_hESCs_DDA_N6 	|  mouse  	|     NS     	|   15   	|
| 210115_hESCs_DDA_N7 	|  mouse  	|     NS     	|   16   	|
| 210115_hESCs_DDA_N8 	|  mouse  	|     NS     	|   17   	|
| 210118_hESCs_DDA_N9 	|  mouse  	|     NS     	|   18   	|

Options 
---------

- **File separator**: which file separator separates the different columns in your file? Check the right one.


Metadata in SDRF format
==========================

If you have a metadata file in SDRF format, you can check the box to indicate this and the file will be read in with the right parameters.
For more information about the SDRF format, see the following paper:
https://doi.org/10.1038/s41467-021-26111-3


Generate a metadata input file in SDRF format
===============================================

If you do not have a metadatafile yet for your dataset, you can click this button to generate one.
It will then look at your intensity file and extract the names of the columns containing the intensities, i.e. the raw file names.
Hence, it is imperative that you first upload a peptide intensity file and fill out the intensity columns identifier.
The file will then be formatted according to the SDRF format, i.e. a tab separated file with the correct column names.
You will then have to fill out the information for each raw file yourself, after downloading the file to your computer.
At spots to fill out, the file will say "to fill out". This will include information about technical replicates, and the different factor values of the experiment. If there is more than one factor value, you can add an additional column yourself.

After completing the file, you can upload it and check the "metadata in SDRF" box.
