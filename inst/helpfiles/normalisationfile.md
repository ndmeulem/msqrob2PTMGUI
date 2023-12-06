Usage (normalisation) options 
=========================

In the preprocessing step, users have the choice to perform a global normalisation, i.e. on the entire dataset. 

In the data visualisation tab, you can again have the option to normalise. 
If "relative abundances", this option will apply a normalisation to the dataset that contains only the protein you have selected in the dropdown menu. This is, in other words, an extra normalisation step within the protein.
This option determines how the intensities can be perceived. This will only impact the visualisation in this tab.

The choices are:
- **absolute abundances**: this will not perform an extra normalisation, giving the absolute abundances of the peptides.
- **relative abundances**: this will perform an extra normalisation, giving the relative abundances of the peptides. There is a choice between center.median or center.mean. They will respectively scale the samples so that they have the same median or mean. 


How many digits to display
=============================
How many digits after the decimal point should be displayed in the protein table?
This will also be the number of digits that will be displayed in the downloaded file, if you choose to download the table.