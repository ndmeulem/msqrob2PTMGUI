Non-enriched dataset input
===============================

Here we expect a non-enriched counterpart of the PTM dataset. It should contain the exact same samples, but representing the background, parent protein status (i.e. non-enriched).

The expectations and options to read in the file are the exact same as for the peptides file.

It is important to note that the intensity column names should be the same as in the peptides file (as they should come from the same sample).

When assessing differential usage, this file will be used to calculate the intensities of the parent proteins. When this file is not provided, but differential usage is selected, the peptides file itself will be used to calculate the intensities of the parent proteins.