Preprocessing options
======================

Logarithmic transformation
----------------------------
Oftentimes, intensity data show a very skewed distribution. Logarithmic transformation can help reduce this skewness, and can help to resolve a multiplicative error structure, resulting in a distribution a bit closer to the Gaussian one. This can allow us to then fit a linear model to the data.
There are four options to choose from:
- **none**: no logarithmic transformation will be performed. For example, if the data have already been log transformed previously.
- **log2, log10 or natural**: which base with respect to which logarithms are computed. natural: e = exp(1)

In proteomics, base 2 is often used. When dealing with fold changes, a log2-transformed value of 1 corresponds to a twofold change in the original scale. This can be convenient when interpreting changes in terms of doubling or halving.

Minimum number of non zero columns
------------------------------------
How many intensity values there should at least be present for each peptidoform, i.e. in at least how many samples a quantitative value has been picked up.
Defaults to at least two samples. Any lower than this value might lead to too many missing values in the computations.

Normalisation
---------------
Normalisation to be performed on the entire dataset.  
There are four options to choose from:
- **none**: no normalisation will be performed.
- **center.mean**: mean centering will be performed, i.e. the mean of a variable will be subtracted from all observations on that variable such that the new mean of the variable is zero.
- **center.median**: median centering will be performed, i.e. the median of a variable will be subtracted from all observations on that variable such that the new median of the variable is zero.
- **against protein**: peptidoforms are normalised against their parent protein, to correct for changes occuring in the parent protein. Indeed, changes in the overall protein abundance between conditions can trigger the associated PTM(s) to be detected as differentially abundant. To infer on PTM(s) or peptidoforms for which the effect of the treatment differs from that of the overall protein, we include this normalisation. When this is selected, a median centering to limit technical variation will also be performed first.

Mean and median centering are performed to limit the effect of technical variability in the dataset. Median centering is hereby more robust to outliers than mean centering. When a normalisation against protein is chosen, it allows to assess differential usage (DU), otherwise, differential abundance (DA) can be assessed.

Summarisation
---------------
Should the data be summarised into PTM level? If so, for each unique PTM (= modification, parent protein, location in parent protein combination), the intensities of its peptidoforms will be summarised into one intensity value per sample via robust summarisation. If a peptidoform contains multiple modifications, this peptidform will be used multiple times. Note that it is important that the modification column contains the location information for each modification (e.g. Phospho 214, Phospho 236). In that way, we can correctly distinguish the different PTMs.

When all options are filled out, the user can click the preprocess button and all preprocessing functions will be carried out, as well as the calculation of a density plot and boxplot of the data after preprocessing. There is also a boxplot and density plot present of the original data as it was uploaded, for comparison.
Do note that the preprocessing should be carried out before the user can go on to the next tab.