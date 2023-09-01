Results table 
=============================

The results table displays the peptidoforms or PTMs in the dataset together with their statistics. The user can choose to display only the significant features in the table or not with the button above. The significant features are the ones with an adjPval <= significance level, the latter can also be chosen above, with a default value of 5%.

When a peptidoform level analysis is carried out, only peptidoforms are displayed. When a PTM-level analysis is performed, the PTMs are displayed. Upon clicking on a PTM, all peptidoforms associated with that PTM will also be displayed in a table below. Their statistics are also shown, which are calculated by performing a peptidoform-level analysis with all the same options (except summarisation will be set to false).

When clicking a feature in the table, the corresponding feature is highlighted in the volcano plot and a line figure is plotted in the next tab. This will display the PTM (green) together with its parent protein (log transformed intensity values, blue) and its associated peptidoforms. The latter are displayed both with original (log transformed) intensities (dark grey) and normalised intensities (light grey). Peptidoforms (at the normalised level) that were also signficant in a peptidoform level analysis are shown in pink instead of light grey.
For the detail plots to be calculated correctly, it is important to only select one feature at a time.