Model formula
======================

Here, the statistical model formula needs to be specified, which is a symbolic representation of the relationship between variables present in the design. It specifies which factors should be taken into account when analysing the data and how they are related.
It will typically look like this:

**response_variable ~ predictor_variable(s)**

The response variable does not need to be written here. This will typically be the PTM/peptidoform abundance.
The predictor variables do need to be specified (and need to be proceeded by the tilde "~" !) according to lm/lme4 R syntax. The variables that can be included are displayed. 

- Multiple predictor variables can be included; when specifying two, the model formula will look like this:
**~variable1 + variable2**

- Interaction effects can be included:
**~variable1*variable2**
This is the same as:
**~variable1 + variable2 + variable1:variable2**, where variable1:variable2 indicates the interaction between the first and second variable, i.e. the effect of variable1 can be altered according to the value of variable2. 

- Random effects can be included:
**~variable1 + (1|randomvariable2)**

Robust regression
===================

When selecting this option, robust regression is performed, if not an OLS fit is performed. Robust regression using an M estimator is used to make the regression more robust to outliers, by adding weights that will down-weigh observations with high residuals. The weights employed here are Huber weights.
