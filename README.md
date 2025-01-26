<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of Quantlet: CPAoutsample

Published in: submitted to N/A

Description: ‘Performs the Penalized Adaptive Method (PAM), a combination of propagation-separation approach and a SCAD penalty, to forecast a model with possible parameter changes on a given dataset. The model is then used for forecasting over a 1-year horizon. The input data are monthly observations of k-year excess bond risk premia, k = 2, 3, 4, 5, forward rates and other pre-defined macro variables. Computes RMSPE and MAPE for the forecasted values. Plots the time series of the predicted vs. observed values of bond risk excess premia over the prediction interval.’

Keywords: ‘linear model, regression, SCAD penalty, bic, time varying, change point, bootstrap, plot, visualization, nonstationary, financial, prediction’

See also: ‘PAMsimLR, PAMsimCP, PAMCocPia, PAMinsam’

Author: Lenka Zboňáková; Xinjue Li

Submitted: 23. Nov. 2019

Datafile: BRP_data.csv, LN_macrodata_transformed.csv

Input: 
- start.date    : Starting year of observations for PAM and CP
- end.date      : Ending year of observations for PAM and CP
- end.datepred  : Ending year for values used for prediction
- end.dateall   : Ending year of the prediction interval
- start.dateLN  : Starting year of observations for LN
- end.dateLN    : Ending year of observations for LN
- end.dateallLN : Ending year of the prediction interval for LN
- Y.real        : Response variable (k-year excess bond risk premia, k = 2, 3, 4, 5)
- n.years       : Number of years as increment between successive subintervals
- n.boot        : Number of bootstrap loops
- mb.type       : Distribution of multipliers (Bound, Exp or Pois)

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/CPAoutsample/master/PAMoutsam-1.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/CPAoutsample/master/PAMoutsam-2.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/CPAoutsample/master/PAMoutsam-3.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/CPAoutsample/master/PAMoutsam-4.png" alt="Image" />
</div>

