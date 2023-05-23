# Metacognitive Judgments during Visuomotor Learning Reflect the Integration of Error History

`.csv` files are original data. `confidence_datastruct_remodel.m` and `confidence_remodel_online.m` transform the `.csv` files into workable `.mat` dataframes, upon which we subsequently run the analyses (these would be dataConf_EXP1.mat, dataConf_EXP2.mat, dataConf_online_gradual.mat, dataConf_online_zeromean.mat).

`.confidenceModeling_CLEAN.m` and `confidenceModeling_online.m`house the model fitting code. Note the other `.mat` files in this repository are saved data frames from this code. Figures were generated using this code.

`func_conf_*.m` files are the model code.

We also included some further analysis code. `consistency_analysis.m` includes some extraneous analyses, such as basic learning metrics for each experiment. `confRegression.m` includes our analysis of how hand angle and aim changed in response to errors, and whether these changes were sensitive to confidence levels. `makeConfFigureCSVs.m` made the CSVs which were ultimately used for figure generation and also perfoms some aggregation of the confidence model parameters in the form of `param_comp.mat` (FDR-corrected paramter comparisons across experiments) and `param_smry.mat` (summary statistics of each confidence model parameter across experiments. For ease of use, we provide several `.mat` and `.csv` files which are the output of these analysis files including `conf_expX.csv` where X is any number 1-4 corresponding to the experiment number which includes the winning confidence model for each experiment and `confReg.mat` and `confReg.csv` for the regression analysis.

Finally, we would like to thank the makers of the Matlab Violin Plot functions for facilitating the creation of Figure 4.

Bechtold, Bastian, 2016. Violin Plots for Matlab, Github Project
[https://github.com/bastibe/Violinplot-Matlab](https://github.com/bastibe/Violinplot-Matlab), DOI: 10.5281/zenodo.4559847
