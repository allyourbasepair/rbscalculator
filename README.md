### RBS Calculator

For an explaination of what the RBS calculator is and how it works, please see the comments in `rbs_calculator.go`.

To develop a model of your own, please have a look at the comments in the `model.go` file in the `model` subpackage.



### Statistics compared to the orignal Salis Lab RBS Calc v2.1

The resulting values from our model differs from the original model for the following reasons:
* We may have incorrectly understood some part of the original model (our model was developed through reading research papers and converting them to code)
* Our model uses the [Linearfold](https://github.com/LinearFold/LinearFold) version of [ViennaRNA](https://github.com/ViennaRNA/ViennaRNA) (to speed up the folding process) while the orignial model doesn't. This leads to different folding structures, and hence different free-energy values.
* To find the shine-dalgarno binding site, the original model co-folds the mRNA sequence with the 16s rRNA sequence. The resulting free-energy value includes the hybridization energy of the mRNA-rRNA binding site (the shine-dalgarno site) and the free energies of the intramolecular folding of the regions upstream and downstream of the binding site on the five prime untranslated region of the mRNA. We achieve the same value by accessing the relevant hybridization free-energy value from a lookup table (which contains the hybridization free-energy values of every possilbe mRNA-rRNA combination), and adding the free-energy values that occurs from folding the upstream and downstream regions on the mRNA sequence. We're not sure if our method gives the same results as co-folding two sequences.




#### Statistics from the original Salis Lab RBS Calc v2.1
| SUBGROUP             | 10-fold Error | 2-fold Error | 5-fold Error | AUROC    | CV       | KL-divergence | Max KL-divergence | N-outliers | N-states | Normalized KL-divergence | Pearson R-squared     | Pearson p       | RIG      | RMSE     | Spearman R-squared     | Spearman p       | intercept | slope   | Sequence entropy | MC       |
| -------------------- | ------------- | ------------ | ------------ | -------- | -------- | ------------- | ----------------- | ---------- | -------- | ------------------------ | --------------------- | --------------- | -------- | -------- | ---------------------- | ---------------- | --------- | ------- | ---------------- | -------- |
| Kosuri_PNAS_2013     | 0.993671      | 0.736777     | 0.924164     | 0.774965 | 1.948457 | 0.979609853   | 4.60517019        | 660        | 11       | 0.212719577              | 0.279536              | 0               | 0.351527 | 1.010964 | 0.279536               | 0                | 10.25464  | \-0.175 | 1209.7078        | 3742.154 |
| Goodman_Science_2013 | 0.944269      | 0.453411     | 0.711638     | 0.53305  | 2.019236 | 0.980949194   | 4.60517019        | 185        | 6        | 0.213010411              | 0.006964              | 9.81E-14        | 0.093707 | 1.610614 | 0.006964               | 9.81E-14         | 10.8653   | \-0.037 | 1404.5946        | 526.4802 |
| ALL                  | 0.97032       | 0.602837     | 0.823708     | 0.738169 | 1.981913 | 0.900072073   | 4.60517019        | 2166       | 11       | 0.195448167              | 0.16712               | 0               | 0.221241 | 1.32865  | 0.16712                | 0                | \-        | \-      | 1447.568         | 2818.294 |

#### Statistics from our version of the Salis Lab RBS Calc v2.1
| SUBGROUP                                |   10-fold Error |   2-fold Error |   5-fold Error |    AUROC |       CV |   KL-divergence |   Max KL-divergence |   N-outliers |   N-states |   Normalized KL-divergence |   Pearson R-squared |    Pearson p |        RIG |     RMSE |   Spearman R-squared |   Spearman p |   intercept |   invalid count |   dataset size |       slope |
|----------------------------------------|---------------- |--------------- |--------------- |--------- |--------- |---------------- |-------------------- |------------- |----------- |--------------------------- |-------------------- |------------- |----------- |--------- |--------------------- |------------- |------------ |---------------- |------- |------------ |
| Kosuri_PNAS_2013                        |        0.992202 |       0.47344  |       0.766162 | 0.656487 | 0.979322 |        0.997283 |             4.60517 |           71 |          4 |                   0.216557 |          0.0924922  | 1.04052e-188 |   0.197295 | 1.40919  |           0.0924922  | 1.04052e-188 |    10.2217  |               0 |   8848 | -0.115384   |
| Goodman_Science_2013                    |        1        |       1        |       1        | 0.534554 | 1.27911  |        1.78756  |             4.60517 |           57 |          1 |                   0.388163 |          0.00535533 | 6.80624e-11  | nan        | 0.482264 |           0.00535533 | 6.80624e-11  |     8.17327 |               0 |   7931 |  0.00634291 |
| Kosuri_PNAS_2013 + Goodman_Science_2013 |        0.998689 |       0.592109 |       0.833065 | 0.855363 | 1.12102  |        1.10263  |             4.60517 |          283 |          4 |                   0.239433 |          0.239092   | 0            |   0.267274 | 1.21933  |           0.239092   | 0            |     9.76937 |               0 |  16779 | -0.0836254  |

Statistics have been computed using the file `./model/stats.py` which has been taken from the [SynBioMTS](https://github.com/reisalex/SynBioMTS) with minor modifications.
