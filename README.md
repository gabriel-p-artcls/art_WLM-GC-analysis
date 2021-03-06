
# WLM-GC (WLM-1)

Photometric analysis of the WLM globular cluster (also referred as WLM-1) located in the [WLM galaxy](https://en.wikipedia.org/wiki/Wolf%E2%80%93Lundmark%E2%80%93Melotte).
The photometric data used in the analysis comes from the [Hodge et al. (1999)](http://adsabs.harvard.edu/abs/1999ApJ...521..577H) article, and corresponds to observations in the `V,I` filters.

**Fig 1 from Hodge et al (1999).**

![](fig1_hodge.png)



## Processed data

Hodge et al. comment about this data:

> "*For the color-magnitude diagram (CMD) and luminosity function analyses below, roughly the central 20% of the image was analyzed to maximize the signal from the globular while minimizing the background star contamination. In that region, the effect of background stars is negligible.*"

Because of this the entire frame was processed as part of the GC, with no membership analysis applied.

There is also an analysis of the cluster in [Larsen et al. (2014)](https://www.aanda.org/articles/aa/abs/2014/05/aa22672-13/aa22672-13.html),
and a morphological analysis in [Stephens et al. (2005)](https://arxiv.org/abs/astro-ph/0511502).


## Structural analysis

The center of the cluster is `(0:01:49.48, -15:27:30.7), (0.4561667, -15.4585277) deg` in J2000.0 according to [Billet et al. (2001)](http://iopscience.iop.org/article/10.1086/339181/fulltext/). See the GC at
[Aladin](http://aladin.unistra.fr/AladinLite/?target=00+01+49.480-15+27+30.70&fov=0.10&survey=P%2FDSS2%2Fcolor).

The cluster's radii are (using the scale 0.045 arcsec/pixel for [HST-WFC3 (UVIS)](http://www.stsci.edu/hst/wfc3/ins_performance/detectors) observations):

    r_c = 1.09'' (~24 px)
    r_t = 31''   (~690 px)

according to Hodge et al.

We used the entire frame as the cluster region, ie: `r_cl~200 px`.


## Fundamental parameters in literature

According to Hodge et al., its fundamental parameters are:

> "*A best fit to theoretical isochrones indicates that this cluster has a metallicity of [Fe/H]=-1.52+-0.08 and an age of 14.8+-0.6 Gyr, thus indicating that it is similar to normal old halo globulars in our Galaxy. From the fit we also find that the distance modulus of the cluster is 24.73+-0.07 and the extinction is AV=0.07+-0.06.*"

The logarithmic age is thus log(age)\~10.17, and the metallicity z\~0.00046 (using z_0=0.0152).


## Analyis with ASteCA

The binary fraction was fixed to 0. as the quality of the photometry is not enough to fit this parameter. The entire frame was assumed to be composed of cluster stars  (ie: no decontamination process was attempted), and I imposed a maximum magnitude cut at 26.5 to minimize the photometric incompleteness effect. About this Hodge et al says:

> "*The cutoff magnitudes of 27 in V and 26 in I were chosen to minimize the corrections required due to incompleteness*"


### Limited E_BV<0.05

The five fitted parameters were given ranges as follows:

```
Param      min     max
----------------------
z       0.0001  0.0028
log(a)     9.8    10.2
E(BV)       .0     .05
dm         24.    25.5
M       400000  700000
```

The analysis with these limits gives the following parameter estimates:

```
Param      Mean     MAP     Median     Mode     16th     84th     stddev
------------------------------------------------------------------------
z       0.00014 0.00014    0.00014  0.00013  0.00012   0.00015   0.00004
log(a)    10.01   9.943      10.01    9.949    9.943     10.12   0.07767
E(BV)    0.0401  0.0472     0.0472   0.0465   0.0284    0.0476    0.0132
dm       25.077   25.09      25.09   25.088   25.022    25.101   0.08315
M        433333  400000     406061   406061   400287    475526     51373
```

The output images are stored in the output folder `E_BV_005`.


### Limited E_BV<0.2

The analysis with the same limits given above except `E_BV_max = 0.2`, gives the following parameter estimates:

```
Param      Mean     MAP     Median     Mode     16th     84th     stddev
------------------------------------------------------------------------
z       0.00012 0.00012    0.00012  0.00011  0.00011   0.00014   0.00003
log(a)    10.02   10.01      10.01    10.01    10.01     10.04   0.03959
E(BV)     0.073   0.119     0.0743    0.116   0.0172     0.119    0.0431
dm       24.974  24.818     24.932   24.836   24.818     25.16   0.16743
M        475758  487879     487879   487879   423879    495840     45935
```

The output images are stored in the output folder `E_BV_02`.


### Limited A_V<0.2

Five more runs were made on Oct 2019 where the fitted parameters were given the following ranges:

```
Param      min     max
----------------------
z       0.0001  0.0028
log(a)     9.8    10.2
E(BV)       .0    .645
dm         24.    25.5
M       400000  700000
```

We obtained the following estimates:

```
Run Param    Mean      MAP    Median     Mode     16th     84th     stddev
--------------------------------------------------------------------------
01   z    0.00014 0.000135  0.000135  0.00013  0.00012  0.00017    3.1E-05
01 log(a)   10.15    10.16     10.16    10.16    10.15    10.17    0.02791
01  E(BV)  0.0413   0.0402    0.0402   0.0405   0.0362   0.0516     0.0115
01   dm    25.019   24.987    24.987    24.99   24.973   25.087   0.086314
01   M     445640   406310    421070   413640   406310   496670      52850

02   z    0.00017  0.00019   0.00019  0.00019  0.00012  0.00019    3.9E-05
02 log(a)   10.16    10.18     10.18    10.18    10.14    10.18    0.03454
02  E(BV)  0.0397   0.0434    0.0434   0.0436   0.0265   0.0497     0.0142
02   dm    25.056    25.05     25.05   25.052   24.999   25.077   0.088321
02   M     493250   503540    503540   500910   438410   514850      47871

03   z    0.00013  0.00012   0.00012  0.00012  0.00012  0.00016    3.7E-05
03 log(a)   10.17    10.19     10.19    10.18    10.14    10.19    0.03461
03  E(BV)  0.0275   0.0253    0.0253   0.0248   0.0121    0.049     0.0154
03   dm    25.092   25.112    25.112   25.099   25.018   25.112   0.094611
03   M     468990   454650    454650   453640   436470   508390      48856

04   z    0.00014  0.00011   0.00012  0.00011  0.00011  0.00018    3.6E-05
04 log(a)   10.16    10.16     10.16    10.16    10.13    10.19    0.03291
04  E(BV)   0.025   0.0164    0.0164   0.0162   0.0076   0.0507     0.0182
04   dm    25.118   25.154    25.153    25.15   24.988    25.16    0.10003
04   M     476280   480670    480670   479090   424830   510530      47964

05   z    0.00014  0.00012   0.00013  0.00013  0.00012  0.00017    3.5E-05
05 log(a)   10.15    10.16     10.16    10.16    10.16    10.16    0.02167
05  E(BV)  0.0374   0.0412    0.0412   0.0412   0.0243   0.0444     0.0132
05   dm    25.072   24.993    24.993       25   24.993   25.137    0.15008
05   M     458310   420090    420680   424550   420090   523970      61390
```