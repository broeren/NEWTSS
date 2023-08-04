# NEWTSS
An algorithm to estimate the error in the wave telescope method as a function of wavevector magnitude.

## Methods


## How to Use
Unaltered, the code will read in the supplied spacecraft positional files from an example (pre Phase-B) 9-spacecraft HelioSwarm configuration. It will then select the optimal subset of spacecraft for a range of wavevectors for one of the configurations (the default is the configuration at hour 94 of the mission). The optimal is selected by minimizing the quantity $\mu_{error} + \alpha \sigma_{error}$ (where $\alpha=2$ by default). To alter which hour configuration is optimized or what this $\alpha$ value is, change these lines of code:
```
hour = 94                       # Configuration hour of HS Phase B DRM mission
std_coef = 2                    # see above
```
The code also reports which magnetudes of wavevectors are possible to reconstruct with errors ($\mu_{error} + \alpha \sigma_{error}$) less than a certain tolerance. To change the value of this tolerance, edit this line of code:
```
Er_thres = 20                   # see above
```

Finally, the code produces a series of figures to visuallize its selection of subset at each scale. The concluding figure, found in the 'Optimal_Ngon' folder, will picture the distribution of error at each scale after selection of the best subset. It also show how many spacecraft ($N$) are in this optimal subset, what the size ($L$) of this subset is, and what it's shape parameter ($\chi$) is.

![alt text](NEWTSS_on_HelioSwarm/figures/Optimal_Ngon/HS_hour94.png)

## References
This work was originally published in the 2023 ApJS paper *Data-driven Uncertainty Quantification of the Wave Telescope Technique: General Equations and Demonstration Using HelioSwarm* by Broeren and Klein
[open access link](https://doi.org/10.3847/1538-4365/acc6c7)
