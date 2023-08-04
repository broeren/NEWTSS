# NEWTSS
An algorithm to estimate the error in the wave telescope method as a function of wavevector magnitude.

## Methods
We generated hundreds of random spacecraft configurations who each contain between 4 and 9 spacecraft. We then used each of these configurations to create synthetic timeseries measurements of magnetic fields, as if they were taken from a plane wave. These plane waves varied in direction of travel, velocity, frequency, and phase. Using the magnetic field timeseries, we applied the wave telescope method to estimate the wavevector of the observed plane wave. These estimations had errors that varied with 
- Number of spacecraft in the configuration: $N$
- Shape of spacecraft configuration: $\chi = \sqrt{E^2 + P^2}$ (where $E$ and $P$ are elongation and planarity)
- Relative magnetude of the wavevector: $\bar{k} = kL$ (where $L$ is the size of the spacecraft configuration)

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
