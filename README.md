# SCI Reconstruction via Nonlocal Low Rank Tensor Train and 3-D Total Variation
Code for **Total Variation Regularized Nonlocal Low-rank Tensor Train for Spectral Compressive Imaging**

## Numerical Results
<p align="center">
<img src="https://github.com/han133/NLTT-TV/blob/main/data/results.png?height="300" width="700" raw=true">
</p>
Table 1. Simulation results of HSIs with size of 512×512×31 from CAVE, ICVL and Harvard datasets used in simulation. 

## Usage
- run our algorithm
```
demo.m
```

- run two baseline algorithms: GAP-TV and DeSCI
```
test_desci_cassi.m
```
## References 
- Deep Learning-Based Methods
    - [Convolutional Autoencoder](https://github.com/KAIST-VCLAB/deepcassi)
    - [λ-net](https://github.com/xinxinmiao/lambda-net)
    - [PnP](https://github.com/zsm1211/PnP-CASSI)
