# SCI Reconstruction via Nonlocal Low Rank Tensor Train and 3-D Total Variation
Code for Total Variation Regularized Nonlocal Low-rank Tensor Train for Spectral Compressive Imaging

## Usage
- run our algorithm
```
demo.m
```

- run two baseline algorithms: GAP-TV adn DeSCI
```
# Run plug-and-play gap based on 3d TV denoiser
python pnp_gap_HSI_3dtv.py
# Run plug-and-play gap based on FFDNet
python pnp_gap_HSI_ffdnet.py
```
