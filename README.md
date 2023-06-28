### Python version of the dust models as incorporated in FSPS 

For all two-conponent dust models, the young stars are attenuated by a power-law

Includes the following dust models:
- Calzetti
- Cardelli (MW extinction)
- Power Law
- Kriek and Conroy

All dust models except Calzetti are two-component dust models where the young stars are
attenuated by a power law with separate Av and power law index parameters before 
being combined with the old stars and attenuated by the selected dust model.

This module, once added to your python path, can be used as follows:

```python 
import dustModels as dm
dustySpec, attenuationMags = dm.cardelli(wavelengths, youngSpec, oldSpec,
                             Av, mwr, uvb, AvYoung, dustIndexYoung,
                             fracNoDust, fracNoDustYoung)
```

Where dustySpec is the Cardelli attenuated version of spectrum 
and attenuationMags is the corresponding attenuation curve.
