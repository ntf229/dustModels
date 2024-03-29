### Python versions of dust attenuation models

Includes the following dust models:
- Calzetti
- Cardelli (MW extinction)
- Power Law
- Kriek and Conroy
- TEA

All dust models except Calzetti and TEA are two-component dust models in which the young stars are
attenuated by a power law with separate Av and power law index parameters before 
being combined with the old stars and attenuated by the selected dust model.

This module, once added to your python path, can be used as follows:

```python 
import dustModels as dm
dustySpec, attenuationMags = dm.cardelli(wavelengths, youngSpec, oldSpec,
                             Av, mwr, uvb, AvYoung, dustIndexYoung,
                             fracNoDust, fracNoDustYoung)
```

Where dustySpec is the Cardelli attenuated version of the composite spectrum youngSpec+oldSpec 
and attenuationMags is the corresponding attenuation curve.

Wavelengths should be given in Angstroms and spectra in Janskys.
