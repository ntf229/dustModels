### Python version of the dust models as incorporated in FSPS 

For all two-conponent dust models, the young stars are attenuated by a power-law

This module, once added to your python path, can be used as follows:

```python 
import dustModels as dm
dustySpec, attenuationMags = dm.calzetti(wavelengths, spectrum, Av, fracNoDust)
```
