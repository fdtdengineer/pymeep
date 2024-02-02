from __future__ import division

import os
import sys
import numpy as np
import meep as mp
import matplotlib.pyplot as plt
from meep import mpb

md = mpb.MPBData(rectify=True, periods=3, resolution=32)
eps = ms.get_epsilon()
converted_eps = md.convert(eps)

plt.imshow(converted_eps.T, interpolation='spline36', cmap='binary')
plt.axis('off')
plt.show()
