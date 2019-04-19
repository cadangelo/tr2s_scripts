import numpy as np

from pyne.cccc import Atflux
from pyne.mesh import Mesh, NativeMeshTag

sc = [np.linspace(-600, 1200, 46),
      np.linspace(-800, 500, 31),
      np.linspace(-500, 500, 26)]

m = Mesh(structured=True, structured_coords=sc, mats=None)

at = Atflux("atflux")
at.to_mesh(m, "flux")

m.flux = NativeMeshTag(217, float)
m.save("adj_p.h5m")

