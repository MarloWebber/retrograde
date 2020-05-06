from tools.body import Body
from orbits.orbit import Orbit
import numpy as np
from astropy.time import Time, TimeDelta

body1 = Body(name="Body1",mu = 1, r = 1)

body1r1 = [10.0,10.0,0.0]
body1v1 = [1.0,0.0,0.0]

# orbit1 = Orbit.fromStateVector(body1r1, body1v1, body1, Time('2000-01-01'), "orbit1")



hPa = np.random.uniform(low = 150, high = 1000)
hPe = np.random.uniform(low = 150, high = 1000)
i = 0#np.random.uniform(low=0, high=np.pi)
aPe = np.random.uniform(high = 2*np.pi)
lAn = np.random.uniform(high = 2*np.pi)
tAn = np.random.uniform(high = 2*np.pi)

if hPe > hPa:
    hPe, hPa = hPa, hPe

orbitTest = Orbit.fromApsis(hPe, hPa, i,aPe, lAn,tAn, "test orbit", body1)



print(orbit1)