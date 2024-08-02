# get_acceleration.py
#
# Benjamin Hanson, Summer 2024
#
# GOAL: Call Monte to return the acceleration of a state in the PCR3BP system


import Monte as M
import poincare
import mpylab
import numpy as np

def PCR3BP(x, mu):
    r1 = ((x[0] + mu)**2 + (x[1])**2)**1.5
    r2 = ((x[0] - 1 + mu)**2 + (x[1])**2)**1.5
    v1 = x[2]
    v2 = x[3]
    v3 = 2*x[3]+x[0]-(mu*(x[0]-1+mu)/r2)-((1-mu)*(x[0]+mu)/r1)
    v4 = -2*x[2]+x[1]-(mu*x[1]/r2)-((1-mu)*x[1]/r1)
    return [v1, v2, v3, v4]

# JPL Horizons Database
DATA_BASE = './poincare/cr3bp_catalog.sqlite'

# Load orbit #159 from Jupiter-Europa Low-prograde eastern orbit set
IC = [
   sol for sol in poincare.queryDB(
      {
         'PRIMARY_BODY': 'JUPITER',
         'SECONDARY_BODY': 'EUROPA',
         'TYPE': 'LIBRATION_POINT',
         'FAMILY': 'LPO_EASTERN',
         'JACOBI_CONSTANTmin': 3.00354116777217,
         'JACOBI_CONSTANTmax': 3.00354116777219,
      },
      DATA_BASE
   )
][0]

# load ephemeris information for the primary and secondary bodies
km = M.UnitDbl(1, "km")
sec = M.UnitDbl(1, "sec")
boa = M.BoaLoad()
M.DefaultGm.addAll(boa)
boa = M.BoaLoad( ["./poincare/jup310.boa", "./poincare/de433.boa"] )

# create a CR3B universe in Monte
jupiter_europa_system = poincare.makeCR3BP( boa, 'jupiter', 'europa' )

# create a Monte state associated with the IC
IC_nd = poincare.NONDIMENSIONAL_ST(IC['X'], IC['Y'], IC['Z'], IC['DX'], IC['DY'], IC['DZ'])

coord_nd = [IC_nd.x, IC_nd.y, IC_nd.dx, IC_nd.dy]

print("\nState:                   ", coord_nd)
print("\nAnalytical acceleration: ", PCR3BP(coord_nd, jupiter_europa_system.mu))
print("\nMonte Acceleration:      ")