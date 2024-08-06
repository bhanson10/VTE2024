# get_acceleration.py
#
# Benjamin Hanson, Summer 2024
#
# GOAL: Call Monte to return the acceleration of a state in the PCR3BP system


import Monte as M
import vista
import poincare
import numpy as np
import mpy.io.data as defaultData
import mpy.traj.force.grav.basic as basicGrav

def CR3BP(x, mu):
    r1 = ((x[0] + mu)**2 + (x[1])**2 + (x[2])**2)**1.5
    r2 = ((x[0] - 1 + mu)**2 + (x[1])**2 + (x[2])**2)**1.5
    v1 = x[3]
    v2 = x[4]
    v3 = x[5]
    v4 = 2*x[4]+x[0]-(mu*(x[0]-1+mu)/r2)-((1-mu)*(x[0]+mu)/r1)
    v5 = -2*x[3]+x[1]-(mu*x[1]/r2)-((1-mu)*x[1]/r1)
    v6 = -(mu*x[2]/r2)-((1-mu)*x[2]/r1)
    return [v1, v2, v3, v4, v5, v6]

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

x0_nd = [IC['X'], IC['Y'], IC['Z'], IC['DX'], IC['DY'], IC['DZ']]

# load ephemeris information for the primary and secondary bodies
km = M.UnitDbl(1, "km")
sec = M.UnitDbl(1, "sec")
boa = M.BoaLoad()
M.DefaultGm.addAll(boa)
boa = M.BoaLoad( ["./poincare/jup310.boa", "./poincare/de433.boa"] )

# create a CR3BP universe in Monte
begin_time = M.Epoch('1-jan-2000 00:00:00et')
jupiter_europa_system = poincare.makeCR3BP( boa, 'jupiter', 'europa', epoch = begin_time)

# true state
print("\nState: ", x0_nd)

# analytical derivative
a_adv = CR3BP(x0_nd, jupiter_europa_system.mu)
print("\nAnalytical Derivative: ", a_adv)

# monte derivative
m_x0_nd = poincare.NONDIMENSIONAL_ST(x0_nd[0], x0_nd[1], x0_nd[2], x0_nd[3], x0_nd[4], x0_nd[5])
m_x0_d = poincare.makeDimensional(m_x0_nd, jupiter_europa_system, 'LPO', epoch = M.Epoch('1-jan-2000 00:00:00et'))
poincare.propagate(m_x0_d, IC["PERIOD"]*jupiter_europa_system.Tstar, setMinStep = 0.01 * sec)
jpGrav = M.GravityBoa.read(jupiter_europa_system.boa, 'LPO')
m_accel = ((jpGrav.accel(M.Epoch('1-jan-2000 00:00:00et'), jupiter_europa_system.barycenter, jupiter_europa_system.rotatingFrame)) * km/sec**2)
m_adv = [m_x0_d.vel()[0] * km/sec * jupiter_europa_system.Tstar/jupiter_europa_system.Lstar, m_x0_d.vel()[1] * km/sec * jupiter_europa_system.Tstar/jupiter_europa_system.Lstar, m_x0_d.vel()[2] * km/sec * jupiter_europa_system.Tstar/jupiter_europa_system.Lstar, m_accel[0] * jupiter_europa_system.Tstar**2/jupiter_europa_system.Lstar, m_accel[1] * jupiter_europa_system.Tstar**2/jupiter_europa_system.Lstar, m_accel[2] * jupiter_europa_system.Tstar**2/jupiter_europa_system.Lstar]
print("\nMonte Derivative: ", m_adv)

# # query orbit at varying dates for visualization
# end_time = begin_time + IC["PERIOD"]*jupiter_europa_system.Tstar
# dates = M.Epoch.range(begin_time, end_time, 60 * sec)
# query_LPO = M.TrajQuery(jupiter_europa_system.boa, 'LPO', 'Jupiter Barycenter', 'Jupiter-Europa Rotating Frame')
# states_LPO = []
# for date in dates:
#     state_LPO = query_LPO.state(date)
#     state_LPO = state_LPO.relativeTo(jupiter_europa_system.barycenter, jupiter_europa_system.rotatingFrame)
#     state_LPO = poincare.makeNondimensional(state_LPO, jupiter_europa_system)
#     states_LPO.append(state_LPO)

# with open('jupiter_europa_lpo.txt', 'w') as file:
#     for row in states_LPO:
#         row = [row.x, row.y, row.z, row.dx, row.dy, row.dz]
#         file.write(' '.join(map(str, row)) + '\n')