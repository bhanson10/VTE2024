# get_acceleration.py
#
# Benjamin Hanson, Summer 2024
#
# GOAL: Call Monte to return the acceleration of a state in the PCR3BP system

import Monte as M # type: ignore
import vista # type: ignore
import poincare # type: ignore
import numpy as np

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
DATA_BASE = './ephem/cr3bp_catalog.sqlite'

# Load orbit #159 from Jupiter-Europa Low-prograde eastern orbit set
IC = [
      sol for sol in poincare.queryDB(
      {
         'PRIMARY_BODY': 'SATURN',
         'SECONDARY_BODY': 'ENCELADUS',
         'TYPE': 'LIBRATION_POINT',
         'FAMILY': 'DPO',
         'JACOBI_CONSTANTmin': 3.0000780982052,
         'JACOBI_CONSTANTmax': 3.0000780982053,
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
boa = M.BoaLoad( ["./ephem/sat375l.boa", "./ephem/de430.boa"] )

# create a CR3BP universe in Monte
begin_time = M.Epoch('1-jan-2000 00:00:00et')
saturn_europa_system = poincare.makeCR3BP( boa, 'saturn', 'enceladus', epoch = begin_time)
x0_nd[2] = 0
x0_nd[5] = 0

# true state
print("\nState: ", x0_nd)

# analytical derivative
a_adv = CR3BP(x0_nd, saturn_europa_system.mu)
print("\nAnalytical Derivative: ", a_adv)

# monte derivative
m_x0_nd = poincare.NONDIMENSIONAL_ST(x0_nd[0], x0_nd[1], x0_nd[2], x0_nd[3], x0_nd[4], x0_nd[5])
m_x0_d = poincare.makeDimensional(m_x0_nd, saturn_europa_system, 'LPO', epoch = begin_time)
poincare.propagate(m_x0_d, 0.01 * sec)
query_LPO = M.TrajQuery(saturn_europa_system.boa, 'LPO', 'Saturn Barycenter', 'EME2000')
m_x0_d = query_LPO.pva(begin_time, saturn_europa_system.rotatingFrame)
m_adv = [m_x0_d.vel()[0] * km/sec * saturn_europa_system.Tstar/saturn_europa_system.Lstar, m_x0_d.vel()[1] * km/sec * saturn_europa_system.Tstar/saturn_europa_system.Lstar, m_x0_d.vel()[2] * km/sec * saturn_europa_system.Tstar/saturn_europa_system.Lstar, m_x0_d.acc()[0] * km/sec**2 * saturn_europa_system.Tstar**2/saturn_europa_system.Lstar, m_x0_d.acc()[1] * km/sec**2 * saturn_europa_system.Tstar**2/saturn_europa_system.Lstar, m_x0_d.acc()[2] * km/sec**2 * saturn_europa_system.Tstar**2/saturn_europa_system.Lstar]
print("\nMonte Derivative: ", m_adv)

# # query orbit at varying dates for visualization
# end_time = begin_time + IC["PERIOD"]*saturn_europa_system.Tstar
# dates = M.Epoch.range(begin_time, end_time, 60 * sec)
# query_LPO = M.TrajQuery(saturn_europa_system.boa, 'LPO', 'Jupiter Barycenter', 'EME2000')
# states_LPO = []
# for date in dates:
#     state_LPO = query_LPO.state(date)
#     state_LPO = state_LPO.relativeTo(saturn_europa_system.barycenter, saturn_europa_system.rotatingFrame)
#     state_LPO = poincare.makeNondimensional(state_LPO, saturn_europa_system)
#     states_LPO.append(state_LPO)

# with open('saturn_europa_lpo.txt', 'w') as file:
#     for row in states_LPO:
#         row = [row.x, row.y, row.z, row.dx, row.dy, row.dz]
#         file.write(' '.join(map(str, row)) + '\n')