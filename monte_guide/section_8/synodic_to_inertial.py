# synodic_to_inertial.py
#
# Benjamin Hanson, Summer 2024
#
# GOAL: Call Monte to propagate a CR3BP system, and then convert it to inertial frame

import Monte as M # type: ignore
import vista # type: ignore
import poincare # type: ignore
import numpy as np

# load ephemeris information for the primary and secondary bodies
km = M.UnitDbl(1, "km")
sec = M.UnitDbl(1, "sec")
hr = M.UnitDbl(1, "hr")
day = M.UnitDbl(1, "day")
boa = M.BoaLoad()
M.DefaultGm.addAll(boa)
boa = M.BoaLoad( ["./ephem/sat375l.boa", "./ephem/de430.boa"] )

# initial conditions and propagation time
ic = [1.001471995170839, -0.000017518099335, 0, 0.000071987832396, 0.013633926328993, 0]
period = 34 * hr

# create a CR3BP universe in Monte
begin_time = M.Epoch('1-jan-2000 00:00:00et')
saturn_europa_system = poincare.makeCR3BP(boa, 'saturn', 'enceladus', epoch = begin_time)
m_ic = poincare.NONDIMENSIONAL_ST(ic[0], ic[1], ic[2], ic[3], ic[4], ic[5])
m_ic_d = poincare.makeDimensional(m_ic, saturn_europa_system, 'DPO', epoch = begin_time)
poincare.propagate(m_ic_d, period)

# query orbit at varying dates for visualization
end_time = begin_time + period
dates = M.Epoch.range(begin_time, end_time, 60 * sec)
query_DPO = M.TrajQuery(saturn_europa_system.boa, 'DPO', 'Saturn Barycenter', 'EME2000')
query_enc = M.TrajQuery(saturn_europa_system.boa, 'Enceladus', 'Saturn Barycenter', 'EME2000')
rv_DPO_inertial = []
rv_DPO_synodic = []
rv_enc_inertial = []
for date in dates:
    # DPO 
    rv_DPO = query_DPO.state(date)
    rv_DPO_inertial.append(rv_DPO)
    rv_DPO = rv_DPO.relativeTo(saturn_europa_system.barycenter, saturn_europa_system.rotatingFrame)
    rv_DPO = poincare.makeNondimensional(rv_DPO, saturn_europa_system)
    rv_DPO_synodic.append(rv_DPO)

    # Enceladus
    rv_enc = query_enc.state(date)
    rv_enc_inertial.append(rv_enc)

with open('rv_dpo_synodic.txt', 'w') as file:
    for state in rv_DPO_synodic:
        rv= [float(state[0]), float(state[1]), float(state[2]), float(state[3]), float(state[4]), float(state[5])]
        file.write(' '.join(map(str, rv)) + '\n')

with open('rv_dpo_inertial.txt', 'w') as file:
    for state in rv_DPO_inertial:
        rv = [float(state.pos()[0]), float(state.pos()[1]), float(state.pos()[2]), float(state.vel()[0]), float(state.vel()[1]), float(state.vel()[2])]
        file.write(' '.join(map(str, rv)) + '\n')

with open('rv_enc_inertial.txt', 'w') as file:
    for state in rv_enc_inertial:
        rv = [float(state.pos()[0]), float(state.pos()[1]), float(state.pos()[2]), float(state.vel()[0]), float(state.vel()[1]), float(state.vel()[2])]
        file.write(' '.join(map(str, rv)) + '\n')

