# LPO_mc.py
#
# Benjamin Hanson, Summer 2024
#
# GOAL: Perform a MC simulation on a Jupiter-Europa LPO


import Monte as M
import poincare
import mpylab
import numpy as np

# JPL Horizons Database
km = M.UnitDbl(1, "km")
sec = M.UnitDbl(1, "sec")
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
boa = M.BoaLoad( ["./poincare/jup310.boa", "./poincare/de433.boa"] )

# create a CR3B universe in Monte
jupiter_europa_system = poincare.makeCR3BP( boa, 'jupiter', 'europa' )

# create a Monte state associated with the IC
IC_nd = poincare.NONDIMENSIONAL_ST(IC['X'], IC['Y'], IC['Z'], IC['DX'], IC['DY'], IC['DZ'])

t = 0 
n = 1000
Tfull = IC['PERIOD']*jupiter_europa_system.Tstar
dt = Tfull/(n-1)
nom_full = []
IC_d = poincare.makeDimensional(IC_nd, jupiter_europa_system, 'LPO')
while t < Tfull:
    dt = min(dt, Tfull-t)
    nom_full.append(IC_d)
    IC_d = poincare.propagate(IC_d, dt, setMinStep = 0.01 * sec)
    t += dt 

t = 0 
T = (1.03784852354)*jupiter_europa_system.Tstar
dt = T/(n-1)
nom = []
IC_d = poincare.makeDimensional(IC_nd, jupiter_europa_system, 'LPO')
while t < T:
    dt = min(dt, T-t)
    nom.append(IC_d)
    IC_d = poincare.propagate(IC_d, dt, setMinStep = 0.01 * sec)
    t += dt 

# Performing MC

# # Non-dimensional IC
# mu = [IC_nd.x, IC_nd.y, IC_nd.z, IC_nd.dx, IC_nd.dy, IC_nd.dz]

# # Covariance corresponding to 100 km uncertainty in position and 10 m/s uncertainty in velocity
# P = [[2.375E-8, 0, 0, 0, 0, 0],
#      [0, 2.375E-8, 0, 0, 0, 0],
#      [0, 0, 0.000000, 0, 0, 0],
#      [0, 0, 0, 5.277E-7, 0, 0], 
#      [0, 0, 0, 0, 5.277E-7, 0],
#      [0, 0, 0, 0, 0, 0.000000]]

# ns = 1000
# samples = np.random.multivariate_normal(mu, P, ns)
# IC_s = []
# x = []
# for i in range(ns):
#     IC_nd = poincare.NONDIMENSIONAL_ST(samples[i][0], samples[i][1], samples[i][2], samples[i][3], samples[i][4], samples[i][5])
#     IC_d = poincare.makeDimensional(IC_nd, jupiter_europa_system, 'LPO')
#     IC_s.append(IC_d)
#     x.append([IC_d])

# t = 0 
# dt = T/(n-1)
# while t < T:
#    dt = min(dt, T-t)
#    for j in range(ns):
#       IC_s[j] = poincare.propagate(IC_s[j], dt, setMinStep = 0.01 * sec)
#       x[j].append(IC_s[j])
#    t += dt 

# mc_x_i = []
# mc_y_i = []
# mc_dx_i = []
# mc_dy_i = []

# mc_x_f = []
# mc_y_f = []
# mc_dx_f = []
# mc_dy_f = []

# for i in range(ns):
#    mc = x[i]
#    mc_x_i.append(mc[0].pos()[0])
#    mc_y_i.append(mc[0].pos()[1])
#    mc_dx_i.append(mc[0].vel()[0])
#    mc_dy_i.append(mc[0].vel()[1])
#    mc_x_f.append(mc[-1].pos()[0])
#    mc_y_f.append(mc[-1].pos()[1])
#    mc_dx_f.append(mc[-1].vel()[0])
#    mc_dy_f.append(mc[-1].vel()[1])
      
# nom_full_x  = [state.pos()[0] for state in nom_full]
# nom_full_y  = [state.pos()[1] for state in nom_full]
# nom_full_dx = [state.vel()[0] for state in nom_full]
# nom_full_dy = [state.vel()[1] for state in nom_full]

# nom_x  = [state.pos()[0] for state in nom]
# nom_y  = [state.pos()[1] for state in nom]
# nom_dx = [state.vel()[0] for state in nom]
# nom_dy = [state.vel()[1] for state in nom]

# europa = [6.6850E5, 0]
# # Get a Monte unit and time system configured figure and axis for the plot.
# fig, ax = mpylab.subplots()

# ax.scatter(mc_x_i, mc_y_i, s = 4, color = 'grey', label = 'MC')
# ax.scatter(europa[0], europa[1], s = 50, color = 'magenta', label = 'Europa')
# ax.scatter(mc_x_f, mc_y_f, s = 4, color = 'grey')
# ax.plot(nom_full_x, nom_full_y, color = 'red', linestyle = 'dashed', label = 'LPO')
# ax.plot(nom_x, nom_y, color = 'red')

# # Add axis labels
# ax.set_xlabel( "X in Jupiter-Europa Synodic Frame (km)", fontsize="14")
# ax.set_ylabel( "Y in Jupiter-Europa Synodic Frame (km)", fontsize="14")

# # Draw the legend
# ax.legend(fontsize="14")
# ax.axis('equal')

# mpylab.show()

# fig2, ax2 = mpylab.subplots()

# ax2.scatter(mc_dx_i, mc_dy_i, s = 4, color = 'grey', label = 'MC')
# ax2.scatter(mc_dx_f, mc_dy_f, s = 4, color = 'grey')
# ax2.plot(nom_full_dx, nom_full_dy, color = 'red', linestyle = 'dashed', label = 'LPO')
# ax2.plot(nom_dx, nom_dy, color = 'red')

# # Add axis labels
# ax2.set_xlabel( "dX in Jupiter-Europa Synodic Frame (km/s)", fontsize="14")
# ax2.set_ylabel( "dY in Jupiter-Europa Synodic Frame (km/s)", fontsize="14")

# # Draw the legend
# ax2.legend(fontsize="14")
# ax2.axis('equal')

# mpylab.show()

# # with open('jupiter_europa_lpo.txt', 'w') as file:
# #     for row in x:
# #         row = [row.pos()[0], row.pos()[1], row.pos()[2], row.vel()[0], row.vel()[1], row.vel()[2]]
# #         file.write(' '.join(map(str, row)) + '\n')
