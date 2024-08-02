# NRHO_mc.py
#
# Benjamin Hanson, Summer 2024
#
# GOAL: Perform a MC simulation on a capstone NRHO

import Monte as M
import poincare
import mpy.traj.force.grav.basic as basicGrav
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from tqdm import tqdm  # If you do not have tqdm installed, comment out this line

km = M.UnitDbl(1, "km")
sec = M.UnitDbl(1, "sec")
boa = M.BoaLoad()
M.DefaultGm.addAll(boa)
boa.load("./poincare/de438.boa")
earth_moon_rot_frame = M.BodyPosDirFrame(boa, 'Earth-Moon Rotating Frame', 'EME2000', ['Distant_Past', 'Distant_Future'], 'Moon', 'Earth', False, '')
earth_moon_system = poincare.makeCR3BP(boa, 'earth', 'moon')

# ==========================================================================
# SET UP THE capstone ORBIT
# ==========================================================================

# Select a name for our capstone spacecraft
sc_capstone = "capstone"

# Define the initial state of the spacecraft using Cartesian elements
ic_capstone = M.State(
   boa, sc_capstone, 'Moon',
   M.Cartesian.x(-8.905163842295311e+03 * km),
   M.Cartesian.y(3.976089073032184e+04 * km),
   M.Cartesian.z(-4.031472411377732e+04 * km),
   M.Cartesian.dx(3.700837104575062e-02 * km/sec),
   M.Cartesian.dy(-7.003840793039799e-02 * km/sec),
   M.Cartesian.dz(2.068511731364004e-01 * km/sec)
   )

# Define which forces will act on the spacecraft during propagation.
forces = [M.GravityForce(boa, sc_capstone)]
basicGrav.add(boa, sc_capstone, ["Sun", "Earth", "Moon"])

# Set 8 day propagation period
begin_time = M.Epoch("08-NOV-2023 05:47:12.2950 TAI")
end_time   = M.Epoch("16-NOV-2023 05:47:12.2950 TAI")
dates = M.Epoch.range(begin_time, end_time, 60 * sec)

# Add the initial state to the "IntegState"
integ_state_capstone = M.IntegState(
   boa,           # Model database used in integration
   begin_time,     # Start time
   end_time,       # End time
   [],            # Events to trigger integration end (none)
   sc_capstone,    # Spacecraft name
   'Moon',        # Center body
   'EME2000',     # Input frame
   'EME2000',     # Integration frame
   ic_capstone, # State initial conditions
   forces,        # Forces which act on state
   False,         # Integrate only partial derivatives (false)
   [],            # Parameters to be used in partial derivative calculations (none)
   []             # Partials tolerance scale factors (allows different partial derivatives to have different integration tolerances, none)
   )

# Add state to our propagation manager "IntegSetup"
integ = M.IntegSetup(boa)
integ.add(integ_state_capstone)
prop = M.DivaPropagator(boa, "DIVA", integ)
prop.create(boa, begin_time, end_time)

# Query orbit at varying dates for visualization
query_capstone = M.TrajQuery(boa, sc_capstone, 'Moon', 'EME2000')
nom_capstone = []
for date in dates:
    state_capstone = query_capstone.state(date)
    state_capstone = state_capstone.relativeTo(earth_moon_system.secondary, earth_moon_system.rotatingFrame)
    state = [state_capstone.pos()[0], state_capstone.pos()[1], state_capstone.pos()[2], state_capstone.vel()[0], state_capstone.vel()[1], state_capstone.vel()[2]]
    nom_capstone.append(state)

# ==========================================================================
# SET UP THE MC SIMULATION
# ==========================================================================
mu = [ic_capstone.pos()[0], ic_capstone.pos()[1], ic_capstone.pos()[2], ic_capstone.vel()[0], ic_capstone.vel()[1], ic_capstone.vel()[2]]

cov_idx = 2

# (sigma_r = 10 km, sigma_v = 1 cm/s)
P = [[9.8659878032566640e+01, -2.605376858616041e+00, 6.8952691844135290e-01, -8.694186035449139e-06, -8.793380508256887e-06, 6.0361176221931510e-06],
     [-2.605376858616041e+00, 9.4249661475048440e+01, 1.0166194646785320e+00, -8.518851092953021e-06, 5.4858797506923830e-06, -2.411998094464790e-05],
     [6.8952691844135290e-01, 1.0166194646785320e+00, 9.9338326392522290e+01, 4.8650538306351200e-06, -2.549743749704787e-05, 5.4037720279180660e-06],
     [-8.694186035449139e-06, -8.518851092953021e-06, 4.8650538306351200e-06, 1.0149199725619390e-10, -8.257339400475887e-13, 1.9555216325606350e-12],
     [-8.793380508256887e-06, 5.4858797506923830e-06, -2.549743749704787e-05, -8.257339400475887e-13, 1.0753326839292990e-10, -3.151759320093024e-12],
     [6.0361176221931510e-06, -2.411998094464790e-05, 5.4037720279180660e-06, 1.9555216325606350e-12, -3.151759320093024e-12, 1.0674614595469450e-10]]
if cov_idx == 1:
    pass
elif cov_idx == 2:
    # (sigma_r = 100 km, sigma_v = 0.1 m/s)
    P = list(np.array(P)*4)
elif cov_idx == 3:
    # (sigma_r = 1 km, sigma_v = 1 mm/s)
    P = list(np.array(P)*1e-2)

ns = 5000
ic_samples = np.random.multivariate_normal(mu, P, ns)
mc_dates_str = ["08-NOV-2023 05:47:12.2950 TAI", "09-NOV-2023 20:30:00.0000 TAI", "16-NOV-2023 05:47:12.2950 TAI"]
mc_dates = []
for i in range(len(mc_dates_str)):
    mc_dates.append(M.Epoch(mc_dates_str[i]))

samples_capstone = []
for i in tqdm(range(ns)): # If you do not have tqdm installed, comment out this line and replace with "for i in range(ns):"

    integ.eraseState(sc_capstone)
    M.TrajLegBoa.eraseAll(boa, sc_capstone)
    M.DivaPropagatorBoa.erase(boa, "DIVA")

    state_sample = M.State(
        boa, sc_capstone, 'Moon',
        M.Cartesian.x(ic_samples[i][0] * km),
        M.Cartesian.y(ic_samples[i][1] * km),
        M.Cartesian.z(ic_samples[i][2] * km),
        M.Cartesian.dx(ic_samples[i][3] * km/sec),
        M.Cartesian.dy(ic_samples[i][4] * km/sec),
        M.Cartesian.dz(ic_samples[i][5] * km/sec)
    )

    integ_state_capstone.setIntegState(state_sample)
    integ.add(integ_state_capstone)
    prop = M.DivaPropagator(boa, "DIVA", integ)
    prop.create(boa, begin_time, end_time)

    query_samples = M.TrajQuery(boa, sc_capstone, 'Moon', 'EME2000')
    sample_capstone = []
    for date in mc_dates:
        state_sample = query_samples.state(date)
        state_sample = state_sample.relativeTo(earth_moon_system.secondary, earth_moon_system.rotatingFrame)
        state = [state_sample.pos()[0], state_sample.pos()[1], state_sample.pos()[2], state_sample.vel()[0], state_sample.vel()[1], state_sample.vel()[2]]
        sample_capstone.append(state)

    samples_capstone.append(sample_capstone)

samples_capstone_temp = []
for i in range(len(mc_dates)):
    x = []
    y = []
    z = []
    dx = []
    dy = []
    dz = []
    for j in range(ns):
        x.append(samples_capstone[j][i][0])
        y.append(samples_capstone[j][i][1])
        z.append(samples_capstone[j][i][2])
        dx.append(samples_capstone[j][i][3])
        dy.append(samples_capstone[j][i][4])
        dz.append(samples_capstone[j][i][5])
    
    samples_capstone_temp.append([x, y, z, dx, dy, dz])

samples_capstone = samples_capstone_temp

# Figure 1 - Full Position Trajectory
figp = plt.figure()
axp = figp.add_subplot(111, projection='3d')
axp.scatter(0, 0, 0, s = 40, color = 'grey', label = "Moon")
for i in range(len(mc_dates)):
    if i==0:
        axp.scatter(samples_capstone[i][0], samples_capstone[i][1], samples_capstone[i][2], s = 2, color = 'black', label = "MC")
    else:
        axp.scatter(samples_capstone[i][0], samples_capstone[i][1], samples_capstone[i][2], s = 2, color = 'black')
axp.plot([row[0] for row in nom_capstone], [row[1] for row in nom_capstone], [row[2] for row in nom_capstone], color = 'red', label = "CAPSTONE NRHO")
axp.set_xlabel("X relative to Moon, Synodic Frame (km)" )
axp.set_ylabel("Y relative to Moon, Synodic Frame (km)" )
axp.set_zlabel("Z relative to Moon, Synodic Frame (km)" )
axp.legend()

# Figure 2 - Full Velocity Trajectory
figv = plt.figure()
axv = figv.add_subplot(111, projection='3d')
for i in range(len(mc_dates)):
    if i==0:
        axv.scatter(samples_capstone[i][3], samples_capstone[i][4], samples_capstone[i][5], s = 2, color = 'black', label = "MC")
    else:
        axv.scatter(samples_capstone[i][3], samples_capstone[i][4], samples_capstone[i][5], s = 2, color = 'black')
axv.plot([row[3] for row in nom_capstone], [row[4] for row in nom_capstone], [row[5] for row in nom_capstone], color = 'red', label = "CAPSTONE NRHO" )
axv.set_xlabel("dX relative to Moon, Synodic Frame (km/s)" )
axv.set_ylabel("dY relative to Moon, Synodic Frame (km/s)" )
axv.set_zlabel("dZ relative to Moon, Synodic Frame (km/s)" )
axv.legend()

figmc={}
axmc={}
count = 0
for i in range(len(mc_dates)):
    figmc[count]=plt.figure()
    axmc[count]=figmc[count].add_subplot(111, projection='3d')
    axmc[count].scatter(samples_capstone[i][0], samples_capstone[i][1], samples_capstone[i][2], alpha = 0.5, s = 2, color = 'black', label = "MC")
    axmc[count].plot([row[0] for row in nom_capstone], [row[1] for row in nom_capstone], [row[2] for row in nom_capstone], color = 'red')
    axmc[count].set_xlabel("X relative to Moon, Synodic Frame (km)" )
    axmc[count].set_ylabel("Y relative to Moon, Synodic Frame (km)" )
    axmc[count].set_zlabel("Z relative to Moon, Synodic Frame (km)" )
    axmc[count].set_title(mc_dates_str[i])
    axmc[count].axis("equal")
    axmc[count].set_xlim3d(min(samples_capstone[i][0]), max(samples_capstone[i][0]))
    axmc[count].set_ylim3d(min(samples_capstone[i][1]), max(samples_capstone[i][1]))
    axmc[count].set_zlim3d(min(samples_capstone[i][2]), max(samples_capstone[i][2]))
    axmc[count].legend()

    count+=1 

    figmc[count]=plt.figure()
    axmc[count]=figmc[count].add_subplot(111, projection='3d')
    axmc[count].scatter(samples_capstone[i][3], samples_capstone[i][4], samples_capstone[i][5], alpha = 0.5, s = 2, color = 'black', label = "MC")
    axmc[count].plot([row[3] for row in nom_capstone], [row[4] for row in nom_capstone], [row[5] for row in nom_capstone], color = 'red')
    axmc[count].set_xlabel("dX relative to Moon, Synodic Frame (km/s)" )
    axmc[count].set_ylabel("dY relative to Moon, Synodic Frame (km/s)" )
    axmc[count].set_zlabel("dZ relative to Moon, Synodic Frame (km/s)" )
    axmc[count].set_title(mc_dates_str[i])
    axmc[count].axis("equal")
    axmc[count].set_xlim3d(min(samples_capstone[i][3]), max(samples_capstone[i][3]))
    axmc[count].set_ylim3d(min(samples_capstone[i][4]), max(samples_capstone[i][4]))
    axmc[count].set_zlim3d(min(samples_capstone[i][5]), max(samples_capstone[i][5]))
    axmc[count].legend()

    count+=1

plt.show()