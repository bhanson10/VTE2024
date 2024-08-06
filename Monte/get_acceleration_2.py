# get_acceleration.py
#
# Benjamin Hanson, Summer 2024
#
# GOAL: Call Monte to return the acceleration of a state in the PCR3BP system


import Monte as M
import vista
import numpy as np
import poincare
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

# load ephemeris information for the primary and secondary bodies
km = M.UnitDbl(1, "km")
sec = M.UnitDbl(1, "sec")
boa = M.BoaLoad()
defaultData.loadInto(boa, ["body/gm/Jupiter/Jupiter", "body/gm/Jupiter/Europa"])
jupiter_europa_rot_frame = M.BodyPosDirFrame(boa, 'Jupiter-Europa Rotating Frame', 'EME2000', ['Distant_Past', 'Distant_Future'], 'Europa', 'Jupiter Barycenter', False, '')

begin_time = M.Epoch('1-jan-2000 00:00:00et')
Lstar = 668519 * km
Tstar = 48562 * sec
mu_j  = M.GmBoa.read(boa, "Jupiter")
mu_e  = M.GmBoa.read(boa, "Europa")
mu    = mu_e.gm()/(mu_j.gm() + mu_e.gm())
w = np.sqrt((mu_e.gm() + mu_j.gm())/(Lstar**3))

p_j = [(-mu) * Lstar, 0 * km, 0 * km]
v_j = [0 * km/sec, -w * p_j[0], 0 * km/sec]
Jupiter_IC = M.State(p_j, v_j)
p_e = [(1-mu) * Lstar, 0 * km, 0 * km]
v_e = [0 * km/sec, w * p_e[0], 0 * km/sec]
Europa_IC = M.State(p_e, v_e)

M.TwoBodyTraj(boa, "Jupiter", "Jupiter Barycenter", "EME2000", begin_time, Jupiter_IC)
M.TwoBodyTraj(boa, "Europa", "Jupiter Barycenter", "EME2000", begin_time, Europa_IC)
# M.BaryShiftTraj(boa, "Jupiter", "Jupiter Barycenter", "Jupiter-Europa Rotating Frame", ["Europa"])

sc_name = "LPO"
forces = [M.GravityForce(boa, "LPO")]
basicGrav.add(boa, "LPO", ["Jupiter", "Europa"])

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
x0_d = [x0_nd[0]*Lstar, x0_nd[1]*Lstar, x0_nd[2]*Lstar, x0_nd[3]*Lstar/Tstar, x0_nd[4]*Lstar/Tstar, x0_nd[5]*Lstar/Tstar]

ic_LPO = M.State(
   boa, sc_name, 'Jupiter Barycenter',
   M.Cartesian.x(x0_d[0]),
   M.Cartesian.y(x0_d[1]),
   M.Cartesian.z(x0_d[2]),
   M.Cartesian.dx(x0_d[3]),
   M.Cartesian.dy(x0_d[4]),
   M.Cartesian.dz(x0_d[5])
   )

# Set propagation period
end_time   = begin_time + IC["PERIOD"]*Tstar
dates = M.Epoch.range(begin_time, end_time, 60 * sec)

# Add the initial state to the "IntegState"
integ_state_LPO = M.IntegState(
   boa,                             # Model database used in integration
   begin_time,                      # Start time
   end_time,                        # End time
   [],                              # Events to trigger integration end (none)
   sc_name,                         # Spacecraft name
   'Jupiter Barycenter',            # Center body
   'Jupiter-Europa Rotating Frame', # Input frame
   'EME2000',                       # Integration frame
   ic_LPO,                          # State initial conditions
   forces,                          # Forces which act on state
   False,                           # Integrate only partial derivatives (false)
   [],                              # Parameters to be used in partial derivative calculations (none)
   []                               # Partials tolerance scale factors (allows different partial derivatives to have different integration tolerances, none)
   )

vista.gui.start(boa)

# Add state to our propagation manager "IntegSetup"
integ = M.IntegSetup(boa)
integ.add(integ_state_LPO)
prop = M.DivaPropagator(boa, "DIVA", integ)
prop.create(boa, begin_time, end_time)