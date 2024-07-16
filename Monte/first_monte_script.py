# first_monte_script.py
#
# Benjamin Hanson, Summer 2024
#
# GOAL: Design a phasing orbit which has a specified resonance with respect to a circular reference orbit

import Monte as M
import mpy.io.data as defaultData
import mpy.traj.force.grav.basic as basicGrav
from mpy.units import *
import mpylab
import matplotlib

# ==========================================================================
# SETUP PROJECT DATA
# ==========================================================================

# Define a project database.
boa = M.BoaLoad()

# Load de405 (planetary ephemeris) and the default body and coordinate systems
# into the project database.
defaultData.loadInto( boa, [
   "frame",
   "body",
   "ephem/planet/de405"
   ] )

# ==========================================================================
# SET UP THE REFERENCE SPACECRAFT ORBIT
# ==========================================================================

# Select a name for our reference spacecraft
scReference = "spacestation"

# Define the radial distance of our reference circular orbit
radiusReference = 6778.0 * km

# Define the initial state of the spacecraft using Conic elements
stateReference = M.State(
   boa, scReference, 'Earth',
   M.Conic.semiMajorAxis( radiusReference ),
   M.Conic.eccentricity( 0.0 ),
   M.Conic.inclination( 5 * deg ),
   M.Conic.argumentOfLatitude( 0 * deg),
   M.Conic.longitudeOfNode( 0 * deg),
   M.Conic.trueAnomaly( 0 * deg)
   )

# Define which forces will act on the spacecraft during propagation.
forces = [
   M.GravityForce( boa, scReference ),
   ]

# Build the gravity nodes connecting the spacecraft to the gravitational
# bodies we want active.
basicGrav.add( boa, scReference, [ "Sun", "Earth", "Moon" ] )

# Set up the beginning and end times for our scenario.
beginTime = M.Epoch( "01-JAN-2000 00:00:00 ET" )
endTime = M.Epoch( "10-JAN-2000 00:00:00 ET" )

# Add the initial state to the "IntegState"
integStateReference = M.IntegState(
   boa,         # Model database used in integration
   beginTime,   # Start time
   endTime,     # End time
   [],          # Events to trigger integration end (none)
   scReference, # Spacecraft name
   'Earth',     # Center body
   'EME2000',   # Input frame
   'EME2000',   # Integration frame
   stateReference, # State initial conditions
   forces,      # Forces which act on state
   False,       # Integrate only partial derivatives (false)
   [],          # Parameters to be used in partial derivative calculations (none)
   []           # Partials tolerance scale factors (allows different partial
                # derivatives to have different integration tolerances, none)
   )

# Add state to our propagation manager "IntegSetup"
integ = M.IntegSetup( boa )
integ.add( integStateReference )

# Set up the propagator.
prop = M.DivaPropagator( boa, "DIVA", integ )

prop.create( boa, beginTime, endTime )

# ==========================================================================
# SET UP PHASING ORBIT
# ==========================================================================

# Select a name for our phasing spacecraft
scPhasing = "capsule"

# Define the resonance we want in our phasing orbit
orbsReference = 50.0
orbsPhasing = 51.0
phaseFactor = orbsReference / orbsPhasing

# Calculate the period that will give our phasing orbit the desired resonance
periodReference = M.Conic.period( stateReference )
periodPhasing = periodReference * phaseFactor

# Define the initial state of the phasing spacecraft using the apoapsis
# range (which is just the range of our circular reference orbit) and
# the solved for period.
statePhasing = M.State(
   boa, scPhasing, 'Earth',
   M.Conic.apoapsisRange( radiusReference ),
   M.Conic.period( periodPhasing ),
   M.Conic.inclination( 5 * deg ),
   M.Conic.argumentOfLatitude( 0 * deg),
   M.Conic.longitudeOfNode( 0 * deg),
   M.Conic.trueAnomaly( -180 * deg) # We set this to -180 so that we start
                                    # propagating from apo instead of peri.
   )

# Build a new IntegState and IntegMass for our phasing orbit
forcesPhasing = [ M.GravityForce( boa, scPhasing ) ]
basicGrav.add( boa, scPhasing, [ "Sun", "Earth", "Moon" ] )
integStatePhasing = M.IntegState(
   boa, beginTime, endTime, [], scPhasing, 'Earth',  'EME2000', 'EME2000',
   statePhasing, forcesPhasing, False, [], []
   )

# Add to the same IntegSetup we previously built. This is so that when we
# run the integrator, it will propagate BOTH the reference and phasing orbits
# at the same time.
integ.add( integStatePhasing )

# Define the end time of our trajectory integration to be periodReference *
# orbsReference ( one complete phasing cycle ).
endTime = beginTime + periodReference * orbsReference

# Update the integration end time span on our reference and phasing
# IntegStates.
integStateReference.setEndTime( endTime )
integStatePhasing.setEndTime( endTime )

# Set up the propagator.
prop = M.DivaPropagator( boa, "DIVA", integ )

# Simultaneously propagate the reference and phasing spacecraft orbits,
# and write the resulting trajectories to our project boa database.
prop.create( boa, beginTime, endTime )

# ==========================================================================
# CHECK TO MAKE SURE THE PHASING ORBITS PERIAPSE ISN'T TOO LOW
# ==========================================================================

# Get a TrajQuery for our phasing orbit.
queryPhasing = M.TrajQuery( boa, scPhasing, 'Earth', 'EME2000' )

# Set up and run our apsis search. We will only run the search over one
# period of our phasing orbit. Since we have no non-conservative forces
# acting on our spacecraft, the periapse will be the same for every orbit
# repetition.
periSearch = M.ApsisEvent( queryPhasing, 'PERIAPSIS' )
searchInterval = M.TimeInterval( beginTime, beginTime + periodPhasing )
searchResults = periSearch.search( searchInterval, 1 * sec )

# Get the value for the first periapse found in the search (there should be
# only one). If it's lower than 200 * km altitude (~6578 * km), issue a
# warning and exit.
foundPeriapse = searchResults[0].value()
if foundPeriapse < 6578 * km:
   print( "Warning: Phasing orbit periapse less "
          "than minimum periapse requirement" )
   print( "   Lower periapse limit: {0}".format( periLimit ) )
   print( "   Phasing orbit periapse: {0}".format( periActual ) )
   exit()

# ==========================================================================
# PLOT THE REFERENCE AND PHASING ORBIT IN EARTH CENTERED EME2000
# ==========================================================================

# Get trajectory query which can be used to get the spacecraft state from
# the reference orbit (we already defined one for the phasing orbit).
queryReference = M.TrajQuery( boa, scReference, 'Earth', 'EME2000' )

# Sample the orbits every 60 seconds for the plot
dates = M.Epoch.range( beginTime, endTime, 60 * sec )
statesReference = []
statesPhasing = []
for date in dates:
   # Get the spacecraft states at current date value.
   stateReference = queryReference.state( date )
   statePhasing = queryPhasing.state( date )
   # Collect the states into their respective lists.
   statesReference.append( stateReference )
   statesPhasing.append( statePhasing )

# Get a Monte unit and time system configured figure and axis for the plot.
fig, ax = mpylab.subplots()

# Extract the EME2000 X and Y values from the reference states using
# Python list comprehensions.
xReference = [ state.pos()[0] for state in statesReference ]
yReference = [ state.pos()[1] for state in statesReference ]

# Plot the reference orbit (looking down the EME2000 Z-axis)
ax.plot( xReference, yReference, color = 'red', label = "Reference Orbit" )

# Plot the phasing orbit (looking down the EME2000 Z-axis)
xPhasing = [ state.pos()[0] for state in statesPhasing ]
yPhasing = [ state.pos()[1] for state in statesPhasing ]
ax.plot( xPhasing, yPhasing, color = 'blue', label = "Phasing Orbit" )

# Add Earth for scale
c = matplotlib.patches.Circle( (0, 0), 6378.0, color = 'lightgrey' )
ax.add_patch( c )
ax.annotate( "Earth", (0, 0) )

# Add axis labels
ax.set_xlabel( "X-EME2000 (Km)" )
ax.set_ylabel( "Y-EME2000 (Km)" )
ax.set_title( "{0:2.0f}-{1:2.0f} Phasing Orbit".format( orbsReference, orbsPhasing  ) )

# Draw the legend
ax.legend()

mpylab.show()

# ==========================================================================
# PLOT THE PHASING ORBIT WRT REFERENCE SC FIXED FRAME
# ==========================================================================

# Set up a body-centered coordinate frame for our reference spacecraft. X is
# defined by the center -> spacecraft direction, Z is the spacecrafts orbital
# angular momentum vector, and Y = Z x X and points in the spacecraft
# downtrack direction.
scFrame = M.BodyPosDirFrame( boa,
   "SpaceStation Frame", # Frame name
   'EME2000', # Base frame
    M.TimeInterval( beginTime, endTime ), # Time interval for frame
    scReference, # Body
    'Earth' # Center
    )

queryScFrame = M.TrajQuery( boa, scPhasing, scReference, "SpaceStation Frame")

statesScFrame = []
# Sample the trajectory and place states into list
for date in dates:
   stateScFrame = queryScFrame.state( date )
   statesScFrame.append( stateScFrame )

fig2, ax2 = mpylab.subplots()

radial = [ s.pos()[0] for s in statesScFrame ]
dtrack = [ s.pos()[1] for s in statesScFrame ]

ax2.plot( dtrack, radial )
ax2.set_xlabel( 'Downtrack - Reference SC Frame (Km)' )
ax2.set_ylabel( 'Radial - Reference SC Frame (Km)' )
ax2.set_title( 'Position of Phasing SC with respect to Reference SC over complete phasing cycle' )

mpylab.show()