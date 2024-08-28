import Monte as M

# ==========================================================================
# CALCULATE LEONARD NIMOY'S AGE AT THE FIRST AIRING OF OF STAR TREK
# ==========================================================================

# Leonard Nimoy's Birthday
nimoyBday = M.Epoch( "26-MAR-1931 ET" )

# Date the first episode of Star Trek aired
startrekAirdate = M.Epoch( "8-SEP-1966 ET" )

# Calculate Age
nimoyAge = startrekAirdate - nimoyBday

# Get integer years and days
years = int( nimoyAge.years() )
days = int( ( nimoyAge.years() - years ) * 365 )

# Format and print output
msg  = ( f"Leonard Nimoy was {years} years, {days} days old "
         f"when the first episode of Star Trek aired!" )
print( msg )
