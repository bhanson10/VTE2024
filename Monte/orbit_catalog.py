import Monte as M
import poincare
import json
import pickle


# DATA_BASE = '/group/monte/delivery/poincare/cr3bp_catalog.sqlite'

# # view/print the full content of the database
# poincare.printSummary( {}, DATA_BASE )

# # view/print the content of the database by specifying a three body
# # system and a family of orbits
# # e.g. view/print all L2 butterfly orbits in the Jupiter-Europa system
# poincare.printSummary(
#    {
#       'PRIMARY_BODY': 'JUPITER',
#       'SECONDARY_BODY': 'EUROPA',
#       'TYPE': 'LIBRATION_POINT',
#       'FAMILY': 'L2_BUTTERFLY'
#    },
#    DATA_BASE
# )

# load ephemeris information for the primary and secondary bodies
boa = M.BoaLoad( [
   "/nav/common/import/ephem/jup310.boa",
   "/nav/common/import/ephem/de433.boa",
   ] )
