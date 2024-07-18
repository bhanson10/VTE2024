# INPUTS:
# 
#   x  ∈ NX1       - true coordinates 
#   dx ∈ NX1       - grid width in N dimensions
#   T  = (σ, b, r) - Lorenz system coefficients 
#
# OUTPUTS
#
#   v ∈ NX1 - advection at coordinates x

def advection(x, dx, T):
    v1 = T[0]*(x[1] - (x[0] + (dx[0]/2.0)))
    v2 = -(x[1] + (dx[1]/2.0)) - x[0]*x[2]
    v3 = -T[1]*(x[2] + (dx[2]/2.0)) + x[0]*x[1] - T[1]*T[2]
    return [v1, v2, v3]