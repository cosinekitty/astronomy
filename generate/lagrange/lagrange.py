#!/usr/bin/env python3
import math

class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return '[{:19.16f}, {:19.16f}, {:19.16f}]'.format(self.x, self.y, self.z)

def Normal(p, v):
    # Compute the cross product of the position and velocity vectors.
    nx = p.y*v.z - p.z*v.y
    ny = p.z*v.x - p.x*v.z
    nz = p.x*v.y - p.y*v.x
    # Convert the normal to a unit vector.
    mag = math.sqrt(nx*nx + ny*ny + nz*nz)
    return Vector(nx/mag, ny/mag, nz/mag)

def ArcminBetween(a, b):
    # Take the dot product of the two vectors.
    dot = a.x*b.x + a.y*b.y + a.z*b.z
    # The inverse cosine tells us the angle between them.
    radians = math.acos(dot)
    # Convert radians to arcminutes.
    return 60.0 * math.degrees(radians)

# Geocentric Moon state vector from JPL Horizons:
# 2458932.750000000 = A.D. 2020-Mar-24 06:00:00.0000 TDB
# X = 2.705953999132032E-03 Y = 1.893525408300819E-04 Z =-1.752774223513111E-04
# VX=-1.991849279045934E-05 VY= 5.148133848757774E-04 VZ= 2.218646176831359E-04
moon_pos = Vector( 2.705953999132032E-03, 1.893525408300819E-04, -1.752774223513111E-04)
moon_vel = Vector(-1.991849279045934E-05, 5.148133848757774E-04,  2.218646176831359E-04)

# Geocentric L4 state vector from JPL Horizons:
# 2458932.750000000 = A.D. 2020-Mar-24 06:00:00.0000 TDB
# X = 1.262170577138596E-03 Y = 2.254660738885763E-03 Z = 8.439158665950467E-04
# VX=-4.935321457947509E-04 VY= 2.245912867279576E-04 VZ= 1.436629421560135E-04
L4_pos = Vector( 1.262170577138596E-03,  2.254660738885763E-03,  8.439158665950467E-04)
L4_vel = Vector(-4.935321457947509E-04,  2.245912867279576E-04,  1.436629421560135E-04)

moon_norm = Normal(moon_pos, moon_vel)
print('moon_norm = ', moon_norm)

L4_norm = Normal(L4_pos, L4_vel)
print('L4_norm   = ', L4_norm)

# The normal vectors should be the same, but they aren't.
# Find the angle between the two vectors.
arcmin = ArcminBetween(moon_norm, L4_norm)
print('Angle between vectors = {:0.6f} arcmin'.format(arcmin))
