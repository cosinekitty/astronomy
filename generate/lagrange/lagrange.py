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

# EMB-centered Moon state vector from JPL Horizons:
#2458932.750000000 = A.D. 2020-Mar-24 06:00:00.0000 TDB 
# X = 2.673075076694347E-03 Y = 1.870517968021612E-04 Z =-1.731476992383482E-04
# VX=-1.967647146277114E-05 VY= 5.085581013947419E-04 VZ= 2.191688329215286E-04
moon_pos = Vector( 2.673075076694347E-03, 1.870517968021612E-04, -1.731476992383482E-04)
moon_vel = Vector(-1.967647146277114E-05, 5.085581013947419E-04,  2.191688329215286E-04)

# EMB-centered L4 state vector from JPL Horizons:
#2458932.750000000 = A.D. 2020-Mar-24 06:00:00.0000 TDB 
# X = 1.229291654700910E-03 Y = 2.252359994857843E-03 Z = 8.460455897080097E-04
# VX=-4.932901244670627E-04 VY= 2.183360032469220E-04 VZ= 1.409671573944061E-04
L4_pos = Vector( 1.229291654700910E-03, 2.252359994857843E-03, 8.460455897080097E-04)
L4_vel = Vector(-4.932901244670627E-04, 2.183360032469220E-04, 1.409671573944061E-04)

moon_norm = Normal(moon_pos, moon_vel)
print('moon_norm = ', moon_norm)

L4_norm = Normal(L4_pos, L4_vel)
print('L4_norm   = ', L4_norm)

# Find the angle between the two normal vectors.
arcmin = ArcminBetween(moon_norm, L4_norm)
print('Angle between vectors = {:0.6f} arcmin'.format(arcmin))
