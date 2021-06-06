#!/usr/bin/env python3
#
#   galeqj_matrix.py  -  Don Cross  -  2021-06-06
#
from astronomy import Pivot, Time, Rotation_EQJ_EQD

def PrintRot(rot):
    for i in range(3):
        print('{:20.16f} {:20.16f} {:20.16f}'.format(rot.rot[i][0], rot.rot[i][1], rot.rot[i][2]))

if __name__ == '__main__':
    # This program is a sanity check that I understand galactic coordinates.
    # Compute the GAL/EQJ rotation matrix from
    # the original paper at:
    # https://watermark.silverchair.com/mnras121-0123.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAqowggKmBgkqhkiG9w0BBwagggKXMIICkwIBADCCAowGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQM-h6li-3fYbw16GYKAgEQgIICXRrk54M4D9ltWJzbTRkLJ4Hb3fA4I1rbZZY7AWQlxDU7rZapZ0X7GfQWFlpffxZ-n1MULG121w6MHMZZuugksNGCqhU_8tBrNh3J1WDdLT0kig__6h72D-2ceJwydiARM-6bYd6S-z_CqScXBH6QsU8LflzPnpLz7OtfgEQQVN28ePPO2Box00NC5Sthh0Ouv06e-CyRb3ig2xP-WWP2QeGElCxJicnaDew8BUyofEHKpebe3d-Io_gXZeV7NHV1gIDTXLEp0UpM_KB9nmh39pCr4AW9zd-0j5qUlkC3TqNRlt40WuwDJn-cxFpY4O3hPtBvPmUKqspvtRP6YiqXeuCbQ6uw4OG9xcf7G39D387TeSjgGoeZRPkc0ZhpryawudLiHpXualt_uw8ws8--tAO5aUWVGL_oM2jzLEl7WVvS8jbYxXFY6GRPNBGZeDpYCP0hfXgogG305sLSWYtejlBdYECf1oJ1h_jC7HHfxI2xykD7EwI5tlRPKmNKcADJWS8Wf8zbwIE50ZB-eV1miTpBEFa6KnpD-m34k4-6dNOX9k0ALSMPOOag3o3-RGoeLdMrvqFtENM84Blo-OOADNMU5Mmo2UCjl83onZ83TdlJvSy-e2Z49i1vLScvFZlePdwJArcitb44dR1xfsbsLIR7ApcDCuF_ke96Sc3FkOOkFqM1waQZP2zuiyOr4Jg_6IFyLV5aO-2ocGoY4UkED5RzeNVSIiVC3r2QyhG51lD90ABO7LQ7RWXnsZV915Bkohh5fKIjvNM0ejAm87uMzsYhlzIIind9U6DqnCde

    # Equatorial B1950.0 coordinates of the Galactic North Pole:
    ra = 12.0 + 49.0/60.0
    dec = +27.4

    # Orientation angle for the zero lon/lat point in galactic coordinates:
    theta = 123.0

    # Calculate the B1950.0 epoch based on:
    # https://en.wikipedia.org/wiki/Epoch_(astronomy)
    b1950 = Time.FromTerrestrialTime(2433282.4235 - 2451545.0)
    print('B1950 = {}'.format(b1950))

    # Create a rotation matrix to convert from J2000 equatorial to B1950 equatorial.
    rot = Rotation_EQJ_EQD(b1950)

    # Pivot using the ra/dec angles for the galactic north pole.
    rot = Pivot(rot, 2, -15*ra)
    rot = Pivot(rot, 1, -(90-dec))
    rot = Pivot(rot, 2, theta-180)
    PrintRot(rot)
