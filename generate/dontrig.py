import math

Debug = False

def xsin(angle:float) -> float:
    angle = math.fmod(angle, 2*math.pi)
    angleSquared = angle * angle
    sum = 0.0
    fact = 1.0
    pow = angle
    n = 1
    while n <= 99:
        # n = 1, 3, 5, ...
        term = pow / fact
        prev = sum
        sum += term
        if Debug:
            print('xsin: n={:02d}  term={:24.16e}  sum={:20.16f}  fact={:g}  diff={:24.16e}'.format(n, term, sum, fact, sum-prev))
        if prev == sum:
            return sum
        pow *= angleSquared
        fact *= -((n+1) * (n+2))
        n += 2
    raise Exception('xsin({:0.16g}): failure to converge'.format(angle))


def xcos(angle:float) -> float:
    angle = math.fmod(angle, 2*math.pi)
    angleSquared = angle * angle
    sum = 0.0
    fact = 1.0
    pow = 1.0
    n = 0
    while n <= 99:
        # n = 0, 2, 4, ...
        term = pow / fact
        prev = sum
        sum += term
        if Debug:
            print('xcos: n={:02d}  term={:24.16e}  sum={:20.16f}   fact={:g}'.format(n, term, sum, fact))
        if prev == sum:
            return sum
        pow *= angleSquared
        fact *= -((n+1) * (n+2))
        n += 2
    raise Exception('xcos({:0.16g}): failure to converge'.format(angle))
