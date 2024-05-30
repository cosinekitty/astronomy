'''dontrig.py
Don Cross <cosinekitty@gmail.com>
https://github.com/cosinekitty/astronomy

These are my own hand-rolled replacements for sine and cosine
on platforms where I'm getting divergent results (order 1e-12).
I want to find out who is correct and who is wrong!
'''
import math

Debug = False

def xsin(angle:float) -> float:
    angle = math.fmod(angle, 2*math.pi)
    numerator = -(angle * angle)
    sum = 0.0
    term = angle
    n = 1
    while n <= 99:
        # n = 1, 3, 5, ...
        prev = sum
        sum += term
        if Debug:
            print('xsin: n={:02d}  term={:24.16e}  sum={:20.16f}  diff={:24.16e}'.format(n, term, sum, sum-prev))
        if prev == sum:
            return sum
        term *= numerator / ((n+1) * (n+2))
        n += 2
    raise Exception('xsin({:0.16g}): failure to converge'.format(angle))


def xcos(angle:float) -> float:
    angle = math.fmod(angle, 2*math.pi)
    numerator = -(angle * angle)
    sum = 0.0
    term = 1.0
    n = 0
    while n <= 99:
        # n = 0, 2, 4, ...
        prev = sum
        sum += term
        if Debug:
            print('xcos: n={:02d}  term={:24.16e}  sum={:20.16f}  diff={:24.16e}'.format(n, term, sum, sum-prev))
        if prev == sum:
            return sum
        term *= numerator / ((n+1) * (n+2))
        n += 2
    raise Exception('xcos({:0.16g}): failure to converge'.format(angle))
