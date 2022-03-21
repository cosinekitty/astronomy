# Astronomy Engine module initializer by Don Cross <cosinekitty@gmail.com>.

# Pull Astronomy Engine's public symbols into module scope.
from .astronomy import *

# Delete the redundant nested module namespace.
del astronomy
