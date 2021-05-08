import sympy as sp
from sympy.matrices import Matrix 

a, a_act, aa, appase = sp.symbols("A A* AA* APpase")  # components of the vector
ppase, kcat, kaf, kar, kacat, kpf, kpr, kpcat = sp.symbols("PPase kcat kaf kar kacat kpf kpr kpcat")  # other symbols
f = Matrix([- kcat*a - kaf*a*a_act + kar*aa + kpcat*appase,
            kcat*a - kaf*a*a_act + kar*aa + 2*kacat*aa - kpf*a_act*ppase + kpr*appase, 
            kaf*a*a_act - kar*aa - kacat*aa, 
            kpf*ppase*a_act - kpr*appase - kpcat*appase])
unit = Matrix([a, a_act, aa, appase])
print(f.jacobian(unit))

