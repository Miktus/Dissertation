import re
import json
import dolang
import subprocess


rr = subprocess.run("pwd")

filename = 'example1.mod'


out = subprocess.check_output(["./dynare_m", f"{filename}", "json=parse", "jsonstdout", "onlyjson"])

print(out.decode("utf8"))


regex = re.compile("\{(.*)\}", re.MULTILINE)

m = regex.search(out.decode('utf8').replace("\n",""))
txt = m.group(0)
data = json.loads(txt)

endo_vars = [e['name'] for e in data['modfile']['endogenous']]
exo_vars = [e['name'] for e in data['modfile']['exogenous']]
parameters = [e['name'] for e in data['modfile']['parameters']]


variables = endo_vars + exo_vars # all variables with time-index

equations = [f"{e['rhs']} - ({e['lhs']})" for e in data['modfile']['model']]
equations = [dolang.symbolic.sanitize(eq, variables) for eq in equations]
equations_strings = [dolang.stringify(eq, variables) for eq in equations]


args = dict([
    ('y_f', [(e,1) for e in endo_vars]),
    ('y', [(e,0) for e in endo_vars]),
    ('y_p', [(e,-1) for e in endo_vars]),
    ('e', [(e,0) for e in exo_vars]),
    ('p', [e for e in parameters])
])

args = dict( [(k,[stringify_symbol(e) for e in v]) for k,v in args.items()] )



from dolang.symbolic import stringify_symbol
from dolang.factory import FlatFunctionFactory

eqs = dict([ (f"equation_{i+1}", eq) for i, eq in enumerate(equations_strings) ])



fff = FlatFunctionFactory(dict(), eqs, args, 'f_dynamic')

from dolang.function_compiler import make_method_from_factory

fun, gufun = make_method_from_factory(fff, vectorize=True, debug=True)


calibration = dict()
data['modfile']['statements']

calibration = dict()
for s in data['modfile']['statements']:
    # parameters
    if s['statementName'] == 'param_init':
        n = s['name']
        v = s['value']
        calibration[n] = float(v)
    elif s['statementName'] == 'init_val':
        for ss in s['vals']:
            n = ss['name']
            v = ss['value']
            calibration[n] = float(v)

import numpy as np
y0 = np.array([calibration[v] for v in endo_vars])
p = np.array([calibration[v] for v in parameters])
e0 = np.array([calibration.get(v,0) for v in exo_vars])

out = np.zeros(len(y0))
res = gufun(y0, y0, y0, e0, p)

from dolo.numeric.serial_operations import numdiff2

ϵ = 1e-9

def numdiff(f,x0):
    y0 = f(x0)
    n = y0.shape[0]
    m = x0.shape[0]
    J = np.zeros((n,m))
    for i in range(m):
        x1 = np.array(x0)
        x1[i] += ϵ
        y1 = f(x1)
        J[:,i] = (y1 - y0)/ϵ
    return J

A = numdiff(lambda x: gufun(x,y0,y0,e0,p), y0)
B = numdiff(lambda x: gufun(y0,x,y0,e0,p), y0)
C = numdiff(lambda x: gufun(y0,y0,x,e0,p), y0)
D = numdiff(lambda x: gufun(y0,y0,y0,x,p), e0)

# these are the matrices from the model:
# A X_{t+1} + B X_t + C X_{t-1} + D E_t
