%cd DSGE_estimation

import tensorflow as tf

import numpy as np


# variables with underscores are outside of the graph (regular python and numpy)
_A = np.load("A_value.npy")
_B = np.load("B_value.npy")
_C = np.load("C_value.npy")
n = _A.shape[0]
_X0 = np.eye(n)
_K = 1000 # max number of iterations



# build tensorflow graph
tf.reset_default_graph()

μ = tf.Variable(0.0, dtype=tf.float64)

A = tf.placeholder(name="A", shape=(n,n), dtype=tf.float64)
B = tf.placeholder(name="B", shape=(n,n), dtype=tf.float64)
C = tf.placeholder(name="C", shape=(n,n), dtype=tf.float64)
X0 = tf.placeholder(name="X0", shape=(n,n), dtype=tf.float64)
k0 = tf.constant(0, dtype=tf.int32)
K = tf.constant(_K, dtype=tf.int32)

# The system we try to solve is A X^2 + B*(1+μ) X + C
# The goal is to find the solution and its derivative w.r.t. μ


def loop(Xi, k):
    F = A@Xi + B*(1+μ)
    Xii = tf.matrix_solve(F, -C)
    return Xii, k+1

def cond(Xi, k):
    return k<=K

X, k = tf.while_loop(cond, loop, [X0, k0])

dX_dμ = tf.gradients(X,μ)

# inputs to the graph
arguments = {A: _A, B: _B, C: _C, X0: _X0}
# g = tf.get_default_graph()

with tf.Session() as sess:
    tf.global_variables_initializer().run()
    _Z = sess.run(X, arguments)
    %time _Z = sess.run(X, arguments)

with tf.Session() as sess:
    tf.global_variables_initializer().run()
    _Z, _dZ = sess.run([X,dX_dμ],  arguments) # warmup
    %time _Z, _dZ = sess.run([X,dX_dμ], arguments)



# Test whether the solution solves the system
test = (_A@_Z@_Z + _B@_Z + _C)

abs(test).max()



# Basic numpy version:
import numpy.linalg
def solveit(A,B,C,X0,K):
    X = X0
    for k in range(K):
        F = A@X0 + B
        X = numpy.linalg.solve(F, -C)
    return X

%time _Z_numpy = solveit(_A, _B, _C, _X0, _K);


## Timings On my computer:
# - basic numpy (no differentiation): 34 ms
# - tensorflow version: 17
# - tensorflow + differentiation w.r.t. \lambda: 67 ms


# TODO:
# - implement stopping rule in the loop
# - evaluate radius
