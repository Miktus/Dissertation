"""
Function implementing the Linear Time Iteration algorithm for solving DSGE models
in the spirit of P. Rendahl (2017)
Final version written by Michal Miktus, April 2019
"""

import numpy as np
from numpy import linalg as npla
from numpy import eye, absolute, amax

def Linear_Time_Iteration(A, B, C, mu, epsilon):
    """
    This function will find the linear time iteration solution to the system of equations in the form of
    AX(-1) + BX + CE[X(+1)] + epsilon = 0
    with a recursive solution in the form of X = FX(-1) + Q*epsilon
    Parameters
    ----------
    A : array_like, dtype=float
        The matrix of coefficients next to endogenous variables entering with a lag
    B : array_like, dtype=float
        The matrix of coefficients next to endogenous, contemporanous variables
    C : array_like, dtype=float
        The matrix of coefficients next to endogenous variables entering with a lead
    mu : number, dtype=float
        Small positive real number to be multiplied by a conformable identity matrix
    epsilon : number, dtype=float
        Threshold value, should be set to a small value like 1e-16
    Returns
    -------
    F : array_like, dtype=float
        The matrix of coefficients next to the endogenous variable in the solution
    Q : array_like, dtype=float
        The matrix of coefficients next to the disturbance term in the solution
    Notes
    -----
    """
    F = 0
    S = 0

    I = eye(*A.shape)*epsilon
    Ch = C
    Bh = (B + 2*C.dot(I))
    Ah = (C.dot(I**2) + B.dot(I) + A)

    if (1/npla.cond(Ah))<1e-16:
        print('Matrix Ah is singular')

    metric = 1

    while metric>epsilon:
        F = -npla.solve((Bh + Ch.dot(F)), Ah)
        S = -npla.solve((Bh + Ah.dot(S)), Ch)
        metric1 = amax(absolute(Ah + Bh.dot(F) + Ch.dot(F.dot(F))))
        metric2 = amax(absolute(Ah.dot(S.dot(S)) + Bh.dot(S) + Ch))
        metric = amax([metric1, metric2])

    eig_F = amax(absolute(npla.eigvals(F)))
    eig_S = amax(absolute(npla.eigvals(S)))

    if (eig_F > 1) or (eig_S > 1) or (mu > 1-eig_S):
        print('Conditions of Proposition 3 violated')

    F = F+I
    Q = -npla.inv(B + C.dot(F))

    return F, Q


