"""
Function implementing the Linear Time Iteration algorithm for solving DSGE models
in the spirit of P. Rendahl (2017)
Final version written by Michal Miktus, April 2019
"""

from torch import eye, abs, max, gesv, inverse, mm, matrix_power, zeros


def Linear_Time_Iteration(A, B, C, F_initial, mu, epsilon):
    """
    This function will find the linear time iteration solution to the system of equations in the form of
    AX(-1) + BX + CE[X(+1)] + epsilon = 0
    with a recursive solution in the form of X = FX(-1) + Q*epsilon
    Parameters
    ----------
    A : torch, array_like, dtype=float
        The matrix of coefficients next to endogenous variables entering with a lag
    B : torch, array_like, dtype=float
        The matrix of coefficients next to endogenous, contemporanous variables
    C : torch, array_like, dtype=float
        The matrix of coefficients next to endogenous variables entering with a lead
    F : torch, array_like, dtype=float
        The initial guess for F
    mu : number, dtype=float
        Small positive real number to be multiplied by a conformable identity matrix
    epsilon : number, dtype=float
        Threshold value, should be set to a small value like 1e-16
    Returns
    -------
    F : torch, array_like, dtype=float
        The matrix of coefficients next to the endogenous variable in the solution
    Q : torch, array_like, dtype=float
        The matrix of coefficients next to the disturbance term in the solution
    Notes
    -----
    """

    F = F_initial
    S = zeros(*A.shape)

    # F.requires_grad_()
    # S.requires_grad_()

    Id = eye(*A.shape) * mu
    Ch = C
    Bh = (B + 2 * mm(C, Id))
    Ah = (mm(C, matrix_power(Id, 2)) + mm(B, Id) + A)

    metric = 1
    iter = 1

    while metric > epsilon:
        if iter % 10000 == 0:
            print(iter)
        F = -gesv(Ah, (Bh + mm(Ch, F)))[0]
        S = -gesv(Ch, (Bh + mm(Ah, S)))[0]
        metric1 = max(abs(Ah + mm(Bh, F) + mm(Ch, (mm(F, F)))))
        metric2 = max(abs(mm(Ah, mm(S, S)) + mm(Bh, S) + Ch))
        metric = max(metric1, metric2)
        iter += 1
        if iter > 1000000:
            break

    # eig_F = max(abs(eig(F)[0]))
    # eig_S = max(abs(eig(S)[0]))
    # eig_stable = max(abs(eig(mm(inverse(mm(Ah, F) + Bh), Ah))[0]))

    # if (eig_F > 1) or (eig_S > 1) or (mu > 1-eig_S):
    #     print('Conditions of Proposition 3 violated')

    # if (eig_F > 1) or (eig_stable > 1):
    #     print('Conditions for stable and unique solution violated')

    F = F + Id
    Q = -inverse(B + mm(C, F))

    return F, Q
