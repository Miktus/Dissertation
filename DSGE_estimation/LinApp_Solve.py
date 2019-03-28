"""
MATLAB version 1.1 written by Kerk Phillips, April 2014
PYTHON version adapted by Yulong Li, November 2015, based on the previous adaptations of
Uhlig's Toolkit (1999) by Spencer Lyon in May 2012 and later Chase Coleman
Final version written by Michal Miktus, March 2019
"""
from __future__ import division
import scipy as sp
import numpy as np
from numpy import hstack, vstack, zeros, dot, eye, kron
from scipy import linalg as la
from numpy import linalg as npla


def _nullSpaceBasis(A):
    """
    This function will find the basis of the null space of the matrix A.
    Parameters
    ----------
    A : array_like, dtype=float
        The matrix you want the basis for
    Returns
    -------
    vecs : array_like, dtype=float
        A numpy matrix containing the vectors as row vectors.
    Notes
    -----
    If A is an empty matrix, an empty matrix is returned.
    """
    if A.any():
        U, s, Vh = la.svd(A)
        vecs = np.array([])
        toAppend = A.shape[1] - s.size
        s = np.append(s, zeros((1, toAppend)))
        for i in range(0, s.size):
            if s[i] == 0:
                vecs = Vh[-toAppend:, :]
        if vecs.size == 0:
            vecs = zeros((1, A.shape[1]))
        return np.mat(vecs)
    else:
        return zeros((0, 0))


def qzswitch(i, A, B, Q, Z):
    """
    Takes upper-triangular matrices A, B, orthonormal matrices Q,Z, interchanges
    diagonal elements i and i+1 of both A and B, while maintaining
    Q'AZ' and Q'BZ' unchanged.  Does nothing if ratios of diagonal elements
    in A and B at i and i+1 are the same.  Aborts if diagonal elements of
    both A and B are zero at either position.
    Parameters
    ----------
    i : number, dtype=int
        Index (>=1) of the diagonal element to be interchanged
    A : array_like, dtype=float
        The upper triangular matrix of which some diagonal elements are to
        be interchanged
    B : array_like, dtype=float
        The other upper triangular matrix of which some diagonal elements are
        to be interchanged
    Q : array_like, dtype=float
        An orthonormal matrix from the QZ decomposition
    Z : array_like, dtype=float
        An orthonormal matrix from the QZ decomposition
    Returns
    -------
    A : array_like, dtype=float
        Altered A matrix
    B : array_like, dtype=float
        Altered A matrix
    Q : array_like, dtype=float
        Altered Q matrix
    Z : array_like, dtype=float
        Altered Z matrix
    Notes
    -----
    Copyright: C.A. Sims, 1996, Yale University.
    """

    a = A[i-1, i-1]
    d = B[i-1, i-1]
    b = A[i-1, i]
    e = B[i-1, i]
    c = A[i, i]
    f = B[i, i]

    wz = hstack((dot(c, e)-dot(f, b), (dot(c, d)-dot(f, a)).conj().T))
    xy = hstack(((dot(b, d)-dot(e, a)).conj().T, (dot(c, d)-dot(f, a)).conj().T))

    n = np.sqrt(dot(wz, wz.conj().T))
    m = np.sqrt(dot(xy, xy.conj().T))

    if n == 0:
        print("qzswitch(): Inputs unchanged!")
        return A, B, Q, Z
    else:
        wz = wz/n
        xy = xy/m
        wz = vstack((wz, hstack((-wz[1].conj().T, wz[0].conj().T))))
        xy = vstack((xy, hstack((-xy[1].conj().T, xy[0].conj().T))))
        A[i-1:i+1, :] = xy.dot(A[i-1:i+1, :])
        B[i-1:i+1, :] = xy.dot(B[i-1:i+1, :])
        A[:, i-1:i+1] = A[:, i-1:i+1].dot(wz)
        B[:, i-1:i+1] = B[:, i-1:i+1].dot(wz)
        Z[:, i-1:i+1] = Z[:, i-1:i+1].dot(wz)
        Q[i-1:i+1, :] = xy.dot(Q[i-1:i+1, :])
    return A, B, Q, Z


def qzdiv(stake, A, B, Q, Z):
    """
    Takes upper-triangular matrices A, B, orthonormal matrices Q,Z, rearranges them
    so that all cases of abs(B(i,i)/A(i,i))>stake are in lower right
    corner, while preserving upper-triangular and orthonormal properties and Q'AZ' and Q'BZ'.

    Parameters
    ----------
    stake : number, dtype=float
    A : array_like, dtype=float
        An upper triangular matrix
    B : array_like, dtype=float
        An upper triangular matrix
    Q : array_like, dtype=float
        An orthonormal matrix from the QZ decomposition
    Z : array_like, dtype=float
        An orthonormal matrix from the QZ decomposition
    Returns
    -------
    A : array_like, dtype=float
        Rearranged A matrix
    B : array_like, dtype=float
        Rearranged B matrix
    Q : array_like, dtype=float
        Rearranged Q matrix
    Z : array_like, dtype=float
        Rearranged Z matrix
    Notes
    -----
    Copyright: C.A. Sims, 1996, Yale University.
    """

    n, jnk = A.shape

    root = abs(vstack((np.diag(A), np.diag(B))).T)
    tmp = (root[:, 0] < 1.e-13).astype(int)
    root[:, 0] = root[:, 0] - tmp * (root[:, 0]+root[:, 1])
    root[:, 1] = root[:, 1]/root[:, 0]
    for i in xrange(n, 0, -1):
        m = 0
        for j in xrange(i, 0, -1):
            if (root[j-1, 1] > stake or root[j-1, 1] < -.1):
                m = j
                break
        if m == 0:
            print("qzdiv(): Inputs unchanged!")
            return A, B, Q, Z
        for k in xrange(m, i, 1):
            A, B, Q, Z = qzswitch(k, A, B, Q, Z)
            tmp = root[k-1, 1]
            root[k-1, 1] = root[k, 1]
            root[k, 1] = tmp
    return A, B, Q, Z


def LinApp_Solve(AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, WWW, TT, NN, Z0, Sylv):
    """
    This code takes Uhlig's original code and puts it in the form of a
    function.  This version outputs the policy function coefficients: PP,
    QQ and UU for X, and RR, SS and VV for Y.
    It solves for the decision rules in a linear system,
    which is assumed to be of the form:
    0 = AA x(t) + BB x(t-1) + CC y(t) + DD z(t)
    0 = E_t [ FF x(t+1) + GG x(t) + HH x(t-1) + JJ y(t+1) + KK y(t) + LL z(t+1) + MM z(t)]
    z(t+1) = NN z(t) + epsilon(t+1) with E_t [ epsilon(t+1) ] = 0,
    where it is assumed that x(t) is the endogenous state vector,
    y(t) the other endogenous variables and z(t) the exogenous state
    vector.  It is assumed that the row dimension of AA is at least as large as
    the dimensionality of the endogenous state vector x(t).
    The program solves for the equilibrium law of motion
    x(t) = PP x(t-1) + QQ z(t)
    y(t) = RR x(t-1) + SS z(t).
    To use this function, define the matrices AA, BB, .., TT.
    Inputs overview:
    The matrices of derivatives: AA - TT.
    The autoregression coefficient matrix NN from the law of motion for Z.
    Z0 is the Z-point about which the linearization is taken.  For
    linearizing about the steady state this is Zbar and normally Zbar = 0.
    Sylv is an indicator variable telling the program to use the built-in
    function sylvester() to solve for QQ and SS, if possible.  Default is
    to use Sylv=1.
    Parameters
    ----------
    AA : array_like, dtype=float, shape=(ny, nx)
        The matrix represented above by :math:`A`. It is the matrix of
        derivatives of the Y equations with respect to :math:`X_t`
    BB : array_like, dtype=float, shape=(ny, nx)
        The matrix represented above by :math:`B`. It is the matrix of
        derivatives of the Y equations with respect to
        :math:`X_{t-1}`.
    CC : array_like, dtype=float, shape=(ny, ny)
        The matrix represented above by :math:`C`. It is the matrix of
        derivatives of the Y equations with respect to :math:`Y_t`
    DD : array_like, dtype=float, shape=(ny, nz)
        The matrix represented above by :math:`C`. It is the matrix of
        derivatives of the Y equations with respect to :math:`Z_t`
    FF : array_like, dtype=float, shape=(nx, nx)
        The matrix represented above by :math:`F`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`X_{t+1}`
    GG : array_like, dtype=float, shape=(nx, nx)
        The matrix represented above by :math:`G`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`X_t`
    HH : array_like, dtype=float, shape=(nx, nx)
        The matrix represented above by :math:`H`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`X_{t-1}`
    JJ : array_like, dtype=float, shape=(nx, ny)
        The matrix represented above by :math:`J`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`Y_{t+1}`
    KK : array_like, dtype=float, shape=(nx, ny)
        The matrix represented above by :math:`K`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`Y_t`
    LL : array_like, dtype=float, shape=(nx, nz)
        The matrix represented above by :math:`L`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`Z_{t+1}`
    MM : array_like, dtype=float, shape=(nx, nz)
        The matrix represented above by :math:`M`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`Z_t`
    WWW : array, dtype=float, shape=(ny,)
        The vector of the numerical errors of first ny characterizing
        equations
    TT : array, dtype=float, shape=(nx,)
        The vector of the numerical errors of the next nx characterizing
        equations following the first ny equations
    NN : array_like, dtype=float, shape=(nz, nz)
        The autocorrelation matrix for the exogenous state vector z.
    Z0 : array, dtype=float, shape=(nz,)
        The Z-point about which the linearization is taken.  For linearizing
        about the steady state this is Zbar and normally Zbar = 0.
        QQ if true.
    Sylv: binary, dtype=int
        An indicator variable telling the program to use the built-in
        function sylvester() to solve for QQ and SS, if possible.  Default is
        to use Sylv=1.
    Returns
    -------
    P : 2D-array, dtype=float, shape=(nx, nx)
        The matrix :math:`P` in the law of motion for endogenous state
        variables described above.
    Q : 2D-array, dtype=float, shape=(nx, nz)
        The matrix :math:`Q` in the law of motion for exogenous state
        variables described above.
    U : array, dtype=float, shape=(nx,)
        The vector of the constant term of the policy function for X,
        the endogenous state variables
    R : 2D-array, dtype=float, shape=(ny, nx)
        The matrix :math:`R` in the law of motion for endogenous state
        variables described above.
    S : 2D-array, dtype=float, shape=(ny, nz)
        The matrix :math:`S` in the law of motion for exogenous state
        variables described above.
    V : array, dtype=float, shape=(ny,)
        The vector of the constant term of the policy function for Y,
        the endogenous non-state variables
    References
    ----------
    .. [1] Uhlig, H. (1999): "A toolkit for analyzing nonlinear dynamic
       stochastic models easily," in Computational Methods for the Study
       of Dynamic Economies, ed. by R. Marimon, pp. 30-61. Oxford
       University Press.
    """

    # Make sure that the numpy.array is used

    AA = np.array(AA)
    BB = np.array(BB)
    CC = np.array(CC)
    DD = np.array(DD)
    FF = np.array(FF)
    GG = np.array(GG)
    HH = np.array(HH)
    JJ = np.array(JJ)
    KK = np.array(KK)
    LL = np.array(LL)
    MM = np.array(MM)
    NN = np.array(NN)
    WWW = np.array(WWW)
    TT = np.array(TT)
    Z0 = np.array(Z0)

    # Tolerance level to use

    TOL = .000001

    # Here we use matrices to get pertinent dimensions

    nx = FF.shape[1] # Number of endogenous state variables
    l_equ = CC.shape[0] # Number of deterministic equations
    ny = CC.shape[1] # Number of endogenous control variables
    nz = min(NN.shape) # Number of exogenous variables

    if npla.matrix_rank(CC) < ny:
        print("Error: rank of CC must be at least equal to the number of endogenous control variables, stopping ...")
    else:

        # The following if and else blocks form the
        # Psi, Gamma, Theta Xi, Delta arrays

        if l_equ == 0:
            if CC.any():
                # This block makes sure you don't throw an error with an empty CC.
                CC_plus = la.pinv(CC)
                CC_0 = _nullSpaceBasis(CC.T)
            else:
                CC_plus = np.mat([])
                CC_0 = np.mat([])
            Psi_mat = FF
            Gamma_mat = -GG
            Theta_mat = -HH
            Xi_mat = np.mat(vstack((hstack((Gamma_mat, Theta_mat)),
                                    hstack((eye(nx), zeros((nx, nx)))))))
            Delta_mat = np.mat(vstack((hstack((Psi_mat, zeros((nx, nx)))),
                                       hstack((zeros((nx, nx)), eye(nx))))))

        else:
            CC_plus = la.pinv(CC)
            CC_0 = _nullSpaceBasis(CC.T)
            if l_equ != ny:
                Psi_mat = vstack((zeros((l_equ - ny, nx)), FF
                                  - dot(dot(JJ, CC_plus), AA)))
                Gamma_mat = vstack((dot(CC_0, AA), dot(dot(JJ, CC_plus), BB)
                                    - GG + dot(dot(KK, CC_plus), AA)))
                Theta_mat = vstack((dot(CC_0, BB), dot(dot(KK, CC_plus), BB) - HH))
            else:
                CC_inv = la.inv(CC)
                Psi_mat = FF - dot(JJ.dot(CC_inv), AA)
                Gamma_mat = dot(JJ.dot(CC_inv), BB) - GG + dot(dot(KK, CC_inv), AA)
                Theta_mat = dot(KK.dot(CC_inv), BB) - HH
            Xi_mat = vstack((hstack((Gamma_mat, Theta_mat)),
                             hstack((eye(nx), zeros((nx, nx))))))
            Delta_mat = vstack((hstack((Psi_mat, np.mat(zeros((nx, nx))))),
                                hstack((zeros((nx, nx)), eye(nx)))))

    # Now we need the generalized eigenvalues/vectors for Xi with respect to
    # Delta. That is eVals and eVecs below.

    eVals, eVecs = la.eig(Xi_mat, Delta_mat)

    if npla.matrix_rank(eVecs) < nx:
        print("Error: Xi is not diagonalizable, stopping ...")

    # Diagonalize Xi, form Lambda/Omega, find P

    else:
        Xi_sortabs = np.sort(abs(eVals))
        Xi_sortindex = np.argsort(abs(eVals))
        Xi_sortedVec = np.array([eVecs[:, i] for i in Xi_sortindex]).T
        Xi_sortval = eVals[Xi_sortindex]
        Xi_select = np.arange(0, nx)
        if np.imag(Xi_sortval[nx - 1]).any():
            if (abs(Xi_sortval[nx - 1] - sp.conj(Xi_sortval[nx])) < TOL):
                drop_index = 1
                cond_1 = (abs(np.imag(Xi_sortval[drop_index-1])) > TOL)
                cond_2 = drop_index < nx
                while cond_1 and cond_2:
                    drop_index += 1
                if drop_index >= nx:
                    print("There is an error: too many complex eigenvalues."
                          + "Try increasing the dimension of your state space."
                          + "Stopping ...")
                else:
                    print("Dropping the lowest real eigenvalue to get real PP. Beware of sunspots!")
                    Xi_select = np.array([np.arange(0, drop_index - 1),
                                          np.arange(drop_index, nx + 1)])

        if max(abs(Xi_sortval[Xi_select])) > 1 + TOL:
            print("Potential unstable roots: it might not work ...")
        if abs(max(abs(Xi_sortval[Xi_select])) - 1) < TOL:
            print("Matrix PP contains a unit root. Check the model to make sure of unique steady state." +
                  "There is no convergence back to steady state after a shock" +
                  "Do not trust long simulations.")
        Lambda_mat = np.diag(Xi_sortval[Xi_select])
        Omega_mat = Xi_sortedVec[nx:2 * nx, Xi_select]

        if npla.matrix_rank(Omega_mat) < nx:
            print("Omega matrix is not invertible, Can't solve for P." +
                  "Proceed with the alternative, QZ-method, to get P.")

            # QZ-method codes from SOLVE_QZ #

            # It employs the QZ - method, due to C.A.Sims(1989, 1996)
            # and P. Klein (1997), adapted to the method of undetermined coefficients
            # The QZ-method is perhaps numerically more stable than the generalized
            # eigenvalue method employed in the standard solve algorithm.
            # More importantly, the QZ-method works, even if PP is
            # not diagonalizable.
            # SOLVE_QZ.M in turn makes use of the routines
            # QZDIV.M and QZSWITCH.M, written by C.A. Sims (1996).

            # Theory: with the matrix definitions in the paper,
            # find the QZ-decomposition for Delta and Xi, i.e.
            # unitary matrices U and V and upper triangular matrices
            # Delta_up and Xi_up so that
            # U * Delta * V = Delta_up
            # U * Xi * V = Xi_up
            # and such that the ratios of the diagonal entries
            # | Xi_up(i,i) / Delta_up(i,i) | are ordered in ascend order.
            # Let U[a,b] be the m x m sub-matrices of U, where a,b = 1,2
            # and likewise for V', where V' is the complex conjugate transpose.
            # If V'[2,1] and U[2,1] are invertible
            # then it can be shown that P = - inv(V'[2,1])*V'[2,2]
            # is a solution to the matrix quadratic equation with
            # the most stable roots selected.

            Delta_up, Xi_up, UUU, VVV = la.qz(Delta_mat, Xi_mat, output='complex')
            UUU = UUU.T
            Xi_eigval = np.diag(np.diag(Xi_up)/np.maximum(np.diag(Delta_up), TOL))
            Xi_sortabs = np.sort(abs(np.diag(Xi_eigval)))
            Xi_sortindex = np.argsort(abs(np.diag(Xi_eigval)))
            Xi_sortval = Xi_eigval[Xi_sortindex, Xi_sortindex]
            Xi_select = np.arange(0, nx)
            stake = max(abs(Xi_sortval[Xi_select])) + TOL

            # Sorting the eigenvalues using the Sims (1996) code

            Delta_up, Xi_up, UUU, VVV = qzdiv(stake, Delta_up, Xi_up, UUU, VVV)

            # Check conditions

            if np.imag(Xi_sortval[nx - 1]).any():
                if (abs(Xi_sortval[nx - 1] - sp.conj(Xi_sortval[nx])) < TOL):
                    print("Problem: there are complex eigenvalues!" +
                          "PP matrix will contain complex numbers by this method.")

                drop_index = 1
                cond_1 = (abs(np.imag(Xi_sortval[drop_index-1])) > TOL)
                cond_2 = drop_index < nx
                while cond_1 and cond_2:
                    drop_index += 1
                if drop_index >= nx:
                    print("There is an error: too many complex eigenvalues."
                          + "Try increasing the dimension of your state space."
                          + "Stopping ...")
                else:
                    print("Dropping the lowest real eigenvalue to get real PP. Beware of sunspots!")
                    for i in xrange(drop_index, nx+1):
                        Delta_up, Xi_up, UUU, VVV = qzswitch(i, Delta_up, Xi_up, UUU, VVV)
                    Xi_select1 = np.arange(0, drop_index-1)
                    Xi_select = np.append(Xi_select1, np.arange(drop_index, nx+1))

            if Xi_sortval[max(Xi_select)] < 1 - TOL:
                print("There are stable roots NOT used for PP derivation." +
                      "Proceeding with the smallest root." +
                      "Potential better solution: " +
                      "Move the time index of some endogenous variables back by one " +
                      "and turn them into (predetermined) state variables")
            if max(abs(Xi_sortval[Xi_select])) > 1 + TOL:
                print("Potential unstable roots: it might not work ...")
            if abs(max(abs(Xi_sortval[Xi_select])) - 1) < TOL:
                print("Matrix PP contains a unit root. Check the model to make sure of unique steady state." +
                      "There is no convergence back to steady state after a shock" +
                      "Do not trust long simulations.")

            # End of checking conditions

            # Proceeding with the calculations

            VVV = VVV.conj().T
            VVV_2_1 = VVV[nx: 2*nx, 0: nx]
            VVV_2_2 = VVV[nx: 2*nx, nx:2*nx]
            UUU_2_1 = UUU[nx: 2*nx, 0: nx]
            VVV = VVV.conj().T

            if abs(la.det(UUU_2_1)) < TOL:
                print("One necessary condition for computing P is NOT satisfied")
            if abs(la.det(VVV_2_1)) < TOL:
                print("VVV_2_1 matrix, used to compute for P, is not invertible")

            PP = np.array(la.solve(- VVV_2_1, VVV_2_2))
            PP_imag = np.imag(PP)
            PP = np.real(PP)
            if (sum(sum(abs(PP_imag))) / sum(sum(abs(PP))) > .000001).any():
                print("PP is complex." +
                      "Continue with the real part, hoping not to lose too much information.")

            # End of QZ-method #

        # This follows the original uhlig.py file

        else:
            PP = dot(dot(Omega_mat, Lambda_mat), la.inv(Omega_mat))
            PP_imag = np.imag(PP)
            PP = np.real(PP)
            if (sum(sum(abs(PP_imag))) / sum(sum(abs(PP))) > .000001).any():
                print("PP is complex." +
                      "Continue with the real part, hoping not to lose too much information.")

    # The code from here to the end was from the Uhlig file calc_qrs.m

    # The if and else below make RR and VV depending on the model's setup

    if l_equ == 0:
        RR = zeros((0, nx))
        VV = hstack((kron(NN.T, FF) + kron(eye(nz),
                                           (dot(FF, PP) + GG)), kron(NN.T, JJ) + kron(eye(nz), KK)))

    else:
        RR = - dot(CC_plus, (dot(AA, PP) + BB))
        VV = sp.vstack((hstack((kron(eye(nz), AA),
                                kron(eye(nz), CC))), hstack((kron(NN.T, FF) +
                                                             kron(eye(nz), dot(FF, PP) + dot(JJ, RR) + GG),
                                                             kron(NN.T, JJ) + kron(eye(nz), KK)))))

    # Use LL, NN, RR, VV to get the QQ, RR, SS, VV matrices

    # First try using Sylvester equation solver

    if ny > 0:
        PM = (FF - JJ.dot(la.solve(CC, AA)))
        if npla.matrix_rank(PM) < nx+ny:
            Sylv = 0
            print("Sylvester equation solver condition is not satisfied." +
                  "Proceed with the original method.")
    else:
        if npla.matrix_rank(FF) < nx:
            Sylv = 0
            rint("Sylvester equation solver condition is not satisfied." +
                 "Proceed with the original method.")

    if Sylv:
        print("Using Sylvester equation solver.")
        if ny > 0:
            Anew = la.solve(PM, (FF.dot(PP)+GG+JJ.dot(RR) -
                                 la.solve(KK.dot(CC), AA)))
            Bnew = NN
            Cnew1 = la.solve(JJ.dot(CC), DD.dot(NN))+la.solve(KK.dot(CC), DD) -\
                LL.dot(NN)-MM
            Cnew = la.solve(PM, Cnew1)
            QQ = la.solve_sylvester(Anew, Bnew, Cnew)
            SS = la.solve(-CC, (AA.dot(QQ)+DD))
        else:
            Anew = la.solve(FF, (FF.dot(PP)+GG))
            Bnew = NN
            Cnew = la.solve(FF, (-LL.dot(NN)-MM))
            QQ = la.solve_sylvester(Anew, Bnew, Cnew)
            SS = np.zeros((0, nz))

    # Otherwise the Uhlig's way

    else:
        if (npla.matrix_rank(VV) < nz * (nx + ny)):
            print("V is not invertible. Can't solve for Q and S.")

        LL = sp.mat(LL)
        NN = sp.mat(NN)
        LLNN_plus_MM = dot(LL, NN) + MM

        if DD.any():
            DD = np.array(DD).T.ravel().reshape((-1, 1))
            LLNN_plus_MM = np.array(LLNN_plus_MM).reshape((-1, 1))
            impvec = np.vstack((DD, LLNN_plus_MM))
        else:
            impvec = LLNN_plus_MM.flatten()

        QQSS_vec = np.array(la.solve(-VV, impvec))

        if (max(abs(QQSS_vec)) == sp.inf).any():
            print("There are issues with Q and S. Entries are undefined." +
                  " Probably because V is not invertible.")

        # Build QQ and SS

        QQ = np.reshape(np.array(QQSS_vec[0:nx * nz, 0]),
                        (nx, nz), 'F')

        SS = np.reshape(QQSS_vec[(nx * nz):((nx + ny) * nz), 0],
                        (ny, nz), 'F')

    # Build WW - WW has the property [x(t)',y(t)',z(t)']=WW [x(t)',z(t)']

    WW = sp.vstack((
        hstack((eye(nx), zeros((nx, nz)))),
        hstack((dot(RR, la.pinv(PP)), (SS - dot(dot(RR, la.pinv(PP)), QQ)))),
        hstack((zeros((nz, nx)), eye(nz)))))

    # Find constant terms
    # Redefine matrices to be 2D-arrays for generating vector UU and VVV

    AA = np.array(AA)
    CC = np.array(CC)
    FF = np.array(FF)
    GG = np.array(GG)
    JJ = np.array(JJ)
    KK = np.array(KK)
    LL = np.array(LL)
    NN = np.array(NN)
    RR = np.array(RR)
    QQ = np.array(QQ)
    SS = np.array(SS)

    if ny > 0:
        UU1 = -(FF.dot(PP)+GG+JJ.dot(RR)+FF-(JJ+KK).dot(la.solve(CC, AA)))
        UU2 = (TT+(FF.dot(QQ)+JJ.dot(SS)+LL).dot(NN.dot(Z0)-Z0) -
               (JJ+KK).dot(la.solve(CC, WWW)))
        UU = la.solve(UU1, UU2)
        VVV = la.solve(- CC, (WWW+AA.dot(UU)))
    else:
        UU = la.solve(-(FF.dot(PP)+FF+GG), (TT+(FF.dot(QQ)+LL).dot(NN.dot(Z0)-Z0)))
        VVV = np.array([])

    return np.array(PP), np.array(QQ), np.array(UU), np.array(RR), np.array(SS),\
        np.array(VVV)
