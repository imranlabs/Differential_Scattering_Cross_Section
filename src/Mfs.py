#foldy_lax
import numpy as np
from scipy.optimize import curve_fit

# function to compute Green's function

def ComputeG( k, Rd ):
    """
    This function computes the whole space Green's function given a wavenumber k and a distance Rd.
    
    The output of this function is the evaluation of the whole space Green's function.
    """
    
    # compute Green's function
    
    G = np.exp( 1j * k * Rd ) / ( 4 * np.pi * Rd )
      
    return G;
    
    # function to compute the normal derivative of Green's function

def ComputeDνG( k, Rd, CosTheta ):
    """
    This function computes the normal derivative of the whole space Green's function 
    given a wavenumber k, a distance Rd, and the cosine of the angle made between the
    difference vector and the unit normal.
    
    The output of this function is the evaluation of the normal derivative of the whole 
    space Green's function.
    """
    
    # compute Green's function
    
    G = ComputeG( k, Rd )
    
    # compute the normal derivative of Green's function
    
    DνG = CosTheta * ( 1j * k - 1 / Rd ) * G
    
    return DνG

#Compute the MFS expansion coefficients and Foldy-Lax exciting fields.

def ComputeExpansionCoefficientsExcitingFields( k0, k1, a0, ν, ρ_bdy, ρ_int, ρ_ext, α, R_N, M, N ):
    
    
    # interior points

    indx, jndx = np.meshgrid( np.arange(0,M), np.arange(0,M) )

    R_int = np.sqrt( ( ( ρ_bdy[indx] - ρ_int[jndx] ) ** 2 ).sum( axis = 2 ) )
    μ_int = np.divide(( ν[indx] * ( ρ_bdy[indx] - ρ_int[jndx] ) ).sum( axis = 2 ), R_int)

    # exterior points

    R_ext = np.sqrt( ( ( ρ_bdy[indx] - ρ_ext[jndx] ) ** 2 ).sum( axis = 2 ) )
    μ_ext = np.divide(( ν[indx] * ( ρ_bdy[indx] - ρ_ext[jndx] ) ).sum( axis = 2 ) , R_ext)
    
    # Foldy-Lax points
    
    indx, jndx = np.meshgrid( np.arange(0,M), np.arange(0,N) )

    ρ_bdy2RN = np.sqrt( ( ( ρ_bdy[indx] - R_N[jndx] ) ** 2 ).sum( axis = 2 ) )
    μ_RN     = np.divide(( ν[indx] * ( ρ_bdy[indx] - R_N[jndx] ) ).sum( axis = 2 ), ρ_bdy2RN)
    RN2ρ_ext = np.sqrt( ( ( ρ_ext[indx] - R_N[jndx] ) ** 2 ).sum( axis = 2 ) )

    indx, jndx = np.meshgrid( np.arange(0,N), np.arange(0,N) )

    RN2RN    = np.sqrt( ( ( R_N[indx] - R_N[jndx] ) ** 2 ).sum( axis = 2 ) )
    
    # compute the incident field and its normal derivative on the M boundary points
    
    Ψ_inc   = np.exp( 1j * k0 * ρ_bdy[:,2] )
    DνΨ_inc = 1j * k0 * np.multiply( ν[:,2], Ψ_inc )

    # compute the incident field on the N point scatterers
    
    Ψ_inc4FL = np.exp( 1j * k0 * R_N[:,2] )
    
    # compute the matrix blocks
    
    A11 =  ComputeG( k1, R_int )
    A12 = - ComputeG( k0, R_ext )
    A13 = -α * ComputeG( k0, ρ_bdy2RN ).T
    
    A21 =  ComputeDνG( k1, R_int, μ_int )
    A22 = - ComputeDνG( k0, R_ext, μ_ext )
    A23 = - α * ComputeDνG( k0, ρ_bdy2RN, μ_RN ).T
    
    A31 = np.zeros( ( N, M ), dtype ='complex' )
    A32 = - ComputeG( k0, RN2ρ_ext )
    
    A33 = np.zeros( ( N, N ), dtype = 'complex')
    offdiags = np.where( RN2RN != 0 )
    A33[offdiags] = -α * ( np.exp( 1j * k0 * RN2RN[offdiags]  ) / ( 4 * np.pi * RN2RN[offdiags] ) )
    #A33 = -α * np.divide( np.exp( 1j * k0 * RN2RN ), ( 4 * np.pi * RN2RN ) )
    A33 += np.eye( N, dtype = 'complex' )
   
    # form the linear system of equations
    
    A = np.block( [ [ A11, A12, A13 ], [ A21, A22, A23 ], [ A31, A32, A33 ] ] )
    b = np.block( [ Ψ_inc, DνΨ_inc, Ψ_inc4FL ] )
    
    # solve the linear system
    
    c = np.linalg.solve( A, b )
    
    # parse the solution
    
    c_int = c[0:M]
    c_ext = c[M:2*M]
    Ψ_E   = c[2*M:2*M*N]
    
    return c_int, c_ext, Ψ_E



def ComputeMFSPoints( a0, M ):
    """
    This function computes the set of points needed for the method of fundamental solutions (MFS).
    In particular, given the radius of a sphere, a0, and the number of points M, this function computes
    the Fibonnaci lattice on the unit sphere and stores them as the unit normal vectors ν. Using ν, we
    then compute ρ_bdy = a0 * ν, ρ_int = ( a0 + ℓ ) * ν, and ρ_ext = ( a0 - ℓ ) * ν.
    
    This function outputs four vectors: ν, ρ_bdy, ρ_int, and ρ_ext.
    """
    
    # allocate memory for the Fibonacci lattice points on the unit sphere

    ν = np.full( ( M, 3 ), float( 'nan' ) )

    # compute the "golden angle"

    golden_angle = np.pi * ( 3 - np.sqrt( 5 ) )

    # compute the points on the unit sphere

    ν[:,2] = ( 1 - 1 / M ) * ( 1 - 2 * np.arange( 0, M ) / ( M - 1 ) )

    ρ = np.sqrt( 1 - ν[:,2] ** 2 )
    θ = golden_angle * np.arange( 0, M )

    ν[:,0] = ρ * np.cos( θ )
    ν[:,1] = ρ * np.sin( θ )

    # compute the boundary points, interior points, and exterior points

    ℓ = 0.25 * a0

    ρ_bdy = a0 * ν
    ρ_int = ( a0 + ℓ ) * ν
    ρ_ext = ( a0 - ℓ ) * ν
    
    return ν, ρ_bdy, ρ_int, ρ_ext


# B. Compute the MFS expansion coefficients 

def ComputeMFSExpansionCoefficients( k0, k1, a0, ν, ρ_bdy, ρ_int, ρ_ext, M ):
    """
    This function solves the 2M x 2M system of equations for the MFS expansion coefficients.
    
    This code requires the results from ComputeMFSPoints, namely ν, ρ_bdy, ρ_int, and ρ_ext, in 
    addition to the two wavenumbers k0 and k1, the sphere radius, a0, and the number of MFS points, 
    M.
    
    The output from this code are the 2 M-vectors, c_int and c_sca, corresponding to the MFS
    expansions for the interior and scattered fields, respectively.
    """
    
    # interior points

    indx, jndx = np.meshgrid( np.arange(0,M), np.arange(0,M) )

    R_int = np.sqrt( ( ( ρ_bdy[indx] - ρ_int[jndx] ) ** 2 ).sum( axis = 2 ) )
    μ_int = ( ν[indx] * ( ρ_bdy[indx] - ρ_int[jndx] ) ).sum( axis = 2 ) / R_int

    # exterior points

    R_ext = np.sqrt( ( ( ρ_bdy[indx] - ρ_ext[jndx] ) ** 2 ).sum( axis = 2 ) )
    μ_ext = ( ν[indx] * ( ρ_bdy[indx] - ρ_ext[jndx] ) ).sum( axis = 2 ) / R_ext
    
    # compute the incident field and its normal derivative on the M boundary points
    
    Ψ_inc   = np.exp( 1j * k0 * ρ_bdy[:,2] )
    DνΨ_inc = 1j * k0 * np.multiply( ν[:,2], Ψ_inc )
    
    # compute the matrix blocks
    
    A11 =  ComputeG( k1, R_int )
    A12 = -ComputeG( k0, R_ext )
    
    A21 =  ComputeDνG( k1, R_int, μ_int )
    A22 = -ComputeDνG( k0, R_ext, μ_ext )
   
    # form the linear system of equations
    
    A = np.block( [ [ A11, A12 ], [ A21, A22 ] ] )
    b = np.block( [ Ψ_inc, DνΨ_inc ] )
    
    # solve the linear system
    
    c = np.linalg.solve( A, b )
    
    # parse the solution
    
    c_int0 = c[0:M]
    c_ext0 = c[M:2*M]
    
    return c_int0, c_ext0

# Computing the total cross section


def ComputeTotalCrossSection( k0, k1, a0, ν, ρ_bdy, ρ_int, ρ_ext, α, R_N, M, N ):
    """
    This function computes the total cross-section by evaluating the Optical Theorem using
    the results from the MFS approximation for the scattered field.
    """
    
    # compute the MFS expansion coefficients
    
    c_int, c_ext, Ψ_E = ComputeExpansionCoefficientsExcitingFields( k0, k1, a0, ν, ρ_bdy, ρ_int, ρ_ext, α, R_N, M, N )
    
    # compute the scattering amplitude
    
    f1 = (0.25 / np.pi)  * np.exp( -1j * k0 * ρ_ext[:,2] ).T @ c_ext
    f2 = (0.25 / np.pi ) * α * np.exp( -1j * k0 * R_N[:,2] ).T @ Ψ_E
    
    # compute the scattering cross-section
    
    σ_t = 4 * np.pi * np.imag( f1 + f2 ) / k0
    
    return σ_t

def ComputeTotalCrossSectionMFS( k0, k1, a0, ν, ρ_bdy, ρ_int, ρ_ext, M ):
    """
    This function computes the total cross-section by evaluating the Optical Theorem using
    the results from the MFS approximation for the scattered field.
    """
    
    # compute the MFS expansion coefficients
    
    c_int, c_sca = ComputeMFSExpansionCoefficients( k0, k1, a0, ν, ρ_bdy, ρ_int, ρ_ext, M )
    
    # compute the scattering amplitude
    
    f = 0.25 / np.pi * np.exp( -1j * k0 * ρ_ext[:,2] ).T @ c_sca
    
    # compute the scattering cross-section
    
    σ_t = 4 * np.pi * np.imag( f ) / k0
    
    return σ_t
    
    
# Compute Gauss-Legendre points for integration

def GaussLegendre(N):
    """
    Compute the Gauss Legendre quadrature points and weights using
    the method given in Spectral Methods in MATLAB by L. N. Trefethen (2000).

    Direct implementation by A. D. Kim.
    Converted to Python.

    Parameters
    ----------
    N : int
        Number of quadrature points.

    Returns
    -------
    x : numpy.ndarray
        Quadrature points (roots).
    w : numpy.ndarray
        Quadrature weights.
    """
    # beta = 0.5 / sqrt(1 - (2*(1:N-1)).^(-2))
    # Use np.arange(1, N) to get [1, 2, ..., N-1]
    indices = np.arange(1, N)
    beta = 0.5 / np.sqrt(1.0 - (2.0 * indices) ** (-2))

    # Build the symmetric tridiagonal matrix T
    # Use np.diag to create a matrix with 'beta' on the +1 and -1 diagonals
    T = np.diag(beta, k=1) + np.diag(beta, k=-1)

    # Compute the eigenvalues and eigenvectors of T
    D, V = np.linalg.eig(T)  # Note: np.linalg.eig returns eigenvalues first, then eigenvectors in columns of V.

    # Sort eigenvalues (quadrature points) and eigenvectors
    # np.linalg.eig returns eigenvalues in no particular order, so we sort them.
    i = np.argsort(D)  # Get the indices that would sort the array of eigenvalues
    x = D[i]           # Sorted eigenvalues are the quadrature points 'x'
    w = 2 * (V[0, i] ** 2)  # The weights are 2 times the square of the first component of each eigenvector

    return x, w
    
    
# Define the Henyey-Greenstein scattering law function
def HG_function(mu, albedo, g):
    """
    Henyey-Greenstein scattering phase function
    
    Parameters
    ----------
    mu : array_like
        Cosine of scattering angle (cosθ)
    albedo : float
        Single scattering albedo
    g : float
        Anisotropy factor
    
    Returns
    -------
    array_like
        HG scattering law values
    """
    return (albedo / (4 * np.pi)) * (1 - g**2) / (1 + g**2 - 2 * g * mu)**(3/2)
