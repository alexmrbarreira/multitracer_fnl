from numpy import *
from scipy import *
from pylab import *
import pylab
rcParams.update({'text.usetex': True, 'mathtext.fontset': 'stix'}) 
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

## ======================================================================
## General parameters
## ======================================================================

lpath = '/home/barreira/a.Other/d-github_workrepos/d-multitracer_fnl/'

Gnewton = 4.302e-9 #Mpc (km/s)^2 / Msun
c_light = 299792.458 #km/s

# Cosmology
fNL    = 5.
h      = 0.6774 # 0.6774 for consistency with lookup tables
H0     = 100. * h #km/s/Mpc
Om0    = 0.3089 # 0.3089 for consistency with lookup tables
Ob0    = 0.022301   /h**2. # 0.022301/0.6774**2 for consistency with lookup tables
rho_m0 = Om0 * 3.* (100.)**2. / (8.*pi*Gnewton) #h^2 Msun/Mpc^3

# Multiply matter transfer function data from CAMB by numfac_matter to get transfer function in (km/s)^2/h^2 and with definition delta_m(k,z) = 2/(5 Omega_m H0^2) * k^2 * R(k) * T_m(k)
numfac_matter = h**2./(2./5./Om0/H0**2.) 

# Critical sph. collapse density
dc = 1.686

# Survey setup
z    = 1.0
Vs   = 100. * 1.0e9 # Mpc^3 / h^3 
kF   = pi/Vs**(1./3)
kmax = 0.2 #h/Mpc
Nk   = 40
kk_edges = exp(linspace(log(kF),log(kmax),Nk+1,endpoint=True))
kk  = (kk_edges[1::] + kk_edges[0:-1])/2.
DDk =  kk_edges[1::] - kk_edges[0:-1]

print ('')
print ('Fiducial and survey spec. values:')
print ('fNL    = ', fNL)
print ('z      = ', z)
print ('Vs     = ', Vs*1.0e-9, '[Gpc^3/h^3]')
print ('kF     = ', kF, 'h/Mpc')
print ('kmax   = ', kmax, 'h/Mpc')
print ('Nkbins = ', len(kk))

# ================================================================================ #
# Load linear matter power spectra and transfer functions; create 2D interpolators
# Needs data to be produced in lookup_tables/ first
# ================================================================================ #
# Create the interpolators that are needed; it does trick of setting to zero (or linear result for responses) for k = 0 (to prevent crashes/problems when k<k_min)
# This needs to be adapted if the format of the files changes
def make_interpolator_forP(path1, path2, value_at_k0):
    filetmp  = loadtxt(path1, skiprows = 1)
    ztmp     = loadtxt(path2)
    ktmp     = append(0., filetmp[:,0])
    Ptmp     = vstack([value_at_k0*ones(len(ztmp)), filetmp[:,1:]])
    return interpolate.interp2d(ktmp, ztmp, transpose(Ptmp), kind='linear')

def make_interpolator_forT(path1, path2, value_at_k0, numfactor):
    filetmp  = loadtxt(path1, skiprows = 1)
    ztmp     = loadtxt(path2)
    ktmp     = append(0., filetmp[:,0])
    Ptmp     = vstack([value_at_k0*ones(len(ztmp)), filetmp[:,1:]])
    return interpolate.interp2d(ktmp, ztmp, transpose(Ptmp*numfactor), kind='linear')

# P_L(k, z) w/ Camb
Plin_int  = make_interpolator_forP(lpath+'lookup_tables/P_L_camb.dat'  , lpath+'lookup_tables/P_L_camb_zvalues.dat', 0.0)
# T_m(k, z) w/ Camb
T_m_int   = make_interpolator_forT(lpath+'lookup_tables/T_m_camb.dat'  , lpath+'lookup_tables/P_L_camb_zvalues.dat', 0.0, numfac_matter)
# M(k, z)
def M(k,z):
    return (2./3) * k**2. * T_m_int(k,z) / Om0 / H0**2.

# ==================================================================== #
# Define model, model derivatives and covariance
# ==================================================================== #

# For single tracer
def model_singletracer(z, k, fNL, b1, bphi, ng):
    Pk = Plin_int(k,z)[0]
    Mk = M(k,z)[0]
    return (b1 + bphi*fNL/Mk)**2. * Pk + 1./ng

def model_derivative_singletracer(z, k, fNL, b1, bphi):
    Pk = Plin_int(k,z)[0]
    Mk = M(k,z)[0]
    return 2.*(b1 + bphi*fNL/Mk) * (bphi/Mk) * Pk

def model_derivative_singletracer_fNLbphi(z, k, fNL, b1, bphi):
    Pk = Plin_int(k,z)[0]
    Mk = M(k,z)[0]
    return 2.*(b1 + bphi*fNL/Mk) * (1./Mk) * Pk

def covariance_singletracer(z, k, fNL, b1, bphi, ng, Vs, Dk):
    Vk_factor   = (4.*pi) * k**2. * Dk
    mult_factor = 2.*(2.*pi)**3. / Vs / Vk_factor
    return mult_factor * model_singletracer(z, k, fNL, b1, bphi, ng)**2.

# For multi tracer
def model_multitracer(z, k, fNL, b1_A, b1_B, bphi_A, bphi_B, ng_A, ng_B):
    Pk   = Plin_int(k,z)[0]
    Mk   = M(k,z)[0]
    P_AA = (b1_A + bphi_A*fNL/Mk)**2.                              * Pk + 1./ng_A
    P_AB = (b1_A + bphi_A*fNL/Mk)     * (b1_B + bphi_B*fNL/Mk)     * Pk
    P_BB =                              (b1_B + bphi_B*fNL/Mk)**2. * Pk + 1./ng_B
    return array([P_AA, P_AB, P_BB])

def model_derivative_multitracer(z, k, fNL, b1_A, b1_B, bphi_A, bphi_B):
    Pk   = Plin_int(k,z)[0]
    Mk   = M(k,z)[0]
    d_AA = 2.*(b1_A + bphi_A*fNL/Mk) * (bphi_A/Mk) * Pk 
    d_AB = ( (b1_A + bphi_A*fNL/Mk)*(bphi_B/Mk) + (b1_B + bphi_B*fNL/Mk)*(bphi_A/Mk) ) * Pk
    d_BB = 2.*(b1_B + bphi_B*fNL/Mk) * (bphi_B/Mk) * Pk
    return array([d_AA, d_AB, d_BB])

def model_derivative_multitracer_wrt_fNLbphiA(z, k, fNL, b1_A, b1_B, bphi_A, bphi_B):
    Pk   = Plin_int(k,z)[0]
    Mk   = M(k,z)[0]
    d_AA = 2.*(b1_A + bphi_A*fNL/Mk) * (1./Mk) * Pk
    d_AB =    (b1_B + bphi_B*fNL/Mk) * (1./Mk) * Pk
    d_BB =                                  0. * Pk
    return array([d_AA, d_AB, d_BB])

def model_derivative_multitracer_wrt_fNLbphiB(z, k, fNL, b1_A, b1_B, bphi_A, bphi_B):
    Pk   = Plin_int(k,z)[0]
    Mk   = M(k,z)[0]
    d_AA =                                  0. * Pk
    d_AB =    (b1_A + bphi_A*fNL/Mk) * (1./Mk) * Pk
    d_BB = 2.*(b1_B + bphi_B*fNL/Mk) * (1./Mk) * Pk
    return array([d_AA, d_AB, d_BB])

def covariance_multitracer(z, k, fNL, b1_A, b1_B, bphi_A, bphi_B, ng_A, ng_B, Vs, Dk):
    Pk   = Plin_int(k,z)[0]
    Mk   = M(k,z)[0]
    P_AA = (b1_A + bphi_A*fNL/Mk)**2.                              * Pk + 1./ng_A
    P_AB = (b1_A + bphi_A*fNL/Mk)     * (b1_B + bphi_B*fNL/Mk)     * Pk
    P_BB =                              (b1_B + bphi_B*fNL/Mk)**2. * Pk + 1./ng_B
    Vk_factor   = (4.*pi) * k**2. * Dk
    mult_factor = 2.*(2.*pi)**3. / Vs / Vk_factor
    cov_matrix  = zeros([3,3])
    cov_matrix[0,0] = mult_factor * P_AA**2.
    cov_matrix[1,0] = mult_factor * P_AB * P_AA
    cov_matrix[2,0] = mult_factor * P_AB**2.
    cov_matrix[0,1] = cov_matrix[1,0]
    cov_matrix[1,1] = mult_factor * (P_AA*P_BB + P_AB**2.)/2.
    cov_matrix[2,1] = mult_factor * P_BB*P_AB
    cov_matrix[0,2] = cov_matrix[2,0]
    cov_matrix[1,2] = cov_matrix[2,1]
    cov_matrix[2,2] = mult_factor * P_BB**2.
    return cov_matrix 


