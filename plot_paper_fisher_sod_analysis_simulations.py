from commons import *

## =========================================================================
## Some functions and parameters 
## ========================================================================= 

# List of redshifts that can be selected
list_redshift_strings   = ['z00', 'z05', 'z10', 'z20', 'z30']
list_redshift_values    = [0.0, 0.5, 1.0, 2.0, 3.0]
annotate_list_redshifts = [r"$z = 0$", r"$z = 0.5$", r"$z = 1$", r"$z = 2$", r"$z = 3$"]

# List of number densities that can be selected
ng_targets = [5.0e-4, 2.0e-4, 1.0e-4] # h^3/Mpc^3
ng_strings = ['nbar_5.0e-4', 'nbar_2.0e-4', 'nbar_1.0e-4']
ng_labels  = [r'$\bar{n}_g = 5\times 10^{-4}$', r'$\bar{n}_g = 2\times 10^{-4}$', r'$\bar{n}_g = 1\times 10^{-4}$']

# Select here the desired redshift and number density (index of the above lists)
z_index  = 2
ng_index = 1

pathtosimdata = 'data_simulations_bias/' 

# Function that loads bias parameters for galaxies in IllustrisTNG
def load_bphib1_v1(zind, nind, split_type):
    b1_val_high = load(pathtosimdata+'data_b1/'+list_redshift_strings[zind]+'/data_b1_val_tng300_2_hydro_totmass_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_high_'+split_type+'.npy')
    b1_err_high = load(pathtosimdata+'data_b1/'+list_redshift_strings[zind]+'/data_b1_err_tng300_2_hydro_totmass_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_high_'+split_type+'.npy')
    b1_val_loww = load(pathtosimdata+'data_b1/'+list_redshift_strings[zind]+'/data_b1_val_tng300_2_hydro_totmass_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_loww_'+split_type+'.npy')
    b1_err_loww = load(pathtosimdata+'data_b1/'+list_redshift_strings[zind]+'/data_b1_err_tng300_2_hydro_totmass_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_loww_'+split_type+'.npy')
    bphi_val_high = load(pathtosimdata+'data_bphi/'+list_redshift_strings[zind]+'/data_bphi_val_tng300_2_hydro_totmass_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_high_'+split_type+'.npy')
    bphi_err_high = load(pathtosimdata+'data_bphi/'+list_redshift_strings[zind]+'/data_bphi_err_tng300_2_hydro_totmass_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_high_'+split_type+'.npy')
    bphi_val_loww = load(pathtosimdata+'data_bphi/'+list_redshift_strings[zind]+'/data_bphi_val_tng300_2_hydro_totmass_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_loww_'+split_type+'.npy')
    bphi_err_loww = load(pathtosimdata+'data_bphi/'+list_redshift_strings[zind]+'/data_bphi_err_tng300_2_hydro_totmass_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_loww_'+split_type+'.npy')
    return b1_val_high, b1_err_high, b1_val_loww, b1_err_loww, bphi_val_high, bphi_err_high, bphi_val_loww, bphi_err_loww

# Function that loads bias parameters for halos in a gravity-only simulation
def load_bphib1_v2(zind, nind, split_type):
    b1_val_high = load(pathtosimdata+'data_b1/'+list_redshift_strings[zind]+'/data_b1_val_tng800_2_dmo_M200_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_high_'+split_type+'.npy')
    b1_err_high = load(pathtosimdata+'data_b1/'+list_redshift_strings[zind]+'/data_b1_err_tng800_2_dmo_M200_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_high_'+split_type+'.npy')
    b1_val_loww = load(pathtosimdata+'data_b1/'+list_redshift_strings[zind]+'/data_b1_val_tng800_2_dmo_M200_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_loww_'+split_type+'.npy')
    b1_err_loww = load(pathtosimdata+'data_b1/'+list_redshift_strings[zind]+'/data_b1_err_tng800_2_dmo_M200_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_loww_'+split_type+'.npy')
    bphi_val_high = load(pathtosimdata+'data_bphi/'+list_redshift_strings[zind]+'/data_bphi_val_tng800_2_dmo_M200_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_high_'+split_type+'.npy')
    bphi_err_high = load(pathtosimdata+'data_bphi/'+list_redshift_strings[zind]+'/data_bphi_err_tng800_2_dmo_M200_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_high_'+split_type+'.npy')
    bphi_val_loww = load(pathtosimdata+'data_bphi/'+list_redshift_strings[zind]+'/data_bphi_val_tng800_2_dmo_M200_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_loww_'+split_type+'.npy')
    bphi_err_loww = load(pathtosimdata+'data_bphi/'+list_redshift_strings[zind]+'/data_bphi_err_tng800_2_dmo_M200_' + list_redshift_strings[zind] + '_' + ng_strings[nind] + '_loww_'+split_type+'.npy')
    return b1_val_high, b1_err_high, b1_val_loww, b1_err_loww, bphi_val_high, bphi_err_high, bphi_val_loww, bphi_err_loww

# Function that computes the significance-of-detection for single-tracer
def compute_sod_st(b1, bphi, ng):
    Finfo = zeros(len(kk))
    for i in range(len(kk)):
        Finfo[i] = (model_derivative_singletracer_fNLbphi(z, kk[i], fNL, b1, bphi)**2. / covariance_singletracer(z, kk[i], fNL, b1, bphi, ng, Vs, DDk[i]))
    sigma_fNLbphi = 1./sqrt(sum(Finfo))
    sod_fNLbphi   = fNL*bphi/sigma_fNLbphi
    return sod_fNLbphi

# Function that computes the significance-of-detection for multi-tracer
def compute_sod_mt(b1_A, b1_B, bphi_A, bphi_B, ng_A, ng_B):
    Finfo_AA = zeros(len(kk))
    Finfo_AB = zeros(len(kk))
    Finfo_BB = zeros(len(kk))
    for i in range(len(kk)):
        Finfo_AA[i] = dot(model_derivative_multitracer_wrt_fNLbphiA(z, kk[i], fNL, b1_A, b1_B, bphi_A, bphi_B),
                      dot(linalg.inv(covariance_multitracer(z, kk[i], fNL, b1_A, b1_B, bphi_A, bphi_B, ng_A, ng_B, Vs, DDk[i])),
                          model_derivative_multitracer_wrt_fNLbphiA(z, kk[i], fNL, b1_A, b1_B, bphi_A, bphi_B)))

        Finfo_AB[i] = dot(model_derivative_multitracer_wrt_fNLbphiA(z, kk[i], fNL, b1_A, b1_B, bphi_A, bphi_B),
                      dot(linalg.inv(covariance_multitracer(z, kk[i], fNL, b1_A, b1_B, bphi_A, bphi_B, ng_A, ng_B, Vs, DDk[i])),
                          model_derivative_multitracer_wrt_fNLbphiB(z, kk[i], fNL, b1_A, b1_B, bphi_A, bphi_B)))

        Finfo_BB[i] = dot(model_derivative_multitracer_wrt_fNLbphiB(z, kk[i], fNL, b1_A, b1_B, bphi_A, bphi_B),
                      dot(linalg.inv(covariance_multitracer(z, kk[i], fNL, b1_A, b1_B, bphi_A, bphi_B, ng_A, ng_B, Vs, DDk[i])),
                          model_derivative_multitracer_wrt_fNLbphiB(z, kk[i], fNL, b1_A, b1_B, bphi_A, bphi_B)))

    F_matrix = array([[sum(Finfo_AA), sum(Finfo_AB)],[sum(Finfo_AB), sum(Finfo_BB)]])
    Mean     = array([fNL*bphi_A, fNL*bphi_B])
    invCov   = F_matrix
    sod      = sqrt(dot(Mean, dot(invCov, Mean)))
    return sod

def print_bias_st(b1, bphi, split_type):
    print ('')
    print ('Bias parameter for split by', split_type)
    print ('b1   = ', b1)
    print ('bphi = ', bphi)
    print ('fNLbphi = ', fNL*bphi)
    return 0

def print_bias_mt(b1_A, b1_B, bphi_A, bphi_B, split_type):
    print ('')
    print ('Bias parameter for split by', split_type)
    print ('b1_A    = ', b1_A)
    print ('b1_B    = ', b1_B)
    print ('bphi_A  = ', bphi_A)
    print ('bphi_B  = ', bphi_B)
    print ('fNLbphi_A = ', fNL*bphi_A)
    print ('fNLbphi_B = ', fNL*bphi_B)
    return 0

## =========================================================================
## Load bias parameters  and compute their sigma_fNL
## ========================================================================= 

print ('')
print ('Results for')
print ('z = ', annotate_list_redshifts[z_index])
print ('ng = ', ng_labels[ng_index] + r'$\left[h^3/{\rm Mpc^3}\right]$')
print ('')

ng_st = ng_targets[ng_index]
ng_A   = ng_st/2.
ng_B   = ng_st/2.

bb1 = linspace(0., 8., 100)

xmin = 1.0
xmax = 3.9
ymin = -4.
ymax = 15.

# Full sample of galaxies (same so use any)
b1_galaxy_val   = load(pathtosimdata+'data_b1/'+list_redshift_strings[z_index]+'/data_b1_val_tng300_2_hydro_totmass_'   + list_redshift_strings[z_index] + '_' + ng_strings[ng_index] + '.npy')
b1_galaxy_err   = load(pathtosimdata+'data_b1/'+list_redshift_strings[z_index]+'/data_b1_err_tng300_2_hydro_totmass_'   + list_redshift_strings[z_index] + '_' + ng_strings[ng_index] + '.npy')
bphi_galaxy_val = load(pathtosimdata+'data_bphi/'+list_redshift_strings[z_index]+'/data_bphi_val_tng300_2_hydro_totmass_' + list_redshift_strings[z_index] + '_' + ng_strings[ng_index] + '.npy')
bphi_galaxy_err = load(pathtosimdata+'data_bphi/'+list_redshift_strings[z_index]+'/data_bphi_err_tng300_2_hydro_totmass_' + list_redshift_strings[z_index] + '_' + ng_strings[ng_index] + '.npy')
print_bias_st(b1_galaxy_val, bphi_galaxy_val, 'full galaxy sample')
sigma_fNL_galaxy = compute_sod_st(b1_galaxy_val, bphi_galaxy_val, ng_st)

# Galaxies split by grcolor
b1_galaxy_val_grcolorA, b1_galaxy_err_grcolorA, b1_galaxy_val_grcolorB, b1_galaxy_err_grcolorB, bphi_galaxy_val_grcolorA, bphi_galaxy_err_grcolorA, bphi_galaxy_val_grcolorB, bphi_galaxy_err_grcolorB = load_bphib1_v1(z_index, ng_index, 'grcolor')
print_bias_mt(b1_galaxy_val_grcolorA, b1_galaxy_val_grcolorB, bphi_galaxy_val_grcolorA, bphi_galaxy_val_grcolorB, 'galaxies split by grcolor')
sigma_fNL_galaxy_grcolor = compute_sod_mt(b1_galaxy_val_grcolorA, b1_galaxy_val_grcolorB, bphi_galaxy_val_grcolorA, bphi_galaxy_val_grcolorB, ng_A, ng_B)

# Galaxies split by bhmmdot
b1_galaxy_val_bhmmdotA, b1_galaxy_err_bhmmdotA, b1_galaxy_val_bhmmdotB, b1_galaxy_err_bhmmdotB, bphi_galaxy_val_bhmmdotA, bphi_galaxy_err_bhmmdotA, bphi_galaxy_val_bhmmdotB, bphi_galaxy_err_bhmmdotB = load_bphib1_v1(z_index, ng_index, 'bhmmdot')
print_bias_mt(b1_galaxy_val_bhmmdotA, b1_galaxy_val_bhmmdotB, bphi_galaxy_val_bhmmdotA, bphi_galaxy_val_bhmmdotB, 'galaxies split by bhmmdot')
sigma_fNL_galaxy_bhmmdot = compute_sod_mt(b1_galaxy_val_bhmmdotA, b1_galaxy_val_bhmmdotB, bphi_galaxy_val_bhmmdotA, bphi_galaxy_val_bhmmdotB, ng_A, ng_B)

# Galaxies split by ngalax8
b1_galaxy_val_ngalax8A, b1_galaxy_err_ngalax8A, b1_galaxy_val_ngalax8B, b1_galaxy_err_ngalax8B, bphi_galaxy_val_ngalax8A, bphi_galaxy_err_ngalax8A, bphi_galaxy_val_ngalax8B, bphi_galaxy_err_ngalax8B = load_bphib1_v1(z_index, ng_index, 'ngalax8')
print_bias_mt(b1_galaxy_val_ngalax8A, b1_galaxy_val_ngalax8B, bphi_galaxy_val_ngalax8A, bphi_galaxy_val_ngalax8B, 'galaxies split by ngalax8')
sigma_fNL_galaxy_ngalax8 = compute_sod_mt(b1_galaxy_val_ngalax8A, b1_galaxy_val_ngalax8B, bphi_galaxy_val_ngalax8A, bphi_galaxy_val_ngalax8B, ng_A, ng_B)

# Galaxies split by stemass
b1_galaxy_val_stemassA, b1_galaxy_err_stemassA, b1_galaxy_val_stemassB, b1_galaxy_err_stemassB, bphi_galaxy_val_stemassA, bphi_galaxy_err_stemassA, bphi_galaxy_val_stemassB, bphi_galaxy_err_stemassB = load_bphib1_v1(z_index, ng_index, 'stemass')
print_bias_mt(b1_galaxy_val_stemassA, b1_galaxy_val_stemassB, bphi_galaxy_val_stemassA, bphi_galaxy_val_stemassB, 'galaxies split by stemass')
sigma_fNL_galaxy_stemass = compute_sod_mt(b1_galaxy_val_stemassA, b1_galaxy_val_stemassB, bphi_galaxy_val_stemassA, bphi_galaxy_val_stemassB, ng_A, ng_B)

# Galaxies split by totmass
b1_galaxy_val_totmassA, b1_galaxy_err_totmassA, b1_galaxy_val_totmassB, b1_galaxy_err_totmassB, bphi_galaxy_val_totmassA, bphi_galaxy_err_totmassA, bphi_galaxy_val_totmassB, bphi_galaxy_err_totmassB = load_bphib1_v1(z_index, ng_index, 'totmass')
print_bias_mt(b1_galaxy_val_totmassA, b1_galaxy_val_totmassB, bphi_galaxy_val_totmassA, bphi_galaxy_val_totmassB, 'galaxies split by totmass')
sigma_fNL_galaxy_totmass = compute_sod_mt(b1_galaxy_val_totmassA, b1_galaxy_val_totmassB, bphi_galaxy_val_totmassA, bphi_galaxy_val_totmassB, ng_A, ng_B)

# Full sample of halos (same so use any)
b1_halo_val   = load(pathtosimdata+'data_b1/'+list_redshift_strings[z_index]+'/data_b1_val_tng800_2_dmo_M200_'   + list_redshift_strings[z_index] + '_' + ng_strings[ng_index] + '.npy')
b1_halo_err   = load(pathtosimdata+'data_b1/'+list_redshift_strings[z_index]+'/data_b1_err_tng800_2_dmo_M200_'   + list_redshift_strings[z_index] + '_' + ng_strings[ng_index] + '.npy')
bphi_halo_val = load(pathtosimdata+'data_bphi/'+list_redshift_strings[z_index]+'/data_bphi_val_tng800_2_dmo_M200_' + list_redshift_strings[z_index] + '_' + ng_strings[ng_index] + '.npy')
bphi_halo_err = load(pathtosimdata+'data_bphi/'+list_redshift_strings[z_index]+'/data_bphi_err_tng800_2_dmo_M200_' + list_redshift_strings[z_index] + '_' + ng_strings[ng_index] + '.npy')
print_bias_st(b1_halo_val, bphi_halo_val, 'full halo sample')
sigma_fNL_halo = compute_sod_st(b1_halo_val, bphi_halo_val, ng_st)

# Halos split by c200
b1_halo_val_c200A, b1_halo_err_c200A, b1_halo_val_c200B, b1_halo_err_c200B, bphi_halo_val_c200A, bphi_halo_err_c200A, bphi_halo_val_c200B, bphi_halo_err_c200B = load_bphib1_v2(z_index, ng_index, 'c200')
print_bias_mt(b1_halo_val_c200A, b1_halo_val_c200B, bphi_halo_val_c200A, bphi_halo_val_c200B, 'galaxies split by c200')
sigma_fNL_halo_c200 = compute_sod_mt(b1_halo_val_c200A, b1_halo_val_c200B, bphi_halo_val_c200A, bphi_halo_val_c200B, ng_A, ng_B)

## ======================================================================
## Make plot
## ======================================================================

labelsize   = 36
ticksize    = 36
ticklength_major  = 10.
ticklength_minor  = 5.
tickwidth   = 1.5
tickpad     = 6.
title_font  = 30
text_font   = 30
legend_font = 22
alpha_c     = 0.4

ln_def  = 2.0
ln_def2 = 0.2
ls_def  = 'solid'
s_def   = 14

c_alll    = 'k'
c_grcolor = 'm'
c_bhmmdot = 'g'
c_stemass = 'r'
c_ngalax8 = 'c'
c_totmass = 'b'
c_c200    = 'darkorange'

m_alll = 'o'
m_A = '^'
m_B= 'v'

label_grcolor = r'$g-r$'
label_bhmmdot = r'$\dot{M}_{\rm BH}$'
label_stemass = r'$M_*$'
label_ngalax8 = r'$n_{g,8}$'
label_totmass = r'$M_{\rm t}$'
label_c200    = r'$c_{200}$'+'\n'+'(halos)'

fig1 = plt.figure(1, figsize=(17., 6.))
fig1.subplots_adjust(left=0.09, right=0.98, top=0.92, bottom=0.20, wspace = 0.33, hspace = 0.30)

# Plot bphi(b1)
panel = fig1.add_subplot(1,2,1)
title('Galaxy bias from multi-tracer splits', fontsize = title_font)
# grcolor
errorbar(b1_galaxy_val_grcolorA, bphi_galaxy_val_grcolorA, xerr=b1_galaxy_err_grcolorA, yerr=bphi_galaxy_err_grcolorA, linewidth=ln_def, linestyle=ls_def, c = c_grcolor, marker=m_A, markersize=s_def)
errorbar(b1_galaxy_val_grcolorB, bphi_galaxy_val_grcolorB, xerr=b1_galaxy_err_grcolorB, yerr=bphi_galaxy_err_grcolorB, linewidth=ln_def, linestyle=ls_def, c = c_grcolor, marker=m_B, markersize=s_def)
plot([b1_galaxy_val_grcolorA, b1_galaxy_val_grcolorB], [bphi_galaxy_val_grcolorA, bphi_galaxy_val_grcolorB], linewidth=ln_def, linestyle=ls_def, c = c_grcolor, label = label_grcolor)
# bhmmdot
errorbar(b1_galaxy_val_bhmmdotA, bphi_galaxy_val_bhmmdotA, xerr=b1_galaxy_err_bhmmdotA, yerr=bphi_galaxy_err_bhmmdotA, linewidth=ln_def, linestyle=ls_def, c = c_bhmmdot, marker=m_A, markersize=s_def)
errorbar(b1_galaxy_val_bhmmdotB, bphi_galaxy_val_bhmmdotB, xerr=b1_galaxy_err_bhmmdotB, yerr=bphi_galaxy_err_bhmmdotB, linewidth=ln_def, linestyle=ls_def, c = c_bhmmdot, marker=m_B, markersize=s_def)
plot([b1_galaxy_val_bhmmdotA, b1_galaxy_val_bhmmdotB], [bphi_galaxy_val_bhmmdotA, bphi_galaxy_val_bhmmdotB], linewidth=ln_def, linestyle=ls_def, c = c_bhmmdot, label = label_bhmmdot)
# ngalax8
errorbar(b1_galaxy_val_ngalax8A, bphi_galaxy_val_ngalax8A, xerr=b1_galaxy_err_ngalax8A, yerr=bphi_galaxy_err_ngalax8A, linewidth=ln_def, linestyle=ls_def, c = c_ngalax8, marker=m_A, markersize=s_def)
errorbar(b1_galaxy_val_ngalax8B, bphi_galaxy_val_ngalax8B, xerr=b1_galaxy_err_ngalax8B, yerr=bphi_galaxy_err_ngalax8B, linewidth=ln_def, linestyle=ls_def, c = c_ngalax8, marker=m_B, markersize=s_def)
plot([b1_galaxy_val_ngalax8A, b1_galaxy_val_ngalax8B], [bphi_galaxy_val_ngalax8A, bphi_galaxy_val_ngalax8B], linewidth=ln_def, linestyle=ls_def, c = c_ngalax8, label = label_ngalax8)
# stemass
errorbar(b1_galaxy_val_stemassA, bphi_galaxy_val_stemassA, xerr=b1_galaxy_err_stemassA, yerr=bphi_galaxy_err_stemassA, linewidth=ln_def, linestyle=ls_def, c = c_stemass, marker=m_A, markersize=s_def)
errorbar(b1_galaxy_val_stemassB, bphi_galaxy_val_stemassB, xerr=b1_galaxy_err_stemassB, yerr=bphi_galaxy_err_stemassB, linewidth=ln_def, linestyle=ls_def, c = c_stemass, marker=m_B, markersize=s_def)
plot([b1_galaxy_val_stemassA, b1_galaxy_val_stemassB], [bphi_galaxy_val_stemassA, bphi_galaxy_val_stemassB], linewidth=ln_def, linestyle=ls_def, c = c_stemass, label = label_stemass)
# totmass
errorbar(b1_galaxy_val_totmassA, bphi_galaxy_val_totmassA, xerr=b1_galaxy_err_totmassA, yerr=bphi_galaxy_err_totmassA, linewidth=ln_def, linestyle=ls_def, c = c_totmass, marker=m_A, markersize=s_def)
errorbar(b1_galaxy_val_totmassB, bphi_galaxy_val_totmassB, xerr=b1_galaxy_err_totmassB, yerr=bphi_galaxy_err_totmassB, linewidth=ln_def, linestyle=ls_def, c = c_totmass, marker=m_B, markersize=s_def)
plot([b1_galaxy_val_totmassA, b1_galaxy_val_totmassB], [bphi_galaxy_val_totmassA, bphi_galaxy_val_totmassB], linewidth=ln_def, linestyle=ls_def, c = c_totmass, label = label_totmass)
# c200
errorbar(b1_halo_val_c200A, bphi_halo_val_c200A, xerr=b1_halo_err_c200A, yerr=bphi_halo_err_c200A, linewidth=ln_def, linestyle=ls_def, c = c_c200, marker=m_A, markersize=s_def)
errorbar(b1_halo_val_c200B, bphi_halo_val_c200B, xerr=b1_halo_err_c200B, yerr=bphi_halo_err_c200B, linewidth=ln_def, linestyle=ls_def, c = c_c200, marker=m_B, markersize=s_def)
plot([b1_halo_val_c200A, b1_halo_val_c200B], [bphi_halo_val_c200A, bphi_halo_val_c200B], linewidth=ln_def, linestyle=ls_def, c = c_c200, label = label_c200)
# Cosmetics
annotate(annotate_list_redshifts[z_index]         , xy = (0.02,0.85), xycoords='axes fraction', color = 'k', fontsize = text_font-6)
annotate(ng_labels[ng_index]+r'$h^3/{\rm Mpc}^3$' , xy = (0.02,0.75), xycoords='axes fraction', color = 'k', fontsize = text_font-10)
#plot(bb1, 2.*dc*(bb1 - 1.), linestyle = 'dashed', linewidth = 2., c = 'grey')
#annotate(r'Universality', xy = (0.15, 0.07), xycoords='axes fraction', color = 'grey', fontsize = text_font-2, rotation = 0.)
locator_params(axis='x', nbins=6)
locator_params(axis='y', nbins=5)
tick_params(length=ticklength_major, width=tickwidth , bottom=True, top=True, left=True, right=True, direction = 'in', which = 'major', pad = tickpad, labelsize = ticksize)
xlabel(r'$b_1$'  , fontsize = labelsize)
ylabel(r'$b_{\phi}$', fontsize = labelsize-2)
params = {'legend.fontsize': legend_font}; pylab.rcParams.update(params); legend(loc = 'lower right', ncol = 1)

# Plot sigma_fNL
panel = fig1.add_subplot(1,2,2)
title(r'Improvements to detect $f_{\rm NL} \neq 0$', fontsize = title_font)
# Add data to plot
ratio_sigma_fNL = [sigma_fNL_galaxy_grcolor/sigma_fNL_galaxy, 
                   sigma_fNL_galaxy_stemass/sigma_fNL_galaxy, 
                   sigma_fNL_halo_c200     /sigma_fNL_halo  , 
                   sigma_fNL_galaxy_ngalax8/sigma_fNL_galaxy, 
                   sigma_fNL_galaxy_bhmmdot/sigma_fNL_galaxy, 
                   sigma_fNL_galaxy_totmass/sigma_fNL_galaxy]
x_ticknames = [label_grcolor, label_stemass, label_c200, label_ngalax8, label_bhmmdot, label_totmass]
x_ticks     = range(len(x_ticknames))
colors      = [c_grcolor, c_stemass, c_c200, c_ngalax8, c_bhmmdot, c_totmass]
for i in range(len(x_ticks)):
    scatter(x_ticks[i], ratio_sigma_fNL[i], marker = 'o', c = colors[i], s = 240)
plot(x_ticks, ratio_sigma_fNL, linewidth = 2., linestyle = 'solid', c = 'k')
# Cosmetics
ylabel(r"$\frac{{\rm SoD}\ f_{\rm NL}b_\phi^{A,B} \neq 0\ (MT)}{{\rm SoD}\ f_{\rm NL}b_\phi \neq 0\ (ST)}$", fontsize = labelsize+4)
xlim(-0.5, len(x_ticks)-1+0.5)
xticks(x_ticks, x_ticknames, fontsize = labelsize-6, rotation = 20.)
for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), colors):
    ticklabel.set_color(tickcolor)
tick_params(axis='y', length=ticklength_major, width=tickwidth,bottom=True,left=True,direction='in', which = 'major', pad = tickpad, labelsize = ticksize, labelrotation=00., labelbottom=True, labelleft=True)
annotate(r'$b_1^B b_\phi^A - b_1^A b_\phi^B$'                                                                                , xy=(0.70,0.85), xycoords='axes fraction', color = 'k'      , fontsize = text_font-8)
annotate(r'$ = %.1f$ ' % (b1_galaxy_val_grcolorB*bphi_galaxy_val_grcolorA - b1_galaxy_val_grcolorA*bphi_galaxy_val_grcolorB) , xy=(0.70,0.75), xycoords='axes fraction', color = c_grcolor, fontsize = text_font-6)
annotate(r'$ = %.1f$ ' % (b1_galaxy_val_stemassB*bphi_galaxy_val_stemassA - b1_galaxy_val_stemassA*bphi_galaxy_val_stemassB) , xy=(0.70,0.65), xycoords='axes fraction', color = c_stemass, fontsize = text_font-6)
annotate(r'$ = %.1f$ ' % (  b1_halo_val_c200B   *  bphi_halo_val_c200A    -   b1_halo_val_c200A   *  bphi_halo_val_c200B)    , xy=(0.70,0.55), xycoords='axes fraction', color = c_c200   , fontsize = text_font-6)
annotate(r'$ = %.1f$ ' % (b1_galaxy_val_ngalax8B*bphi_galaxy_val_ngalax8A - b1_galaxy_val_ngalax8A*bphi_galaxy_val_ngalax8B) , xy=(0.70,0.45), xycoords='axes fraction', color = c_ngalax8, fontsize = text_font-6)
annotate(r'$ = %.1f$ ' % (b1_galaxy_val_bhmmdotB*bphi_galaxy_val_bhmmdotA - b1_galaxy_val_bhmmdotA*bphi_galaxy_val_bhmmdotB) , xy=(0.70,0.35), xycoords='axes fraction', color = c_bhmmdot, fontsize = text_font-6)
annotate(r'$ = %.1f$ ' % (b1_galaxy_val_totmassB*bphi_galaxy_val_totmassA - b1_galaxy_val_totmassA*bphi_galaxy_val_totmassB) , xy=(0.70,0.25), xycoords='axes fraction', color = c_totmass, fontsize = text_font-6)


fig1.savefig('fig_store/fig_paper_fisher_sod_analysis_simulations.png')

#show()

