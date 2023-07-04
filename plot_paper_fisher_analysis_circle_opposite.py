from commons import *

## ======================================================================
## Compute sigma_fNL on a circle bphi-b1 for multi tracer study
## ======================================================================

# Choose desired number densities here
ng_F1   = 2.0e-04 # h^3/Mpc^3
ng_F2   = 2.0e-00 # h^3/Mpc^3
label_mt1  = r'$\bar{n}_g = 2 \times 10^{-4}$\ ; $\bar{n}_g^A,\bar{n}_g^B=\bar{n}_g/2$'
label_mt2  = r'$\bar{n}_g = 2$\ ; $\bar{n}_g^A,\bar{n}_g^B=\bar{n}_g/2$'

# Choose desired bias parameters of the full sample
b1_F   = 2.
bphi_F = 2.*dc*(b1_F - 1.)

# Choose desired number density of the two subsamples
ng_A1 = ng_F1/2.
ng_B1 = ng_F1/2.
ng_A2 = ng_F2/2.
ng_B2 = ng_F2/2.

# Choose details of the circle of multi-tracer splits
Nttheta = 50
ttheta = linspace(-pi/2., pi/2., Nttheta)
Radius = 1.

# Compute the multi-tracer Fisher matrix constraints
sigma_fNL1_mt = zeros(Nttheta)
sigma_fNL2_mt = zeros(Nttheta)
bb1_A   = zeros(Nttheta)
bbphi_A = zeros(Nttheta)
bb1_B   = zeros(Nttheta)
bbphi_B = zeros(Nttheta)
for i1 in range(Nttheta):
    theta_now = ttheta[i1]
    b1_A   = b1_F   + Radius*cos(theta_now+pi)
    bphi_A = bphi_F + Radius*sin(theta_now+pi)
    b1_B   = b1_F   + Radius*cos(theta_now)
    bphi_B = bphi_F + Radius*sin(theta_now)
    Finfo1_now = zeros(len(kk))
    Finfo2_now = zeros(len(kk))
    for i2 in range(len(kk)):
        Finfo1_now[i2] = dot(model_derivative_multitracer(z, kk[i2], fNL, b1_A, b1_B, bphi_A, bphi_B),
                         dot(linalg.inv(covariance_multitracer(z, kk[i2], fNL, b1_A, b1_B, bphi_A, bphi_B, ng_A1, ng_B1, Vs, kF)),
                             model_derivative_multitracer(z, kk[i2], fNL, b1_A, b1_B, bphi_A, bphi_B)))
        Finfo2_now[i2] = dot(model_derivative_multitracer(z, kk[i2], fNL, b1_A, b1_B, bphi_A, bphi_B),
                         dot(linalg.inv(covariance_multitracer(z, kk[i2], fNL, b1_A, b1_B, bphi_A, bphi_B, ng_A2, ng_B2, Vs, kF)),
                             model_derivative_multitracer(z, kk[i2], fNL, b1_A, b1_B, bphi_A, bphi_B)))
    sigma_fNL1_mt[i1] = 1./sqrt(sum(Finfo1_now))
    sigma_fNL2_mt[i1] = 1./sqrt(sum(Finfo2_now))
    bb1_A[i1]        = b1_A 
    bbphi_A[i1]      = bphi_A
    bb1_B[i1]        = b1_B
    bbphi_B[i1]      = bphi_B

# Compute the single-tracer Fisher matrix constraints
Finfo1_now = zeros(len(kk))
Finfo2_now = zeros(len(kk))
for i3 in range(len(kk)):
    Finfo1_now[i3] = (model_derivative_singletracer(z, kk[i3], fNL, b1_F, bphi_F)**2. / covariance_singletracer(z, kk[i3], fNL, b1_F, bphi_F, ng_F1, Vs, kF))
    Finfo2_now[i3] = (model_derivative_singletracer(z, kk[i3], fNL, b1_F, bphi_F)**2. / covariance_singletracer(z, kk[i3], fNL, b1_F, bphi_F, ng_F2, Vs, kF))
sigma_fNL1_F = 1./sqrt(sum(Finfo1_now))
sigma_fNL2_F = 1./sqrt(sum(Finfo2_now))

## ======================================================================
## Run analysis and make plot
## ======================================================================

labelsize   = 30
ticksize    = 30
ticklength_major  = 10.
ticklength_minor  = 5.
tickwidth   = 1.5
tickpad     = 6.
title_font  = 30
text_font   = 20
legend_font = 20
alpha_c     = 0.4

m_def = 'o'
s_def = 6

bb1 = linspace(b1_F - Radius*5, b1_F + Radius*5, 100)

fig1 = plt.figure(1, figsize=(16., 6.))
fig1.subplots_adjust(left=0.06, right=0.93, top=0.90, bottom=0.15, wspace = 0.30, hspace = 0.30)

# One number density
panel = fig1.add_subplot(1,2,2)
title(label_mt1, fontsize= title_font)
panel.set_aspect('equal')
scatter(b1_F, bphi_F, s = 80, c = 'k')
scatter(bb1_A, bbphi_A, c = sigma_fNL1_mt/sigma_fNL1_F, s = 80, cmap = 'jet')
scatter(bb1_B, bbphi_B, c = sigma_fNL1_mt/sigma_fNL1_F, s = 80, cmap = 'jet')
plot(bb1, 2.*dc*(bb1-1), linestyle = 'dashed', c = 'grey', linewidth = 2., label = r'Universality')
cb = colorbar()
cb.set_label(r'$\sigma_{f_{\rm NL}}^{\rm MT}/\sigma_{f_{\rm NL}}^{\rm ST}$', fontsize = labelsize-2)
cb.ax.tick_params(labelsize=ticksize-4)
xlim(b1_F - Radius*1.5, b1_F + Radius*1.5)
ylim(bphi_F - Radius*1.5, bphi_F + Radius*1.5)
xlabel(r'$b_1$'  , fontsize = labelsize)
ylabel(r'$b_{\phi}$', fontsize = labelsize)
annotate(r'Full sample', xy = (0.52,0.48), xycoords='axes fraction', color = 'k', fontsize = text_font)
annotate(r'Samples A'  , xy = (0.22,0.87), xycoords='axes fraction', color = 'k', fontsize = text_font)
annotate(r'Samples B'  , xy = (0.70,0.15), xycoords='axes fraction', color = 'k', fontsize = text_font)
annotate(r'The $A,B$ pairs are opposed in the circle'   , xy = (0.01,0.05), xycoords='axes fraction', color = 'k', fontsize = text_font-1)
annotate(r'Universality', xy = (0.640, 0.89), xycoords='axes fraction', color = 'grey', fontsize = text_font, rotation = 0.)
tick_params(length=ticklength_major, width=tickwidth , bottom=True, top=True, left=True, right=True, direction = 'in', which = 'major', pad = tickpad, labelsize = ticksize)

# The other number density
panel = fig1.add_subplot(1,2,1)
title(label_mt2, fontsize= title_font)
panel.set_aspect('equal')
scatter(b1_F, bphi_F, s = 80, c = 'k')
scatter(bb1_A, bbphi_A, c = log10(sigma_fNL2_mt/sigma_fNL2_F), s = 80, cmap = 'jet')
scatter(bb1_B, bbphi_B, c = log10(sigma_fNL2_mt/sigma_fNL2_F), s = 80, cmap = 'jet')
plot(bb1, 2.*dc*(bb1-1), linestyle = 'dashed', c = 'grey', linewidth = 2.)
cb = colorbar()
cb.set_label(r'${\rm log}_{10}\ \sigma_{f_{\rm NL}}^{\rm MT}/\sigma_{f_{\rm NL}}^{\rm ST}$', fontsize = labelsize-2)
cb.ax.tick_params(labelsize=ticksize-4)
xlim(b1_F - Radius*1.5, b1_F + Radius*1.5)
ylim(bphi_F - Radius*1.5, bphi_F + Radius*1.5)
xlabel(r'$b_1$'  , fontsize = labelsize)
ylabel(r'$b_{\phi}$', fontsize = labelsize)
annotate(r'Full sample', xy = (0.52,0.48), xycoords='axes fraction', color = 'k', fontsize = text_font)
annotate(r'Samples A'  , xy = (0.22,0.87), xycoords='axes fraction', color = 'k', fontsize = text_font)
annotate(r'Samples B'  , xy = (0.70,0.15), xycoords='axes fraction', color = 'k', fontsize = text_font)
annotate(r'The $A,B$ pairs are opposed in the circle'   , xy = (0.01,0.05), xycoords='axes fraction', color = 'k', fontsize = text_font-1)
annotate(r'Universality', xy = (0.640, 0.89), xycoords='axes fraction', color = 'grey', fontsize = text_font, rotation = 0.)
tick_params(length=ticklength_major, width=tickwidth , bottom=True, top=True, left=True, right=True, direction = 'in', which = 'major', pad = tickpad, labelsize = ticksize)

fig1.savefig('fig_store/fig_paper_fisher_fNL_analysis_circle_opposite.png')

#show()


