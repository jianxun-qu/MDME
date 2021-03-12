
clear all
close all

csf_m0 = 1;
csf_t1 = 4000;
csf_t2 = 2000;

gm_m0 = 0.9;
gm_t1 = 1450;
gm_t2 = 80;

TR = 4260;
TD = [150, 580, 2000, 4130];
TE = [22, 100];

alpha = 90/180*pi; % excitation flip angle (90)
theta = 120/180*pi; % saturation flip angle (120)
b1scale = 1; % rf field scaling

% % Partial Volume

pvr_arr = 0:0.02:1;

mix_m0_arr = zeros(size(pvr_arr));
mix_t1_arr = zeros(size(pvr_arr));
mix_t2_arr = zeros(size(pvr_arr));

% % [Sig1, ..., Sig8] <---> [qT1, qT2]

cos_theta = cos(theta);
cos_alpha = cos(alpha);

fopts = fitoptions(...
    'Method', 'NonlinearLeastSquares',...
    'Lower', [0, 100, 10],...
    'Upper', [2, 4000, 4000],...
    'Startpoint', [1.0, 1000, 100],...
    'DiffMinChange', 1.0e-6,...
    'DiffMaxChange', 0.01,...
    'MaxIter', 2000);

ftype = fittype(...
    'm0*(1-(1-cos_theta)*exp(-Td/T1)-cos_theta*exp(-TR/T1))/(1-cos_theta*cos_alpha*exp(-TR/T1))*exp(-TE/T2)',...
    'dependent', {'s'},...
    'independent', {'Td', 'TE'},...
    'coefficients', {'m0', 'T1', 'T2'},...
    'problem', {'TR', 'cos_theta', 'cos_alpha'});

TE_4fit = [TE', TE', TE', TE'];
TD_4fit = [TD; TD];

gm_sig_full = auxil_mdme_sig(gm_m0, gm_t1, gm_t2, TD, TE, TR, alpha, theta, b1scale);
csf_sig_full = auxil_mdme_sig(csf_m0, csf_t1, csf_t2, TD, TE, TR, alpha, theta, b1scale);

for pvr_idx = 1: length(pvr_arr)
    
    pvr = pvr_arr(pvr_idx);
    
    gm_sig = (1-pvr) * gm_sig_full;
    csf_sig = pvr * csf_sig_full;
    mix_sig = gm_sig + csf_sig;
    
    fobj = fit([TD_4fit(:), TE_4fit(:)],  mix_sig(:), ftype, fopts, 'problem', {TR, cos_theta, cos_alpha});
    
    mix_m0 = fobj.m0;
    mix_t1 = fobj.T1;
    mix_t2 = fobj.T2;
    
    mix_m0_arr(pvr_idx) = mix_m0;
    mix_t1_arr(pvr_idx) = mix_t1;
    mix_t2_arr(pvr_idx) = mix_t2;
    
end

% % Syn FLAIR

flair_ti = 1966;
flair_tr = 6000;
flair_te = 200;

csf_flair_arr = [];
gm_flair_arr = [];
mix_flair_arr = [];


csf_flair = auxil_mdme_syflair(flair_ti, flair_te, flair_tr, csf_m0*pvr_arr, csf_t1, csf_t2);
gm_flair = auxil_mdme_syflair(flair_ti, flair_te, flair_tr, gm_m0*(1-pvr_arr), gm_t1, gm_t2);
mix_flair = auxil_mdme_syflair(flair_ti, flair_te, flair_tr, mix_m0_arr, mix_t1_arr, mix_t2_arr); 

csf_flair = csf_flair / gm_flair(1);
gm_flair = gm_flair / gm_flair(1);
mix_flair = mix_flair / mix_flair(1);

plot(pvr_arr, csf_flair + gm_flair, 'LineWidth', 6); hold on
plot(pvr_arr, mix_flair, 'LineWidth', 6);
xlim([0 1])
ylim([0 2])
grid on
set(gcf, 'Color', 'w');
set(gca, 'FontName','Calibri', 'FontSize', 30, 'FontWeight', 'BOLD')
title('FLAIR signal change')
xlabel('CSF Partial Volume Ratio')
ylabel('FLAIR Signal Intensity')
legend('Actual', 'SynFLAIR')



