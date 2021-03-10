
% clear all
% close all

csfm0 = 1;
csft1 = 4000;
csft2 = 2000;

gmm0 = 0.9;
gmt1 = 1450;
gmt2 = 80;

TR = 4260;
TD = [150, 580, 2000, 4130];
TE = [22, 100];

alpha = 90/180*pi; % excitation flip angle (90)
theta = 120/180*pi; % saturation flip angle (120)
b1scale = 1; % rf field scaling

figure(1)

TDi = 0:TR;
csfsig = auxil_mdme_mtd(1, csft1, csft2, TDi, 0, TR, alpha, theta, b1scale);
gmsig = auxil_mdme_mtd(1, gmt1, gmt2, TDi, 0, TR, alpha, theta, b1scale);

plot(TDi, csfsig); xlim([0 TR]), ylim([-1 1]), grid on, hold on
plot(TDi, csfsig * exp(-TE(1)/csft2)); xlim([0 TR]), ylim([-1 1]), grid on, hold on
plot(TDi, csfsig * exp(-TE(2)/csft2)); xlim([0 TR]), ylim([-1 1]), grid on, hold on
plot(TDi, gmsig); xlim([0 TR]), ylim([-1 1]), grid on
plot(TDi, gmsig * exp(-TE(1)/gmt2)); xlim([0 TR]), ylim([-1 1]), grid on
plot(TDi, gmsig * exp(-TE(2)/gmt2)); xlim([0 TR]), ylim([-1 1]), grid on


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

% T2FLAIR
flairti = 1965;
flairtr = 6000;
flairte = 250;

% % T1FLAIR
% flairti = 876;
% flairtr = 2000;
% flairte = 240;

csfflairarr = [];
gmflairarr = [];
mixflairarr = [];

mixm0arr = [];
mixt1arr = [];
mixt2arr = [];

pvrarr = 0:0.02:1;

for pvr = pvrarr
    
    csfsig0 = csfm0 * pvr;
    csfsig = auxil_mdme_mtd(csfsig0, csft1, csft2, TD, 0, TR, alpha, theta, b1scale);
    csfsig = exp(-TE./csft2)' * csfsig;
    
    gmsig0 = gmm0 * (1-pvr);
    gmsig = auxil_mdme_mtd(gmsig0, gmt1, gmt2, TD, 0, TR, alpha, theta, b1scale);
    gmsig = exp(-TE./gmt2)' * gmsig;
    
    csig = gmsig' + csfsig';
    te = [TE; TE; TE; TE];
    td = [TD'; TD'];
    
    fobj = fit([td(:), te(:)],  csig(:), ftype, fopts, 'problem', {TR, cos_theta, cos_alpha});
    
    mixm0 = fobj.m0;
    mixt1 = fobj.T1;
    mixt2 = fobj.T2;
    
    csfflair = csfsig0 * (1-2*exp(-flairti/csft1)+exp(-flairtr/csft1))*exp(-flairte/csft2);
    gmflair = gmsig0 * (1-2*exp(-flairti/gmt1)+exp(-flairtr/gmt1))*exp(-flairte/gmt2);

    mixflair = mixm0 * (1-2*exp(-flairti/mixt1)+exp(-flairtr/mixt1))*exp(-flairte/mixt2);
    
    csfflairarr = [csfflairarr, csfflair];
    gmflairarr = [gmflairarr, gmflair];
    mixflairarr = [mixflairarr, mixflair];
    
    mixm0arr = [mixm0arr, mixm0];
    mixt1arr = [mixt1arr, mixt1];
    mixt2arr = [mixt2arr, mixt2];
    
end

% figure(2)
% plot(csfflairarr); hold on
% plot(gmflairarr); hold on

truesumarr = csfflairarr+gmflairarr;

% plot(pvrarr, truesumarr/truesumarr(1)); hold on
plot(pvrarr, mixflairarr/mixflairarr(1)); hold on
legend

% figure(3),
% plot(pvrarr, mixm0arr/1); hold on
% plot(pvrarr, mixt1arr/csft1);
% plot(pvrarr, mixt2arr/csft2);