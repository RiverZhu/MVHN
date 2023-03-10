% Demo for the MVALSE algorithm in Heteroscedastic Noise Environment
% For details, refer to the paper
% Q. Zhang, J. Zhu*, Y. Gu and Z. Xu, Grid-less variational direction of arrival estimation in heteroscedastic noise environment, 
% IEEE Journal of Oceanic Engineering, vol. 46, no. 4, pp. 1313-1329, 2021. 
% This code is written by Qi Zhang and Jiang Zhu. If you have any problems, please feel free to contact
% jiangzhu16@zji.edu.cn.
clear variables 
rng(123);

%% Generate the signal Y
M   = 50;                           % number of measurements for each observation
K   = 3;                            % number of complex sinusoids
N   = 50;                           % assumed number of complex sinusoids
L   = 10;                           % number of snapshots 
d    = pi/N;                        % minimum separation of the angular frequencies
SNR = 0;                           % signal to noise ratio (dB)
CASE = 0;
delta_dB = 10;
omega = [0.3;0.5;1.2];                 % true theta


%% Generate the signal Y
A   = exp(1i*(0:1:M-1).'*omega.');        % matrix with columns a(omega(1)),..., a(omega(K))
R   = 1 + .2.*randn(K,L);                        % magnitudes of the complex amplitudes
W  = R.*exp(1i*2*pi*rand(K,L));                  % complex amplitudes   
X   = A*W;                                       % original signal
nominal_Pn = 10*log10(mean(mean(abs(X).^2))*10^(-SNR/10));
nominal_Pn_max = nominal_Pn + delta_dB;
nominal_Pn_min = nominal_Pn - delta_dB;
switch CASE
    case 0
        Pn_dB = nominal_Pn;
        Pn  = (10^(Pn_dB/10))*ones(M,L);     % noise power
    case 1
        Pn_dB = (nominal_Pn_max - nominal_Pn_min)*rand(1,L) + nominal_Pn_min;     
        Pn  = repmat(10.^(Pn_dB/10),M,1);      
    case 2
        Pn_dB = (nominal_Pn_max - nominal_Pn_min)*rand(M,1) + nominal_Pn_min;     
        Pn  = repmat(10.^(Pn_dB/10),1,L);     
    case 3
        Pn_dB = (nominal_Pn_max - nominal_Pn_min)*rand(M,L) + nominal_Pn_min;     
        Pn  = 10.^(Pn_dB/10);                 
end
eps = sqrt(0.5*Pn).*(randn(M,L)+1i*randn(M,L));  % complex Gaussian noise
Y   = X + eps;                                   % measured signal
    
%% plot the line spectral signal
h = figure;
set(h,'position',[100 100 600 400]); 
hold on;
axis([-pi pi -inf inf]); xlabel('Frequencies'); ylabel('Magnitudes');
stem(omega,mean(abs(W),2),'*k','markersize',8,'LineWidth',1.5);
ylim([0 max(mean(abs(W)))*1.2]);

%% MVALSE algorithm in Heteroscedastic Noise Environment
% The sixth argument of MVALSE_HN is to determine the model order
% criterion method. 0: The method is adopted in the published IEEE JOE
% paper, and is more likely to activate the signal. 1: The method is
% adopted in the published IEEE TAES paper N. Zhang, J. Zhu and Z. Xu, 
% Gridless multisnapshot variational line spectral estimation from coarsely quantized samples,  
% IEEE Transactions on Aerospace and Electronic Systems.
tic;
outMultiVALSE1 = MVALSE_HN( Y, N, X, CASE, 0, 0 );
toc;
fprintf('Runtime of MultiVALSE: %g s.\n',toc);
stem(outMultiVALSE1.freqs,mean(abs(outMultiVALSE1.amps),2),'or','markersize',8,'LineWidth',1.5); 

%% New 
tic;
outMultiVALSE2 = MVALSE_HN( Y, N, X, CASE, 0, 1 );
toc;
fprintf('Runtime of MultiVALSE: %g s.\n',toc);
stem(outMultiVALSE2.freqs,mean(abs(outMultiVALSE2.amps),2),'sb','markersize',8,'LineWidth',1.5); 
legend('true', 'MVALSE_HN', 'MVALSE_HN_MD')

