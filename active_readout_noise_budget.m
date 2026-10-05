%% Active-readout noise budget
% Design-stage input-referred noise analysis of the active
% voltage-to-current readout discussed in Sec. IV.
%
% Reproduces the numerical noise-budget results reported in the paper for
% the three input full-scale ranges. The calculation includes the intrinsic
% thermal floor of the moving-coil movement and the white and 1/f noise
% contributions of the readout electronics.

%% Physical parameters at the deterministic design point
kB=1.380649e-23; T=300;
J=3.149e-8; K=6e-7; k=6.283e-3; b0=4.17e-8; %#ok<NASGU>
Rc=15.28; Rp=154.0; R=Rc+Rp; %#ok<NASGU>
I_FS=100e-6; I_th_RMS=sqrt(kB*T*K)/k;

%% Critically damped mechanical response
omega_n=sqrt(K/J); f_n=omega_n/(2*pi); B_eq=pi*f_n/4;
H2=@(f) 1./((1-(f/f_n).^2).^2+(2*(f/f_n)).^2);

%% Representative precision OP-AMP noise parameters
e_nw=8e-9; f_ce=3; i_nw=0.2e-12; f_ci=20; GBW=1e6;

%% Readout parameters
R_add=100e3; Rsense_list=[10,100,1000]; T_obs=5; f_L=1/T_obs;

%% Electronics output-current noise PSD
SI_white=@(Rs) e_nw^2/Rs^2+i_nw^2+4*kB*T/Rs;
SI_1f=@(f,Rs) e_nw^2*(1+f_ce./f)/Rs^2+i_nw^2*(1+f_ci./f)+4*kB*T/Rs;

%% Storage for table quantities
nR=numel(Rsense_list); V_FS=zeros(nR,1); V_th_RMS=zeros(nR,1);
V_tot_white=zeros(nR,1); V_tot_1f=zeros(nR,1); DR_dB=zeros(nR,1); noise_ppm=zeros(nR,1);

%% Noise budget
for m=1:nR
    Rs=Rsense_list(m);
    V_FS(m)=I_FS*Rs;
    V_th_RMS(m)=I_th_RMS*Rs;
    I_el_white=sqrt(SI_white(Rs)*B_eq);
    I_tot_white=sqrt(I_th_RMS^2+I_el_white^2);
    V_tot_white(m)=Rs*I_tot_white;
    var_el_1f=integral(@(u) SI_1f(exp(u),Rs).*H2(exp(u)).*exp(u),log(f_L),log(1e8),'AbsTol',1e-40,'RelTol',1e-10);
    I_el_1f=sqrt(var_el_1f);
    I_tot_1f=sqrt(I_th_RMS^2+I_el_1f^2);
    V_tot_1f(m)=Rs*I_tot_1f;
    DR_dB(m)=20*log10(V_FS(m)/V_tot_1f(m));
    noise_ppm(m)=1e6*V_tot_1f(m)/V_FS(m);
end

%% Sanity checks
fprintf('f_n            = %.4f Hz\n',f_n)
fprintf('B_eq           = %.4f Hz\n',B_eq)
fprintf('I_th,RMS       = %.3f pA\n',I_th_RMS*1e12)
fprintf('T_obs          = %.1f s\n',T_obs)
fprintf('f_L            = %.3f Hz\n\n',f_L)
fprintf('Electrical closed-loop bandwidth check:\n')
for m=1:nR
    Rs=Rsense_list(m); NG=1+R_add/Rs;
    fprintf('R_sense = %4g Ohm: NG = %7.1f, GBW/NG = %8.1f Hz\n',Rs,NG,GBW/NG)
end

%% Values reported in the readout noise-budget table
fprintf('\n=== Readout noise budget ===\n')
fprintf(['R_sense   V_FS        V_th,RMS     V_tot white   ' ...
         'V_tot 1/f     DR 1/f     Noise/FS\n'])
for m=1:nR
    fprintf(['%4g Ohm   %7.3f mV   %8.3f nV   %10.2f nV   ' ...
             '%9.2f nV   %6.1f dB   %7.3f ppm\n'], ...
        Rsense_list(m),V_FS(m)*1e3,V_th_RMS(m)*1e9,V_tot_white(m)*1e9, ...
        V_tot_1f(m)*1e9,DR_dB(m),noise_ppm(m))
end
