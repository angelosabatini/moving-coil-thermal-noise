%% Monte Carlo verification of the deflection-angle ACF
% Exact-discrete Monte Carlo simulation of the third-order stochastic
% moving-coil movement model discussed in Sec. III.
%
% Reproduces the numerical checks and the deflection-angle ACF shown in
% Fig. 2. The parameters of the numerical example are reconstructed from
% the movement geometry, and the critical-damping resistance is computed
% from the deterministic model.
%
% The continuous-time Ornstein-Uhlenbeck process is propagated exactly
% between stored samples using the Van Loan discretization. Each
% realization is initialized directly from the analytical equilibrium
% covariance, so no burn-in interval or Euler-Maruyama approximation is
% used.
%
% The script compares the complete third-order Monte Carlo model with the
% reduced critically damped second-order ACF and checks the stationary
% covariance predicted by equipartition.

rng(1234)

%% Geometry of the numerical example
N=100; a=1e-2; ell=1.5708e-2; r=50e-6; B=0.2;
rho_m=8960; rho_e=1.68e-8; mu0=4*pi*1e-7; w=2*a;

%% Movement parameters derived from the geometry
k = 2*N*a*ell*B;
J = 2*N*rho_m*pi*r^2*a^2*(ell + (2/3)*a);
Rc = rho_e*N*(2*ell + 2*w)/(pi*r^2);
Lself = @(x) (mu0/(2*pi))*x*(log(2*x/r)-1);
Mmut = @(x,d) (mu0/(2*pi))*(x*asinh(x/d)-sqrt(x^2+d^2)+d);
L = N^2*(2*Lself(ell)+2*Lself(w)-2*Mmut(ell,w)-2*Mmut(w,ell));

%% Remaining mechanical and thermal parameters
K=6e-7; b0=4.17e-8; T=300; kB=1.380649e-23;

%% Critical-damping electrical termination
R_CDRX = k^2/(2*sqrt(J*K)-b0)-Rc;
R = Rc + R_CDRX;
S=k/K; omega_n=sqrt(K/J); f_n=omega_n/(2*pi);
B_3dB=f_n*sqrt(sqrt(2)-1);

fprintf('=== Parameters derived for the numerical example ===\n')
fprintf('Rc       = %.6f Ohm\n',Rc); fprintf('L        = %.6e H\n',L)
fprintf('k        = %.6e N m/A\n',k); fprintf('J        = %.6e kg m^2\n',J)
fprintf('b0       = %.6e N m s\n',b0); fprintf('K        = %.6e N m/rad\n',K)
fprintf('R_CDRX   = %.6f Ohm\n',R_CDRX); fprintf('R        = %.6f Ohm\n',R)
fprintf('S        = %.6e rad/A\n',S); fprintf('f_n      = %.6f Hz\n',f_n)
fprintf('B_3dB    = %.6f Hz\n\n',B_3dB)

%% Continuous-time stochastic state-space model
A=[0,1,0;-K/J,-b0/J,k/J;0,-k/L,-R/L];
G=[0,0;sqrt(2*kB*T*b0)/J,0;0,sqrt(2*kB*T*R)/L];
D=G*G.';

%% Exact equilibrium covariance
M=diag([K,J,L]); %#ok<NASGU>
Sigma_eq=kB*T*diag([1/K,1/J,1/L]);
theta_var_exact=Sigma_eq(1,1); theta_rms_exact=sqrt(theta_var_exact);
Ith_rms_exact=(K/k)*theta_rms_exact;

fprintf('=== Exact equilibrium result ===\n')
fprintf('<theta^2> = %.6e rad^2\n',theta_var_exact)
fprintf('theta_RMS = %.6e rad\n',theta_rms_exact)
fprintf('I_th,RMS  = %.6e A  (%.3f pA)\n\n',Ith_rms_exact,Ith_rms_exact*1e12)

%% Reduced second-order parameters and time-scale separation
b_em=k^2/R; zeta=(b0+b_em)/(2*sqrt(J*K)); tau_e=L/R; tau_m=1/omega_n;
fprintf('=== Time scales and damping ===\n')
fprintf('omega_n          = %.6f rad/s\n',omega_n); fprintf('f_n              = %.6f Hz\n',f_n)
fprintf('L/R              = %.6e s\n',tau_e); fprintf('1/omega_n        = %.6e s\n',tau_m)
fprintf('(1/omega_n)/(L/R)= %.3e\n',tau_m/tau_e); fprintf('b_em             = %.6e N m s\n',b_em)
fprintf('zeta             = %.12f\n\n',zeta)

%% Exact discrete-time propagation: Van Loan method
h=1e-4; n=size(A,1); Z=[-A,D;zeros(n),A.']*h; E=expm(Z);
F=E(n+1:end,n+1:end).'; Qd=F*E(1:n,n+1:end); Qd=(Qd+Qd.')/2;
stationarity_error=norm(Sigma_eq-F*Sigma_eq*F.'-Qd,'fro')/norm(Sigma_eq,'fro');
fprintf('=== Exact discretization check ===\n')
fprintf('relative stationary-covariance residual = %.3e\n\n',stationarity_error)
Lq=chol(Qd,'lower'); Leq=chol(Sigma_eq,'lower');

%% Monte Carlo settings
tau_slow=1/min(abs(real(eig(A)))); nPerRep=round(40*tau_slow/h); nRep=500;
maxLag=round(8*tau_slow/h); lags=(0:maxLag).'; tLag=lags*h;
fprintf('=== Monte Carlo settings ===\n')
fprintf('slowest decay time = %.6f s\n',tau_slow)
fprintf('record length      = %.6f s per replication\n',nPerRep*h)
fprintf('replications       = %d\n',nRep); fprintf('maximum ACF lag    = %.6f s\n\n',maxLag*h)

%% Exact-discrete Monte Carlo simulation
acfMat=zeros(maxLag+1,nRep); sum_x=zeros(3,1); sum_xx=zeros(3,3); Nstat=0;
for rr=1:nRep
    x=Leq*randn(3,1); theta=zeros(nPerRep,1);
    for m=1:nPerRep
        theta(m)=x(1); sum_x=sum_x+x; sum_xx=sum_xx+x*x.'; Nstat=Nstat+1;
        if m<nPerRep, x=F*x+Lq*randn(3,1); end
    end
    acfFull=xcorr(theta,maxLag,'biased'); acfMat(:,rr)=acfFull(maxLag+1:end);
end

%% Stationary variance and covariance checks
mu_mc=sum_x/Nstat; Sigma_mc=sum_xx/Nstat-mu_mc*mu_mc.';
fprintf('=== Stationary covariance check ===\n'); fprintf('Monte Carlo mean:\n'); disp(mu_mc)
fprintf('Empirical stationary covariance:\n'); disp(Sigma_mc)
fprintf('Exact covariance kB*T*M^(-1):\n'); disp(Sigma_eq)
relErrDiag=100*(diag(Sigma_mc)-diag(Sigma_eq))./diag(Sigma_eq);
fprintf('relative error on diagonal (theta^2, omega^2, i^2) [%%]:\n'); disp(relErrDiag.')
theta_rms_mc=sqrt(Sigma_mc(1,1)); Ith_rms_mc=(K/k)*theta_rms_mc;
fprintf('theta_RMS, exact       = %.6e rad\n',theta_rms_exact)
fprintf('theta_RMS, Monte Carlo = %.6e rad\n',theta_rms_mc)
fprintf('I_th,RMS, exact        = %.6e A  (%.3f pA)\n',Ith_rms_exact,Ith_rms_exact*1e12)
fprintf('I_th,RMS, Monte Carlo  = %.6e A  (%.3f pA)\n\n',Ith_rms_mc,Ith_rms_mc*1e12)

%% Correlation-matrix check
s=sqrt(diag(Sigma_mc)); Corr_mc=Sigma_mc./(s*s.');
fprintf('=== Monte Carlo correlation matrix ===\n'); disp(Corr_mc)
offDiagCorr=Corr_mc-diag(diag(Corr_mc));
fprintf('max |off-diagonal correlation| = %.4f  (%.2f%%)\n\n',max(abs(offDiagCorr(:))),100*max(abs(offDiagCorr(:))))

%% Deflection-angle ACF
acfMean=mean(acfMat,2); acfLow=prctile(acfMat,2.5,2); acfHigh=prctile(acfMat,97.5,2);
tol_zeta=1e-10;
if zeta < 1-tol_zeta
    omega_d=omega_n*sqrt(1-zeta^2);
    acfTheory=(kB*T/K).*exp(-zeta*omega_n*tLag).*(cos(omega_d*tLag)+(zeta*omega_n/omega_d).*sin(omega_d*tLag));
elseif zeta > 1+tol_zeta
    alpha=omega_n*sqrt(zeta^2-1);
    acfTheory=(kB*T/K).*exp(-zeta*omega_n*tLag).*(cosh(alpha*tLag)+(zeta*omega_n/alpha).*sinh(alpha*tLag));
else
    acfTheory=(kB*T/K).*exp(-omega_n*tLag).*(1+omega_n*tLag);
end

%% ACF figure
figure; hold on
fill([tLag(2:end);flipud(tLag(2:end))],[acfLow(2:end);flipud(acfHigh(2:end))],[0.85 0.85 0.85],'EdgeColor','none')
plot(tLag(1:750:end),acfTheory(1:750:end),'k.','MarkerSize',14)
plot(tLag,acfMean,'r','LineWidth',1.5)
xlabel('Lag \tau [s]'); ylabel('Deflection angle ACF [rad^2]')
legend('95% Monte Carlo band','Reduced-model theory','Monte Carlo mean','Location','northeast','Box','off')
grid off; set(gca,'FontSize',14,'LineWidth',1.5); xlim([0 1.5])
