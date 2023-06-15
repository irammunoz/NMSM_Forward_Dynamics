function  [out, m1w_p, m2w_p] = FDcrb_NMSM( model, y, wrapc, e, tau, fe, wrap_flag, f_ext)

% FDcrb  Forward Dynamics via Composite-Rigid-Body Algorithm
% FDcrb(model,q,qd,tau,f_ext)  calculates the forward dynamics of a
% kinematic tree via the composite-rigid-body algorithm.  q, qd and tau are
% vectors of joint position, velocity and force variables; and the return
% value is a vector of joint acceleration variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.

NMUS = length(model.m);
NDOF = model.NB;

q = y(1:NDOF); qd = y(NDOF+1:NDOF*2);
a = y(NDOF*2+1:NDOF*2+NMUS); lnM = y(NDOF*2+NMUS+1:NDOF*2+NMUS*2);

if nargin == 7
  [H,C,Xup,Xa,V,J] = HandC( model, q, qd );
else
  [H,C,Xup,Xa,V,J] = HandC( model, q, qd, f_ext );
end

% Kinematics of the upper limb's wrist
Xee = model.KXtree{1} * Xa{model.kparent};
Jg_n = [eye(3) zeros(3,3);model.KXtree{1}(4:6,1:3) eye(3)]*J{end}; % local jacobian
Jg_e = [Xee(1:3,1:3),zeros(3,3);zeros(3,3),Xee(1:3,1:3)]*Jg_n;     % geometric jacobian

% Use when there are external forces
% Apply external forces to the system
% fe = Xee'\fe;
% tau_ext = Jg_n'*fe;

% Use when there are no external forces
fe = zeros(6,1);
tau_ext = zeros(model.NB,1);

% Compute the viscoelastic passive torques to limit the range of motion
p_tau = passiveTorques(model,q,qd);

% Add the contribution of passive torques and external forces to the 
% generalized forces  
tau = tau + p_tau + tau_ext;

% Unpackage the muscle-wrapping parameters for Muscle 1 "BRA"
u0 = wrapc(1); v0 = wrapc(2); up0 = wrapc(3); vp0 = wrapc(4); lg0 = wrapc(5);
p_l = wrapc(6:8);
q_l = wrapc(9:11);
uq0 = wrapc(12); vq0 = wrapc(13); upq0 = wrapc(14); vpq0 = wrapc(15);

if wrap_flag == 1
    % Broyden method
    [lmt(1),lmtd(1),u(1),v(1),up(1),vp(1),lg(1),p_l(:,1),q_l(:,1),uq(1),vq(1),upq(1),vpq(1),iter(1),event(1),Jm1] = SolveWrapping_BM(model,Xup,V,J,1,u0,v0,up0,vp0,lg0,p_l,q_l,uq0,vq0,upq0,vpq0);
    m1w_p = 0;
else
    % Obstcle-Set Method Analytical solution for the cylinder
    [lmt(1),lmtd(1),event(1),Jm1,m1w_p] = SolveWrapping_OS(model,Xup,V,J,1);
    % These variables are not used
    u(1) = 0; v(1) = 0; up(1) = 0; vp(1) = 0; lg(1) = 0;
    p_l(:,1) = zeros(3,1); q_l(:,1) = zeros(3,1); 
    uq(1) = 0; vq(1) = 0; upq(1) = 0; vpq(1) = 0; iter(1) = 0;
end

% Unpackage the muscle-wrapping parameters for Muscle 2 "TRImed"
u0 = wrapc(16); v0 = wrapc(17); up0 = wrapc(18); vp0 = wrapc(19); lg0 = wrapc(20);
p_l_2 = wrapc(21:23);
q_l_2 = wrapc(24:26);
uq0 = wrapc(27); vq0 = wrapc(28); upq0 = wrapc(29); vpq0 = wrapc(30);

if wrap_flag == 1
    % Broyden method
    [lmt(2),lmtd(2),u(2),v(2),up(2),vp(2),lg(2),p_l(:,2),q_l(:,2),uq(2),vq(2),upq(2),vpq(2),iter(2),event(2),Jm2] = SolveWrapping_BM(model,Xup,V,J,2,u0,v0,up0,vp0,lg0,p_l_2,q_l_2,uq0,vq0,upq0,vpq0);
    m2w_p = 0;
else
    % Obstcle-Set Method Analytical solution for the cylinder
    [lmt(2),lmtd(2),event(2),Jm2,m2w_p] = SolveWrapping_OS(model,Xup,V,J,2);
    % These variables are not used
    u(2) = 0; v(2) = 0; up(2) = 0; vp(2) = 0; lg(2) = 0;
    p_l(:,2) = zeros(3,1); q_l(:,2) = zeros(3,1); 
    uq(2) = 0; vq(2) = 0; upq(2) = 0; vpq(2) = 0; iter(2) = 0;
end

NMUS = length(model.m);

par_Tact = 0.015;    % constant muscle activation time
par_Tdeact = 0.060;  % constant muscle deactivation time

% Muscle parameters, details are provided below in MuscleDynamics function
mp = [987.26,0.0858,0.0535,0,10;             % For Muscle 1 "BRA"
      624.30,0.1138,0.0908,0.15707963,10];   % For Muscle 2 "TRImed"
  
% Compute the muscle activation and contraction dynamics
for i=1:NMUS
    % Mscle activation dynamics
    f_a = 0.5*tanh(0.1*(e(i)-a(i)));
    adot(i) = ((1/(par_Tact*(0.5+1.5*a(i))))*(f_a + 0.5) + ((0.5+1.5*a(i))/par_Tdeact)*(-f_a + 0.5))*(e(i)-a(i));
    % Muscle contraction dynamics
    [lnm(i),ft(i),fact(i),fpas(i),fv(i),FT(i),FM(i),dlM_dt(i)] = MuscleDynamics(lmt(i), lmtd(i), lnM(i), a(i), mp(i,:));
end

% tau_m = -[Jm1;Jm2]'*FM';
R = [Jm1;Jm2];
tau_m1 = -Jm1'*FM(1);
tau_m2 = -Jm2'*FM(2);
tau_m = tau_m1 + tau_m2;
tau = tau + tau_m;

% Compute the generalized accelerations
qdd = H \ (tau - C);

% Compute mechanical power
pwr = qd'*tau_m;
wrk = pwr;

% Group the output variable for integration or visualization
out = [qd;qdd;adot';dlM_dt';wrk;...   % Generalized accelerations, Generalized velocities, d/dt muscle activations, normalized fiber velocities, mechanical work
    u(1);v(1);up(1);vp(1);lg(1);p_l(:,1);q_l(:,1);uq(1);vq(1);upq(1);vpq(1);... % Muscle-wrapping algorithm parameters for Muscle 1 "BRA"
    u(2);v(2);up(2);vp(2);lg(2);p_l(:,2);q_l(:,2);uq(2);vq(2);upq(2);vpq(2);... % Muscle-wrapping algorithm parameters for Muscle 2 "TRImed"
    iter';event';...           % Iterations made by the muscle wrapping algorithm and the bool event inidcator of active wrapping
    ft';fact';fpas';fv';...    % Muscle force-length-velocity behavior
    lmt';lmtd';lnm';...        % Muscle-tendon length and velocity. Normalized fiber lengths
    FT';FM';...                % Tendon and muscle force
    tau_m1;tau_m2;p_tau;...    % Joint torques produced by Muscle 1 "BRA", Muscle 2 "TRImed" and passive torques
    fe;tau_ext;...             % External force at the upper limb's wrist and the resultant external joint torques
    pwr];                      % Mechanical power
end

function [lnM,ft,fact,fpas,fv,FT,FM,dlnM_dt] = MuscleDynamics(lmT, dlmT, lnM, a, mp)
% Maximum isometric force that the fibers can generate
FM_max = mp(1);
% Optimal length of the muscle fibers
lM_opt = mp(2);
% Resting length of the tendon
lT_slack = mp(3);
% Angle between tendon and fibers at optimal fiber length expressed in radians
alpha_opt = mp(4);
% Maximum contraction velocity of the fibers, in optimal
% fiberlengths/second
vM_max = mp(5);

% Muscle length denormalized
lM = lnM * lM_opt;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The muscle equilibrium equation as described in
% M. Millard, T. Uchida, A. Seth, and S. L. Delp, “Flexing computational 
% muscle: modeling and simulation of musculotendon dynamics,”
% Journal of biomechanical engineering, vol. 135, no. 2, p. 021005, 2013.
%
% This fragment of code is based on the OpenSim v4.0 C++ API
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = lM_opt * sin(alpha_opt);

% begin with small tendon force
temp_lt = lT_slack*1.01;
temp_lnt = temp_lt / lT_slack;

% Minimum fiber length
alpha_max = acos(0.1);
if alpha_max > eps
    lm_min = h/sin(alpha_max);
else
    lm_min = lM_opt*0.01;
end

% Minimum fiber length along tendon
ls_min = lm_min*cos(alpha_max);

% Fiber length
temp_ls = lmT - temp_lt;
if temp_ls >= ls_min
    temp_lm = sqrt(h*h + temp_ls*temp_ls);
else
    temp_lm = lm_min;
end

if temp_lm > lm_min
    temp_lm = temp_lm;
else
    temp_lm = lm_min;
end
% Normalized fiber length
temp_lnm = temp_lm / lM_opt;

% Muscle kinematics
[alpha, temp_lt, temp_lnm, temp_lnt] = calc_kinematics(h, alpha_max, lmT, temp_lm, lM_opt, lT_slack);

% Multipliers
[temp_fL, temp_fPE, temp_fSE] = calc_multipliers(temp_lnm, temp_lnt);

% Starting guess at the force-velocity multiplier is static
temp_fV = 1.0;
temp_dlnm  = 0.0;
temp_Fm = FM_max*(a*temp_fL*temp_fV + temp_fPE + 0.1*temp_dlnm);

[DFm_Dlm, DFs_Dlm, DFs_Dls, DFt_Dlt, DFt_Dlm] = calc_partials(FM_max, temp_lnm, lM_opt, a, temp_fV, h, temp_lm, alpha, temp_Fm, temp_lnt, lT_slack);

[temp_vt, temp_dlM, temp_dlnm, temp_fV] = calc_velocity(DFs_Dls, DFt_Dlt, temp_lnt, dlmT, alpha, vM_max, lM_opt);

[temp_Fm, temp_Fs, temp_Ft, temp_Ferr] = calc_forceErr(FM_max, a, temp_fL, temp_fV, temp_fPE, temp_dlnm, alpha, temp_fSE);

ferrPrev = temp_Ferr;
lcePrev = temp_lm;

iter = 1;
aSolTolerance = max(1e-8*FM_max, eps*10);
aMaxIterations = 200;

while( (abs(temp_Ferr) > aSolTolerance) && (iter < aMaxIterations))
        % Compute the search direction
        DFerr_Dlm = DFs_Dlm - DFt_Dlm;
        step = 1.0;
        while (abs(temp_Ferr) >= abs(ferrPrev))
            % Compute the Newton step
            delta_lce = -step*ferrPrev / DFerr_Dlm;
            % Take a Newton Step if the step is nonzero
            if (abs(delta_lce) > eps)
                temp_lm = lcePrev + delta_lce;
            else
                % We've stagnated or hit a limit; assume we are hitting local
                % minimum and attempt to approach from the other direction.
                temp_lm = lcePrev - sign(delta_lce)*sqrt(eps);
                % Force a break, which will update the derivatives of
                % the muscle force and estimate of the fiber-velocity
                step = 0;
            end
            
            if (temp_lm < lm_min)
                temp_lm = lm_min;
            end
            
            [alpha, temp_lt, temp_lnm, temp_lnt] = calc_kinematics(h, alpha_max, lmT, temp_lm, lM_opt, lT_slack);
            [temp_fL, temp_fPE, temp_fSE] = calc_multipliers(temp_lnm, temp_lnt);
            [temp_Fm, temp_Fs, temp_Ft, temp_Ferr] = calc_forceErr(FM_max, a, temp_fL, temp_fV, temp_fPE, temp_dlnm, alpha, temp_fSE);
            
            if (step <= sqrt(eps))
                break;
            else
                step = 0.5*step;
            end
        end
        
        ferrPrev = temp_Ferr;
        lcePrev = temp_lm;
        
        [DFm_Dlm, DFs_Dlm, DFs_Dls, DFt_Dlt, DFt_Dlm] = calc_partials(FM_max, temp_lnm, lM_opt, a, temp_fV, h, temp_lm, alpha, temp_Fm, temp_lnt, lT_slack);
        [temp_vt, temp_dlM, temp_dlnm, temp_fV] = calc_velocity(DFs_Dls, DFt_Dlt, temp_lnt, dlmT, alpha, vM_max, lM_opt);
        iter = iter + 1;
end

lnT = temp_lnt;
cos_alpha = cos(alpha);
lnM = temp_lnm;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The parameters of Muscle force-length-velocity and tendon force-length 
% curves are taken from 
%
% F. De Groote, A. L. Kinney, A. V. Rao, and B. J. Fregly, “Evaluation
% of direct collocation optimal control problem formulations for solving
% the muscle redundancy problem,” Annals of biomedical engineering,
% vol. 44, no. 10, pp. 2922–2936, 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tendon force-length curve
ft = tendonFSE(lnT);

% Active force-length curve
fact = activeFL(lnM);

% Passive force-length curve
fpas = passiveFPE(lnM);

%%
% Tendon force
FT = FM_max * ft;
% Muscle force
FM = FT/cos_alpha;
% Normalized Muscle force
FnM = FM/FM_max;
% Inverse velocity curve
FvM = (FnM - fpas)/(a*fact);
vnM = invforceVel(FvM);
% normalized fiber velocity
dlnM_dt = ((vM_max*lM_opt)/lM_opt)*vnM;

% Force-velocity curve
fv = forceVel(vnM);

end

function fact = activeFL(lnM)
% Muscle active force-length curve
fact = 0;

b = [0.8150671134243542, 0.433004984392647, 0.1;
    1.055033428970575, 0.716775413397760, 1.0;
    0.162384573599574, -0.029947116970696, 0.353553390593274;
    0.063303448465465, 0.200356847296188, 0.0];
 
for i=1:3
    fact = fact + b(1,i).*exp(-(0.5.*(lnM-b(2,i)).^2)./((b(3,i)+b(4,i).*lnM)).^2);
end

end

function fpas = passiveFPE(lnM)
% Muscle passive force-length curve
kpe = 4.0;
e0 = 0.6;
lnm_min = 0.2;

offset = exp((kpe*(lnm_min-1))/e0);

denom = exp(kpe)-offset;

fpas = (exp((kpe*(lnM-1))/e0) - offset)/denom;
end

function ft = tendonFSE(lnT)
% Tendon force-length curve
c1 = 0.200;
c2 = 1.0;
c3 = 0.200;

e0T = 0.049;

kT = (log((1.0+c3)/c1))/(1.0+e0T-c2);

ft = c1*exp(kT*(lnT-c2))-c3;
end

function fv = forceVel(vnM)
% Muscle force-velocity curve
d1 = -0.3211346127989808;
d2 = -8.149;
d3 = -0.374;
d4 = 0.8825327733249912;

tempV = d2 .* vnM + d3;
tempLogArg = tempV + sqrt((tempV).^2 + 1.0);

fv = d1 * log(tempLogArg) + d4;
end

function vnM = invforceVel(FvM)
% Muscle inverse force-velocity curve
d1 = -0.3211346127989808;
d2 = -8.149;
d3 = -0.374;
d4 = 0.8825327733249912;

vnM = -(1/d2)*(sinh((-FvM+d4)/d1) + d3);
end

function [alpha, lt, lnm, lnt] = calc_kinematics(h, alpha_max, lmt, lm, lm_opt, lt_slack)
% Pennation angle
if (h/lm) < sin(alpha_max)
    alpha = asin(h/lm);
else
    alpha = alpha_max;
end

lt = lmt - lm*cos(alpha);
lnm = lm/lm_opt;
lnt = lt/lt_slack;

end

function [fL, fPE, fSE] = calc_multipliers(lnm, lnt)
fL = activeFL(lnm);
fPE = passiveFPE(lnm);
fSE = tendonFSE(lnt);
end

function [DFm_Dlm, DFs_Dlm, DFs_Dls, DFt_Dlt, DFt_Dlm] = calc_partials(FM_max, temp_lnm, lM_opt, a, temp_fV, h, temp_lm, alpha, temp_Fm, temp_lnt, lT_slack)
% (\partial Fm) / (\partial lm)
DFm_Dlm = calc_DFm_Dlm(FM_max, temp_lnm, lM_opt, a, temp_fV);
% (\partial Fs) / (\partial lm)
DFs_Dlm = calc_DFs_Dlm(h, temp_lm, alpha, temp_Fm, DFm_Dlm);
% (\partial Fs) / (\partial ls)
DFs_Dls = calc_DFs_Dls(h, temp_lm, alpha, DFs_Dlm);
% (\partial Ft) / (\partial lt)
DFt_Dlt = calc_DFt_Dlt(FM_max, temp_lnt, lT_slack);
% (\partial Ft) / (\partial lm)
DFt_Dlm = calc_DFt_Dlm(h, temp_lm, alpha, DFt_Dlt);
end

function [vt, dlM, dlnm, fV] = calc_velocity(DFs_Dls, DFt_Dlt, lnt, dlmT, alpha, vM_max, lM_opt)
if ((abs(DFs_Dls + DFt_Dlt) > eps) && (lnt > 1.0))
    vt = DFs_Dls / (DFs_Dls + DFt_Dlt) * dlmT;
else
    vt = dlmT;
end

dlM = (dlmT - vt)*cos(alpha);
dlnm = dlM / (vM_max*lM_opt);
fV = forceVel(dlnm);
end

function [Fm, Fs, Ft, Ferr] = calc_forceErr(FM_max, a, fL, fV, fPE, dlnm, alpha, fSE)
Fm = FM_max*(a*fL*fV + fPE + 0.1*dlnm);
Fs = Fm*cos(alpha);
Ft = FM_max*fSE;
Ferr = Fs - Ft;
end

function out = calc_DFm_Dlm(F0M, lnM, lm_opt, a, fV)

% (\partial FL) / (\partial lm)
DFL_Dlm = calc_DFL_Dlm(lnM, lm_opt);
% (\partial FPE) / (\partial lm)
DFPE_Dlm = calc_DFPE_Dlm(lnM, lm_opt);
% (\partial FPE) / (\partial lm)
out = F0M*(a*DFL_Dlm*fV + DFPE_Dlm);

end

function out = calc_DFs_Dlm(h, lm, alpha, Fm, DFm_Dlm)
    Dalpha_Dlm = calc_Dalpha_Dlm(h, lm);
    Dcosalpha_Dlm = -sin(alpha) * Dalpha_Dlm;
    out = DFm_Dlm * cos(alpha) + Fm * Dcosalpha_Dlm;
end

function out = calc_DFs_Dls(h, lm, alpha, DFs_Dlm)
    Dalpha_Dlm = calc_Dalpha_Dlm(h, lm);
    Dls_Dlm = cos(alpha) - lm*sin(alpha)*Dalpha_Dlm;
    out = DFs_Dlm * (1/Dls_Dlm);
end

function out = calc_DFt_Dlt(F0M, lnT, lt_slack)
c1 = 0.200;
c2 = 1.0;
c3 = 0.200;
e0T = 0.049;

kT = (log((1.0+c3)/c1))/(1.0+e0T-c2);
% (\partial FSE) / (\partial lnt)
DFSE_Dlnt = kT*c1*exp(kT*(lnT-c2));

% (\partial lnt) / (\partial lt)
Dlnt_Dlt = 1/lt_slack;

% (\partial Ft) / (\partial lt)
out = F0M * DFSE_Dlnt * Dlnt_Dlt;

end

function out = calc_DFt_Dlm(h, lm, alpha, DFt_Dlt)
% (\partial lt) / (\partial lm)
Dlt_Dlm = calc_Dlt_Dlm(h, lm, alpha);
% (\partial Ft) / (\partial lm)
out = DFt_Dlt * Dlt_Dlm;
end

function out = calc_Dlt_Dlm(h, lm, alpha)
Dalpha_Dlm = calc_Dalpha_Dlm(h, lm);
out = lm*sin(alpha)*Dalpha_Dlm - cos(alpha);
end

function out = calc_Dalpha_Dlm(h, lm)
    out = - h/((lm^2) * sqrt(1 - (h / lm)^2));
end

function out = calc_DFL_Dlm(lnM, lm_opt)
% (\partial lnm) / (\partial lm)
Dlmn_Dlm = 1/lm_opt;
% (\partial FL) / (\partial lnm)
DFL_Dlnm = 0;

b = [0.8150671134243542, 0.433004984392647, 0.1;
    1.055033428970575, 0.716775413397760, 1.0;
    0.162384573599574, -0.029947116970696, 0.353553390593274;
    0.063303448465465, 0.200356847296188, 0.0];
 
for i=1:3
    B = ((b(2,i)*b(4,i) + b(3,i))*(b(2,i) - lnM))/(b(3,i) + b(4,i).*lnM)^3;
    DFL_Dlnm = DFL_Dlnm + b(1,i).*exp(-(0.5.*(lnM-b(2,i)).^2)./((b(3,i)+b(4,i).*lnM)).^2)*B;
end
% (\partial FL) / (\partial lm)
out = DFL_Dlnm * Dlmn_Dlm;
end

function out = calc_DFPE_Dlm(lnM, lm_opt)
% (\partial lnm) / (\partial lm)
Dlmn_Dlm = 1/lm_opt;
% (\partial FPE) / (\partial lnm)
kpe = 4.0;
e0 = 0.6;
lnm_min = 0.2;

offset = exp((kpe*(lnm_min-1))/e0);

denom = exp(kpe)-offset;

DFPE_Dlnm = (kpe*exp((kpe*(lnM-1))/e0))/(e0*denom);
% (\partial FPE) / (\partial lm)
out = DFPE_Dlnm * Dlmn_Dlm;
end