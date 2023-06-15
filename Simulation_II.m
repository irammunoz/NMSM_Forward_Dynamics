%% Simulation for Condition II 

% Add utility functions to run the simulation
addpath(genpath('spatial_v2'));
addpath(genpath('VTP'));
addpath(genpath('OpenSim_Sim'));

% Load model
model = arm22_model;

ti = 0;    % initial time
tf = 5;    % final time
ts = 0.01; % step size

%% Set the initial conditions
qd_0 = [0; 0];           % initial generalized velocities
q_0 = [0; 0.349066];     % initial generalized positions


% For the muscle-wrapping algorithm fifteen parameters are define per 
% muscle to avoid recomputation of a solution if it is the same as the
% previous one.
% First five elements of correspond to the initial wrapping 
% conditions in eq(18)
%
%   $u_{P^{i}}$, $v_{P^{i}}$, $u'_{P^{i}}$, $v'_{P^{i}}$, $s^{i}$
%
% The next six elements correspond to the initial positions of wrapping
% points (these are computed from the surface parametric equations)
%
%   $P^{i}$ and $Q^{i}$
%
% Finally last four elements correspond to initial wrapping conditions
% similar to eq(8) but for poitnt $Q^{i}$
%
%   $u_{Q^{i}}$, $v_{Q^{i}}$, $u'_{Q^{i}}$, $v'_{Q^{i}}$
%

wrap_0 = [ 0; 0; 0.7071;0.7071;eps; 1.60;0;0; 1.60;0;0; 0; 0; 0.7071;0.7071;
          pi;pi;-0.7071;0.7071;eps;-1.60;0;0;-1.60;0;0;pi;pi;-0.7071;0.7071];

% Initial muscle activations
a_0 = [0.5; 0.001]; 

% Initial normalized fiber lengths
% These are computed from the muscle equilibrium equation as described in
% M. Millard, T. Uchida, A. Seth, and S. L. Delp, “Flexing computational 
% muscle: modeling and simulation of musculotendon dynamics,”
% Journal of biomechanical engineering, vol. 135, no. 2, p. 021005, 2013.
lnM_0 = [1.013501012900025; 0.646461437825719];

% Initial mechanical work
wrk_0 = 0;

% All initial conditions are grouped here
y0 = [q_0; qd_0; a_0; lnM_0; wrk_0];

%% Compute the time-varying external wrench
y0_r = [-0.75*pi;pi/2;0;0];
y_R = ode4_rob(rr_robot_a22,@FDcrb,@controller,ti,ts,tf,y0_r);

%% Integrate the system
wrap_flag = 1; % 1: Solve the system using the SGW method
y_BM = ode4_NMSM_fe(model,@FDcrb_NMSM_fe,@excitation_signal,wrap_0,ti,ts,tf,y0,wrap_flag,y_R(9:14,:));

%%
Plot_Results_II;