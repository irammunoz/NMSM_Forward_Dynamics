function  model = arm22_model

% Definition of degrees of freedom
model.NB = 2;
model.parent = [0,1];

[model.child, model.sister] = FindChildren(model);

% Rotation matrix to transform the reference frame of the "arm26.osim" 
% model's parameters 
rx90 = rx(pi/2)';

% Joints' axes of rotation
model.ar{1} = rx90*[-0.05889802;0.0023;0.99826136];
model.ar{2} = rx90*[0.04940001;0.03660001;0.99810825];

% Joint's limit ranges
model.coordrange{1} = [-90*pi/180 180*pi/180];
model.coordrange{2} = [0*pi/180 130*pi/180];

% Definition of Joints' motion subspace matrices
for i=1:model.NB
    model.L{i} = zeros(6,model.NB);
    [~, S] = Rodrigues(model.ar{i},0);
    model.L{i}(:,i) = S;
end

% base
offset = rx90*[0;0.8;0];                                   % location_in_parent
% arm_r_humerus
r_humerus_p = rx90*[-0.017545;-0.007;0.17];                % location_in_parent
r_humerus_g = r_humerus_p + offset;                        % location_in_global
% arm_r_ulna
r_ulna_radius_hand_p = rx90*[0.0061;-0.2904;-0.0123];      % location_in_parent
r_ulna_radius_hand_g = r_ulna_radius_hand_p + r_humerus_g; % location_in_global

% Joints' spatial coordinates
model.Xtree{1} = xlt(r_humerus_g);
model.Xtree{2} = xlt(r_ulna_radius_hand_p);

% Inertial parameters
model.I{1} = mcI( 1.864572, rx90*[0; -0.180496; 0], rx90'*diag([0.01481,0.004551,0.013193])*rx90 ); 
model.I{2} = mcI( 1.534315, rx90*[0; -0.181479; 0], rx90'*diag([0.019281,0.001571,0.020062])*rx90 );

% Gravity vector
model.gravity = [0; 0; -9.81];

% model appearance
model.appearance.body{1} = stlread('ground_ribs_rot.stl');
model.appearance.body{2} = stlread('ground_spine_rot.stl');
model.appearance.body{3} = stlread('ground_skull_rot.stl');
model.appearance.body{4} = stlread('ground_jaw_rot.stl');
model.appearance.body{5} = stlread('ground_r_clavicle_rot.stl');
model.appearance.body{6} = stlread('ground_r_scapula_rot.stl');
model.appearance.body{7} = stlread('arm_r_humerus_rot.stl');
model.appearance.body{8} = stlread('arm_r_ulna_rot.stl');
model.appearance.body{9} = stlread('arm_r_radius_rot.stl');
model.appearance.body{10} = stlread('arm_r_lunate_rot.stl');
model.appearance.body{11} = stlread('arm_r_scaphoid_rot.stl');
model.appearance.body{12} = stlread('arm_r_pisiform_rot.stl');
model.appearance.body{13} = stlread('arm_r_triquetrum_rot.stl');
model.appearance.body{14} = stlread('arm_r_capitate_rot.stl');
model.appearance.body{15} = stlread('arm_r_trapezium_rot.stl');
model.appearance.body{16} = stlread('arm_r_trapezoid_rot.stl');
model.appearance.body{17} = stlread('arm_r_hamate_rot.stl');
model.appearance.body{18} = stlread('arm_r_1mc_rot.stl');
model.appearance.body{19} = stlread('arm_r_2mc_rot.stl');
model.appearance.body{20} = stlread('arm_r_3mc_rot.stl');
model.appearance.body{21} = stlread('arm_r_4mc_rot.stl');
model.appearance.body{22} = stlread('arm_r_5mc_rot.stl');
model.appearance.body{23} = stlread('arm_r_thumbprox_rot.stl');
model.appearance.body{24} = stlread('arm_r_thumbdist_rot.stl');
model.appearance.body{25} = stlread('arm_r_2proxph_rot.stl');
model.appearance.body{26} = stlread('arm_r_2midph_rot.stl');
model.appearance.body{27} = stlread('arm_r_2distph_rot.stl');
model.appearance.body{28} = stlread('arm_r_3proxph_rot.stl');
model.appearance.body{29} = stlread('arm_r_3midph_rot.stl');
model.appearance.body{30} = stlread('arm_r_3distph_rot.stl');
model.appearance.body{31} = stlread('arm_r_4proxph_rot.stl');
model.appearance.body{32} = stlread('arm_r_4midph_rot.stl');
model.appearance.body{33} = stlread('arm_r_4distph_rot.stl');
model.appearance.body{34} = stlread('arm_r_5proxph_rot.stl');
model.appearance.body{35} = stlread('arm_r_5midph_rot.stl');
model.appearance.body{36} = stlread('arm_r_5distph_rot.stl');

% Define parent of each geometry
%   ground (1) bodies = 6
%   humerus (2) bodies = 1
%   forearm & hand (3) bodies = 29
model.b_index = [ones(6,1)*0;ones(1,1)*1;ones(29,1)*2];

model.tr = [offset,r_humerus_g,r_ulna_radius_hand_p];

%% Define the end-effector transformation
model.kparent = [2];
wrist = rx90*[0.005;-0.265195;0.06];
model.KXtree{1} = xlt(wrist);

%% Define the wrapping surfaces
R = eulerXYZ(0,0,0);
rf = [0;0;0];

% WrapCylinder name="TRI"
model.appearance.wrapCyl{1}.radius = 0.016;
model.appearance.wrapCyl{1}.height = 0.05;
[model.appearance.wrapCyl{1}.vertices, model.appearance.wrapCyl{1}.sideFaces, model.appearance.wrapCyl{1}.bottomFaces] = calcCylinder(rf, R, model.appearance.wrapCyl{1}.radius, model.appearance.wrapCyl{1}.height, 21);
model.appearance.wrapCyl{1}.parent = 1; % r_humerus
model.appearance.wrapCyl{1}.tra = rx90*[0.0028;-0.2919;-0.0069];
model.appearance.wrapCyl{1}.rot = rx90*[-0.14015;-0.00628319;0.154985];
model.appearance.wrapCyl{1}.R = eulerXYZ(model.appearance.wrapCyl{1}.rot(1), model.appearance.wrapCyl{1}.rot(2), model.appearance.wrapCyl{1}.rot(3));
model.appearance.wrapCyl{1}.F = @(u,v,up,vp,r)Cyl_tangent(u,v,up,vp,r);
model.appearance.wrapCyl{1}.G = @(t,y,r)Cyl_geodesic(t,y,r);

%% Muscle Via Points
% PathPoint name="BRA-P1"
model.m{1}.P{1}.location = rx90*[0.0068;-0.1739;-0.0036];
model.m{1}.P{1}.parent = 1;
% PathPoint name="BRA-P2"
model.m{1}.P{2}.location = rx90*[-0.0032;-0.0239;0.0009];
model.m{1}.P{2}.parent = 2;

model.m{1}.wrange = [1 2]; % range of via points where wrapping occurs
model.m{1}.wgeom = 1;      % index of wrapping surface
model.m{1}.wd = 1;         % wrapping direction
model.m{1}.name = 'BRA';   % muscle name

% PathPoint name="TRImed-P1"
model.m{2}.P{1}.location = rx90*[-0.00838;-0.13695;-0.00906];
model.m{2}.P{1}.parent = 1;
% PathPoint name="TRImed-P2"
model.m{2}.P{2}.location = rx90*[-0.02601;-0.15139;-0.0108];
model.m{2}.P{2}.parent = 1;
% PathPoint name="TRImed-P3"
model.m{2}.P{3}.location = rx90*[-0.03184;-0.22637;-0.01217];
model.m{2}.P{3}.parent = 1;
% PathPoint name="TRImed-P4"
model.m{2}.P{4}.location = rx90*[-0.01743;-0.26757;-0.01208];
model.m{2}.P{4}.parent = 1;
% PathPoint name="TRImed-P5"
model.m{2}.P{5}.location = rx90*[-0.0219;0.01046;-0.00078];
model.m{2}.P{5}.parent = 2;

model.m{2}.wrange = [4 5];    % range of via points where wrapping occurs
model.m{2}.wgeom = 1;         % index of wrapping surface
model.m{2}.wd = -1;           % wrapping direction
model.m{2}.name = 'TRImed';   % muscle name

end
