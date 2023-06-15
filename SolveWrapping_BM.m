function [lmt,Ld,u,v,up,vp,lg,p_l,q_l,uq,vq,upq,vpq,iter,event,Jm] = SolveWrapping_BM(model,Xup,V,J,midx,u0,v0,up0,vp0,lg0,p_l,q_l,uq0,vq0,upq0,vpq0)
% Solve the muscle wrapping
% model struct variable which contains the information of the whole
%       biomechanical system
% Xup   6D Kinematic Transforms of the system
% midx  muscle index
% u0    initial condition for the u parameter of the parametric equations
% v0    initial condition for the v parameter of the parametric equations
% lg0   initial condition for the geodesic length

u = u0; v = v0; lg = lg0;
up = up0; vp = vp0;

uq = uq0; vq = vq0; upq = upq0; vpq = vpq0;

widx = model.m{midx}.wgeom;

% vector containing initial conditions
ic = [u;v;up;vp;lg];

% Path points expressed in the wrap geometry frame coordinates
% Transform only via points in wrapping range
for i=1:size(model.m{midx}.P,2)
    if model.m{midx}.P{1,i}.parent == model.appearance.wrapCyl{widx}.parent
        R = model.appearance.wrapCyl{widx}.R;   % Rotation from parent to wrapping geometry
        p = model.appearance.wrapCyl{widx}.tra; % Translation from parent to wrapping geometry
    else
        T = inv(pluho(Xup{model.m{midx}.P{1,i}.parent})); % Transformation from path point's parent to its parent
        R = T(1:3,1:3)'*model.appearance.wrapCyl{widx}.R;      % Rotation from parent to wrapping geometry
        p = T(1:3,1:3)'*(model.appearance.wrapCyl{widx}.tra - T(1:3,4)); % Translation from parent to wrapping geometry
    end
    wp(:,i) = R'*(model.m{midx}.P{1,i}.location - p);
end

% Rotate the path points to mantain the cylinder in 0,0,0 coordinates and
% the heigth aligned with $z$ axis, so the points will be expressed in that
% coordinate frame. The scaling parameter is used to avoid numerical
% instability during the Broyden's iterations
scaling = 100;
p0 = (rx(pi/2)'*wp(:,model.m{midx}.wrange(1)))'*scaling;
q0 = (rx(pi/2)'*wp(:,model.m{midx}.wrange(2)))'*scaling;

% Check if points are in the same plane quadrant to avoid the computation 
% of a wrapping solution
nowrap = sign(p0(1))==sign(q0(1)) && sign(p0(2))==sign(q0(2)) && sign(p0(3))==sign(q0(3));

% initial points on the wrap geometry (cylinder)
r = model.appearance.wrapCyl{widx}.radius*scaling;
h = model.appearance.wrapCyl{widx}.height*scaling;
p = p_l'*scaling;
q = q_l'*scaling;

% integrator tolerances
geod_rel_tol = 1.e-6;
geod_abs_tol = 1.e-12;
options = odeset('RelTol',geod_rel_tol, 'AbsTol', geod_abs_tol);

%% Muscle Wrapping Algorithm
iter = 1;

if nowrap==0
    % Parameters for the Broyden method an line search algorithm
    rho = 0.9;
    sigma1 = 0.001;
    sigma2 = 0.001;
    beta = 0.5;

    alpha = 0.1;
    theta_u = 1;
    theta_l = (1-alpha^(1/3))/(1+alpha^(1/3));
    theta = theta_u;

    % initial wrapping point
    e = p-p0;
    e = e/norm(e);
    t = model.appearance.wrapCyl{widx}.F(u,v,up,vp,r);

    % final wrapping point
    e2 = q-q0;
    e2 = e2/norm(e2);
    t2 = model.appearance.wrapCyl{widx}.F(uq,vq,upq,vpq,r);

    f = @myfun;
    F = [(e-t)';(e2+t2)'];
    nF(1) = norm(F);
    nf = nF(1);
    % define H0
    B = jacobi(f,F,ic,r,p0,q0,model,widx);

    while nF(iter)>1e-10

        B_inv = pinv(B);
        d = -B_inv*F;

        et = eta(iter); % compute \eta_{k}
        % If ||F(x_k + p_k)|| <= \rho ||F(x_k)|| - \sigma_{2} ||p_k||^2,
        % then let \lambda_k = 1
        if norm(f(ic+d,r,p0,q0,model,widx)) <= (rho*nF(iter) - sigma2*norm(d)^2)
            lambda = 1;
        else 
            % Let \lambda_k be determined by Procedure 1  
            % (Approximate Norm Descent Line Search)
            [lambda, ~] = linesearch(f, nF(iter), ic, d, sigma1, beta, et, r, p0, q0, model, widx);
        end

        ic = ic + lambda*d;
        s = d;

        F = f(ic,r,p0,q0,model,widx);
        nF(iter+1) = norm(F);

        B = B + theta*(F*s'/(s'*s));
        iter = iter+1;  

    end

    [~,geodesic]=ode45(@int_gp,[0,ic(5)],ic(1:4),options,r);
    up = geodesic(1,1); vp = geodesic(1,2); 
    uq = geodesic(end,1); vq = geodesic(end,2);
    upq = geodesic(end,3); vpq = geodesic(end,4);
    p = [r*cos(up), r*sin(up), vp];
    q = [r*cos(uq), r*sin(uq), vq];

end

p0 = p0/scaling; q0 = q0/scaling;
p = p/scaling; q = q/scaling;

% This condition is the same as the one presented in
% B. A. Garner and M. G. Pandy, “The obstacle-set method for representing 
% muscle paths in musculoskeletal models,” Computer Methods in Biomechanics 
% and Biomedical Engineering, vol. 3, no. 1, pp. 1–30, 2000.
dt1 = model.m{midx}.wd*(r/scaling)*(p(1)*q(2)-p(2)*q(1)); 

if dt1>0 % if dt1>0 wrapping does occur
    
    % Compute the muscle-tendon length
    lmt = 0; j = 1;
    for i=1:size(model.m{midx}.P,2)
        if model.m{midx}.wrange(1)==i
            lmt = lmt + norm(rx(-pi/2)'*p'-wp(:,j));
            j = j+1;
        elseif model.m{midx}.wrange(2)==i
            lmt = lmt + norm(wp(:,j)-rx(-pi/2)'*q');
        else
            lmt = lmt + norm(wp(:,j+1)-wp(:,j));
            j = j+1;
        end
    end
    
    la = ic(5)/scaling;  % geodesic length
    lmt = lmt + la;
    
    % Does wrapping occurs? 1: yes
    event = 1;
    
    % p in shoulder frame (parent frame of the wrapping geometry)
    p_s = model.appearance.wrapCyl{widx}.tra + model.appearance.wrapCyl{widx}.R*(rx(-pi/2)'*(p'));
    % q in shoulder frame (parent frame of the wrapping geometry)
    q_s = model.appearance.wrapCyl{widx}.tra + model.appearance.wrapCyl{widx}.R*(rx(-pi/2)'*(q'));

    wp = [model.m{midx}.P{1,model.m{midx}.wrange(1)}.location,...
          model.m{midx}.P{1,model.m{midx}.wrange(2)}.location,...
          p_s, q_s];
    wg = model.appearance.wrapCyl{widx}.tra;
    % wp = [ps0,q0_s,p_s,q_s]; % wrapping points are concatenated in columns, fixed points first, contact points last
    % compute wrapping points in world coordinates
    [wp_w, e, vp_w, Jm] = wrapPosVel(model,midx,widx,Xup,wp,wg,V,J);

    % lmt velocity
    lmt_d = dot(e{1},vp_w{3}) - dot(e{1},vp_w{1}) + dot(e{2},vp_w{2}) - dot(e{2},vp_w{3});
    
else 
    % Compute the muscle-tendon length
    lmt = 0;
    for i=1:size(model.m{midx}.P,2)-1
        lmt = lmt + norm(wp(:,i+1)-wp(:,i));
    end
    
    % Does wrapping occurs? 0: no
    event = 0;
    
    wp = [model.m{midx}.P{1,model.m{midx}.wrange(1)}.location,...
          model.m{midx}.P{1,model.m{midx}.wrange(2)}.location];
    
      % compute via points in world coordinates
    [wp_w, e, vp_w, Jm] = nowrapPosVel(model,midx,Xup,wp,V,J);
    
    % lmt velocity
    lmt_d = dot(e{1},vp_w{2}-vp_w{1});
    
end

Ld = lmt_d;

u = ic(1); v = ic(2); up = ic(3); vp = ic(4); lg = ic(5);
p_l = p; q_l = q;

end

function J = jacobi(f,y0,x,r,p0,q0,model,widx)
% Numerical Jacobian for function f at x
% y0: f(x);

    delta = 1e-6*(max(1,sqrt(norm(x))));
%     delta = 1e-10;
    n = length(y0);
    m = length(x);
    J = zeros(n,m);
    for i=1:m
        dx = zeros(m,1);
        dx(i) = delta/2;
        J(:,i) = (f(x+dx,r,p0,q0,model,widx) - f(x-dx,r,p0,q0,model,widx))/delta;
    end
end

function [lambda,lsc] = linesearch(f, normf, x, p, sigma1, beta, et, Radius, p0, q0, model, widx)
% Approximate Norm Descent Line Search
lambda = 1;
i = 1;
lsc = 0;

while norm(f(x+lambda*p,Radius,p0,q0,model,widx)) > ((1 + et)*normf - sigma1*norm(lambda*p)^2)
   lambda = beta*lambda;
   i = i+1;
   lsc = lsc + 1;
end

end

function out = eta(k)
out = 1/(k+1)^2;
end

function ydot = int_gp(~,y,r)
u =  y(1);
v =  y(2);

ru = [-r*sin(u);r*cos(u);0];
rv = [0;0;1];

E = dot(ru,ru); Eu = 2*r^2*cos(u)*sin(u) - 2*r^2*cos(u)*sin(u); Ev = 0;
F = dot(ru,rv); Fu = 0; Fv = 0;
G = dot(rv,rv); Gu = 0; Gv = 0;

C111 = (G*Eu - 2*F*Fu + F*Ev)/(2*(E*G-F^2));
C121 = (G*Ev - F*Gu)/(2*(E*G-F^2));
C221 = (2*G*Fv - G*Gu + F*Gv)/(2*(E*G-F^2));

C112 = (2*E*Fu - E*Ev + F*Eu)/(2*(E*G-F^2));
C122 = (E*Gu - F*Ev)/(2*(E*G-F^2));
C222 = (E*Gv - 2*F*Fv + F*Gu)/(2*(E*G-F^2));

sq = norm(ru*y(3) + rv*y(4));
if sq > 1+eps
    t = ru*y(3) + rv*y(4); 
    par = [y(3);y(4)]./norm(t);
    par = par .* (1/(E*par(1)^2+2*F*par(1)*par(2)+G*par(2)^2));

    y(3) = par(1);
    y(4) = par(2);
end

ydot = zeros(4,1);
ydot(1) = y(3); %p
ydot(2) = y(4); %q
ydot(3) = -C111*ydot(1)^2 - 2*C121*ydot(1)*ydot(2) - C221*ydot(2)^2;
ydot(4) = -C112*ydot(1)^2 - 2*C122*ydot(1)*ydot(2) - C222*ydot(2)^2;
end

function F = myfun(ic,r,p0,q0,model,widx)
geod_rel_tol = 1.e-6;
geod_abs_tol = 1.e-12;
options = odeset('RelTol',geod_rel_tol, 'AbsTol', geod_abs_tol);

u = ic(1); v = ic(2); up = ic(3); vp = ic(4); lg = ic(5);

% initial wrapping point
p = [r*cos(u), r*sin(u), v];
e = p-p0;
e = e/norm(e);
t = model.appearance.wrapCyl{widx}.F(u,v,up,vp,r);

[~,geodesic]=ode45(@int_gp,[0,lg],ic(1:4),options,r);

% final wrapping point
q = [r*cos(geodesic(end,1)), r*sin(geodesic(end,1)), geodesic(end,2)];
e2 = q-q0;
e2 = e2/norm(e2);
t2 = model.appearance.wrapCyl{widx}.F(geodesic(end,1),geodesic(end,2),geodesic(end,3),geodesic(end,4),r);

F = [(e-t)';(e2+t2)'];
end