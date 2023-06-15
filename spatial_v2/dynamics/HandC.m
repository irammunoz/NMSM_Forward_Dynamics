function  [H,C,Xup,Xa,v,J] = HandC( model, q, qd, f_ext )

% HandC  Calculate coefficients of equation of motion
% [H,C]=HandC(model,q,qd,f_ext)  calculates the coefficients of the
% joint-space equation of motion, tau=H(q)qdd+C(d,qd,f_ext), where q, qd
% and qdd are the joint position, velocity and acceleration vectors, H is
% the joint-space inertia matrix, C is the vector of gravity,
% external-force and velocity-product terms, and tau is the joint force
% vector.  Algorithm: recursive Newton-Euler for C, and
% Composite-Rigid-Body for H.  f_ext is an optional argument specifying the
% external forces acting on the bodies.  It can be omitted if there are no
% external forces.  The format of f_ext is explained in the source code of
% apply_external_forces.

a_grav = get_gravity(model);
J = {zeros(6,model.NB)};
Jd = {zeros(6,model.NB)};

for i = 1:model.NB
  [ XJ, S{i} ] = Rodrigues(model.ar{i},q(i));
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    avp{i} = Xup{i} * -a_grav;
	Xa{i} = Xup{i}; % Compute i^X_0:  Plücker transform 0 relative to i, expressed in i coordinates
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    avp{i} = Xup{i}*avp{model.parent(i)} + crm(v{i})*vJ;
	Xa{i} = Xup{i} * Xa{model.parent(i)}; % Compute i^X_0:  Plücker transform 0 relative to i, expressed in i coordinates
  end
  fvp{i} = model.I{i}*avp{i} + crf(v{i})*model.I{i}*v{i};
  % Local Geometric Jacobian
  J{i+1} = Xup{i} * J{model.parent(i)+1} + model.L{i};
  % Time derivative of Local Geometric Jacobian
  Jd{i+1} = Xup{i} * Jd{model.parent(i)+1} - crm(vJ)*J{i+1};
end

if nargin == 4
  fvp = apply_external_forces( model.parent, Xup, fvp, f_ext );
end

for i = model.NB:-1:1
  C(i,1) = S{i}' * fvp{i};
  if model.parent(i) ~= 0
    fvp{model.parent(i)} = fvp{model.parent(i)} + Xup{i}'*fvp{i};
  end
end

IC = model.I;				% composite inertia calculation

for i = model.NB:-1:1
  if model.parent(i) ~= 0
    IC{model.parent(i)} = IC{model.parent(i)} + Xup{i}'*IC{i}*Xup{i};
  end
end

H = zeros(model.NB);

for i = 1:model.NB
  fh = IC{i} * S{i};
  H(i,i) = S{i}' * fh;
  j = i;
  while model.parent(j) > 0
    fh = Xup{j}' * fh;
    j = model.parent(j);
    H(i,j) = S{j}' * fh;
    H(j,i) = H(i,j);
  end
end
