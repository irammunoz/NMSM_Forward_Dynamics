function  out = FDcrb( model, y, tau, fe, f_ext)

% FDcrb  Forward Dynamics via Composite-Rigid-Body Algorithm
% FDcrb(model,q,qd,tau,f_ext)  calculates the forward dynamics of a
% kinematic tree via the composite-rigid-body algorithm.  q, qd and tau are
% vectors of joint position, velocity and force variables; and the return
% value is a vector of joint acceleration variables.  f_ext is an optional
% argument specifying the external forces acting on the bodies.  It can be
% omitted if there are no external forces.  The format of f_ext is
% explained in the source code of apply_external_forces.
NDOF = model.NB;
q = y(1:NDOF); qd = y(NDOF+1:NDOF*2);

if nargin == 4
  [H,C,Xup,Xa,V,J] = HandC( model, q, qd );
else
  [H,C,Xup,Xa,V,J] = HandC( model, q, qd, f_ext );
end

tau = tau - 3*eye(model.NB)*qd;

if exist('fe','var')
    Xee = model.KXtree{1} * Xa{model.kparent};
    Jl_e = [eye(3) zeros(3,3);model.KXtree{1}(4:6,1:3) eye(3)]*J{end};
    
    % transform generalized forces to a wrench expressed in absolute
    % coordinates
    fh = pinv(Jl_e')*tau;
    fh = Xee'*fh;
    
end

qdd = H \ (tau - C);

out = [qd;qdd;fh];

end

