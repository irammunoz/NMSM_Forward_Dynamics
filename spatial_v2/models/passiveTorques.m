function tau = passiveTorques(model,q,qd)
% passive joint moment parameters
par_JointK1			= 1;        % overall joint stiffness (Nm/rad)
par_JointK2			= 10000;	% stiffness at joint limits (Nm/rad^2)
par_JointB			= 1;		% joint damping (Nms/rad), exists always
par_MinAngle(1) 	=  model.coordrange{1}(1);
par_MaxAngle(1) 	=  model.coordrange{1}(2);
par_MinAngle(2) 	=  model.coordrange{2}(1);
par_MaxAngle(2) 	=  model.coordrange{2}(2);

NDOF = model.NB;
tau = zeros(NDOF,1);

for i=1:NDOF
    % Is angle above upper limit of ROM?
    d = q(i) - par_MaxAngle(i);
    if d > 0.0
        tau(i) = -par_JointK2 * d*d;
    end

    % Is angle below lower limit of ROM?
    d = q(i) - par_MinAngle(i);
    if d < 0.0
        tau(i) = par_JointK2 * d*d;
    end

    % Add a small amount of damping and overall stiffness
    tau(i) = tau(i) - par_JointB * qd(i) - par_JointK1 * q(i);

end

end