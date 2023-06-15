function [fp_w, e, v, Jm] = nowrapPosVel(model,midx,Xup,wp,V,J)
%% 
range = model.m{midx}.wrange;
% store muscle points' parents in a vector
mp_p = zeros(1,size(range,2));
for i=1:size(range,2)
    mp_p(i) = model.m{midx}.P{1,range(i)}.parent;
end
% compute forward kinematics for the muscle points' parents
for i=1:max(mp_p)
    if model.parent(mp_p(i)) == 0
        r{i} = -skew(Xup{i}(1:3,1:3)'*Xup{i}(4:6,1:3));
        R{i} = Xup{1}(1:3,1:3)';
    else
        r{i} = r{model.parent(mp_p(i))} + R{model.parent(mp_p(i))}*(-skew(Xup{i}(1:3,1:3)'*Xup{i}(4:6,1:3)));
        R{i} = R{model.parent(mp_p(i))}*Xup{i}(1:3,1:3)';
    end
end
% compute world coordinates for muscle fixed points
% wp = [ps0,q0_s]; % wrapping points are concatenated in columns, fixed points first, contact points last
for i=1:size(mp_p,2)
    fp_w{i} = r{mp_p(i)} + R{mp_p(i)}*wp(:,i);
end

% compute velocities (in world coordinates) for wrapping points
for i=1:size(mp_p,2)
    v{i} = R{mp_p(i)}*(-skew(wp(:,i))*V{mp_p(i)}(1:3)+V{mp_p(i)}(4:6));
end

% compute straigth-line unit vector
e{1} = (fp_w{2}-fp_w{1})/norm(fp_w{2}-fp_w{1});

% Compute the moment arm Jacobian matrix eq (15)
Jm = e{1}'*([-R{mp_p(2)}*skew(wp(:,2)),R{mp_p(2)}]*J{mp_p(2)+1} - ... % insertion point
            [-R{mp_p(1)}*skew(wp(:,1)),R{mp_p(1)}]*J{mp_p(1)+1});     % origin point

end