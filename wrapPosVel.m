function [fp_w, e, v, Jm] = wrapPosVel(model,midx,widx,Xup,wp,wg,V,J)
%% 
range = model.m{midx}.wrange;
geom_p = model.appearance.wrapCyl{widx}.parent;
% store muscle points' parents in a vector
% mp_p = zeros(1,size(model.m{midx}.P,2));
mp_p = zeros(1,size(range,2));
for i=1:size(range,2)%size(model.m{midx}.P,2)
    mp_p(i) = model.m{midx}.P{1,range(i)}.parent;
end
% compute forward kinematics for the muscle points' parents
for i=1:max(mp_p)
    if model.parent(mp_p(i)) == 0
        r{i} = -skew(Xup{i}(1:3,1:3)'*Xup{i}(4:6,1:3)); % position expressed in world frame
        R{i} = Xup{1}(1:3,1:3)';                        % orientation expressed in world frame
    else
        r{i} = r{model.parent(mp_p(i))} + R{model.parent(mp_p(i))}*(-skew(Xup{i}(1:3,1:3)'*Xup{i}(4:6,1:3)));
        R{i} = R{model.parent(mp_p(i))}*Xup{i}(1:3,1:3)';
    end
end
% compute world coordinates for muscle fixed points
% wp = [p0_s,q0_s,p_s,q_s]; % wrapping points are concatenated in columns, fixed points first, contact points last
for i=1:size(mp_p,2)
    fp_w{i} = r{mp_p(i)} + R{mp_p(i)}*wp(:,i);
end
% compute world coordinates for muscle moving points
parent = model.m{midx}.P{1,range(1)}.parent;
for i=size(mp_p,2)+1:size(mp_p,2)+2
    fp_w{i} = r{mp_p(parent)} + R{mp_p(parent)}*wp(:,i);
end

% compute velocities (in world coordinates) for wrapping points
for i=1:size(mp_p,2)
    v{i} = R{mp_p(i)}*(-skew(wp(:,i))*V{mp_p(i)}(1:3)+V{mp_p(i)}(4:6));
end
% compute velocities (in world coordinates) for wrapping geometry
v{size(mp_p,2)+1} = R{geom_p}*(-skew(wg)*V{geom_p}(1:3)+V{geom_p}(4:6));

% compute straigth-line unit vectors
e{1} = (fp_w{3}-fp_w{1})/norm(fp_w{3}-fp_w{1});
e{2} = (fp_w{2}-fp_w{4})/norm(fp_w{2}-fp_w{4});

% Compute the moment arm Jacobian matrix eq (14)
Jm = e{1}'*[-R{geom_p}*skew(wg),R{geom_p}]*J{geom_p+1}...           % origin point
    -e{1}'*[-R{mp_p(1)}*skew(wp(:,1)),R{mp_p(1)}]*J{mp_p(1)+1}...   % geodesic p
    +e{2}'*[-R{mp_p(2)}*skew(wp(:,2)),R{mp_p(2)}]*J{mp_p(2)+1}...   % geodesic q
    -e{2}'*[-R{geom_p}*skew(wg),R{geom_p}]*J{geom_p+1};             % insertion point

end