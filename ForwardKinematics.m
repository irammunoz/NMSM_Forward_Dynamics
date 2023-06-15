function [x, y, z, phi, theta, psi, RL] = ForwardKinematics(model, q)

kp = 0;
if not(isempty(model.KXtree))
    kp = length(model.KXtree);
end

lq = size(q,1);

Xa = cell(1,model.NB);
x = zeros(model.NB+kp, lq);
y = zeros(model.NB+kp, lq);
z = zeros(model.NB+kp, lq);
phi = zeros(model.NB+kp, lq);
theta = zeros(model.NB+kp, lq);
psi = zeros(model.NB+kp, lq);

% ee_ind = 4;

for j=1:lq
    for i=1:model.NB
        [ XJ, S_t ] = Rodrigues(model.ar{i},q(j, i));
        RL{i,j} = XJ(1:3,1:3);
        Xa{i} = XJ * model.Xtree{i};
        if model.parent(i) ~= 0
            Xa{i} = Xa{i} * Xa{model.parent(i)};
        end
        if size(Xa{i},1) == 3		% Xa{i} is a planar coordinate xform
            [theta,r] = plnr(Xa{i});
            X = rotz(theta) * xlt([r;0]);
            T = pluho(X);
        else
            T = pluho(Xa{i});
        end
        Tdisp = inv(T);		% displacement is inverse of coord xform

        x(i, j) = Tdisp(1, 4);
        y(i, j) = Tdisp(2, 4);
        z(i, j) = Tdisp(3, 4);
        [phi(i, j), theta(i, j), psi(i, j)] = R2rpy(Tdisp(1:3,1:3));
    end
    
    % end points
    for k=1:kp
        T_e = pluho(model.KXtree{k} * Xa{model.kparent(k)});
        T_e_inv = inv(T_e);
        x(i+k, j) = T_e_inv(1, 4);
        y(i+k, j) = T_e_inv(2, 4);
        z(i+k, j) = T_e_inv(3, 4);
        [phi(i+k, j), theta(i+k, j), psi(i+k, j)] = R2rpy(T_e_inv(1:3,1:3));
    end
        
end
end
   
    
function [phi, theta, psi] = R2rpy(R)
sy = sqrt(R(1, 1) * R(1, 1) + R(2, 1) * R(2, 1));
singular = sy < 1e-69;

if not(singular)
    x = atan2(R(3, 2), R(3, 3));
    y = atan2(-R(3, 1), sy);
    z = atan2(R(2, 1), R(1, 1));
else
    x = atan2(-R(2, 3), R(2, 2));
    y = atan2(-R(3, 1), sy);
    z = 0;
end

phi = x * 180 / pi;
theta = y * 180 / pi;
psi = z * 180 / pi;
end