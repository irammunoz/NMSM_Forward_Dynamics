function T = get_HT(Xup, n, model)
% Comupte The Homogeneous Transformation of the n-th joint expressed in 
% the inertial frame

X_temp = eye(6);

par(1) = n;
for i=1:model.NB
    if model.parent(par(i))==0
        break
    else
        par(i+1) = model.parent(par(i));
    end
end

for i=size(par,2):-1:1
    X_temp = Xup{par(i)} * X_temp;
end

T = pluho(X_temp);
T = inv(T);

end