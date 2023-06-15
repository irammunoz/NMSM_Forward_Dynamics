function t = Cyl_tangent(u,v,up,vp,r)
rv = [0;0;1];
% t = ru*up + rv*vp;
t = zeros(length(u),3);
for i=1:length(u)
    ru = [-r*sin(u(i));r*cos(u(i));0];
    t(i,:) = ru*up(i) + rv*vp(i);
end

end