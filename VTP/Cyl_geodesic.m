function ydot = Cyl_geodesic(~,y,r)
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

ydot = zeros(4,1);
ydot(1) = y(3); %p
ydot(2) = y(4); %q
ydot(3) = -C111*ydot(1)^2 - 2*C121*ydot(1)*ydot(2) - C221*ydot(2)^2;
ydot(4) = -C112*ydot(1)^2 - 2*C122*ydot(1)*ydot(2) - C222*ydot(2)^2;
end