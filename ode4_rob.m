function [yout] = ode4_rob(model,F,controller,t0,ts,tf,y0)
% ODE4  Classical Runge-Kutta ODE solver.
NDOF = model.NB;
fe = zeros(6,1);

idx = 1;
y = y0;
yout = [y; zeros(12,1)];

for t = t0 : ts : tf-ts
    c_pd = controller(t, y(1:NDOF), y(NDOF+1:NDOF*2));
    
    [out1] = F(model, y, c_pd, fe);
    s1 = out1(1:NDOF*2);

    [out2] = F(model, y+ts*s1/2, c_pd, fe);
    s2 = out2(1:NDOF*2);

    [out3] = F(model, y+ts*s2/2, c_pd, fe);
    s3 = out3(1:NDOF*2);

    [out4] = F(model, y+ts*s3, c_pd, fe);
    s4 = out4(1:NDOF*2);

    y = y + ts*(s1 + 2*s2 + 2*s3 + s4)/6;
    
    yout(length(y)+1:end,idx) = [out1;c_pd];
    yout = [yout [y; zeros(12,1)]];
    
    idx = idx + 1;
    
end

[out1] = F(model, y, controller(t, y(1:NDOF), y(NDOF+1:NDOF*2)), fe);
yout(length(y)+1:end,idx) = [out1;c_pd];