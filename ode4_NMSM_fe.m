function [yout, wrap_points] = ode4_NMSM_fe(model,F,e,wrapc,t0,ts,tf,y0,wrap_flag,fe)
% ODE4  Classical Runge-Kutta ODE solver.
NMUS = length(model.m);
NDOF = model.NB;
tau = zeros(model.NB,1);

idx = 1;
y = y0;
yout = [y; zeros(78,1)];

wrap_points = cell(1 + tf/ts,1);

for t = t0 : ts : tf-ts
    ex_s = e(t);
    
    [out1, m1w_p, m2w_p] = F(model, y, wrapc, ex_s, tau, fe(:,idx), wrap_flag);
    s1 = out1(1:NDOF*2+NMUS*2+1);

    [out2, ~, ~] = F(model, y+ts*s1/2, wrapc, ex_s, tau, fe(:,idx), wrap_flag);
    s2 = out2(1:NDOF*2+NMUS*2+1);

    [out3, ~, ~] = F(model, y+ts*s2/2, wrapc, ex_s, tau, fe(:,idx), wrap_flag);
    s3 = out3(1:NDOF*2+NMUS*2+1);

    [out4, ~, ~] = F(model, y+ts*s3, wrapc, ex_s, tau, fe(:,idx), wrap_flag);
    s4 = out4(1:NDOF*2+NMUS*2+1);

    y = y + ts*(s1 + 2*s2 + 2*s3 + s4)/6;
    
    yout(length(y)+1:end,idx) = [out1;ex_s];
    yout = [yout [y; zeros(78,1)]];
    
    wrapc = out1(NDOF*2+NMUS*2+2:NDOF*2+NMUS*2+2+NMUS*15);
    
    wrap_points{idx} = {m1w_p, m2w_p};
    
    idx = idx + 1;
    
end

[out1, m1w_p, m2w_p] = F(model, y, wrapc, e(t), tau, fe(:,idx), wrap_flag);
yout(length(y)+1:end,idx) = [out1;ex_s];
wrap_points{idx} = {m1w_p, m2w_p};