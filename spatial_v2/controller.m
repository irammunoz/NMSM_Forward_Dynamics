function u = controller(t, q, qd)

ref = [0.125*pi*sin(t) - 0.75*pi;
                          0.5*pi];
dref = [0.125*pi*cos(t);
                      0];
                      
error = ref - q;
derror = dref - qd;

P = 50;
D = P*0.015;

u = P*error + D*derror;