function u = excitation_signal(t)
u = zeros(2,1);

A = [0.58, 0.1];
tau = [1.5, 1.5];
tau_r = [0.5, 0.5];
T = [3, 3];
f0 = [1/T(1), 1/T(2)];

t_offset = [0, 1.5];
u_offset = [0.32, 0.05];

sumf = [0,0];
for i=1:2
    for n=1:1000
        d1 = (sin(n*pi*f0(i)*tau(i))/(n*pi*f0(i)*tau(i)));
        d2 = (sin(n*pi*f0(i)*tau_r(i))/(n*pi*f0(i)*tau_r(i)));

        d3 = pi*n*f0(i);
        d4 = (t+t_offset(i));
        d5 = (tau(i)-tau_r(i));

        d6 = cos(2*d3*d4-d3*d5);
        sumf(i) = sumf(i) + d1*d2*d6;
    end
    u(i,1) = u_offset(i) + A(i)*(tau(i)/T(i)) + 2*A(i)*(tau(i)/T(i))*sumf(i);
end
