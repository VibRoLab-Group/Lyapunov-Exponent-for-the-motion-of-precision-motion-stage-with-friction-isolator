function dx = eqn(t,x)
global zeta sigma_0 sigma_1 sigma_2 mu_s mu_k a mu v_rv k_r alpha kappa
k_i=mu;
v_r=v_rv+x(5,1);
g=mu_s+(mu_s-mu_k)*exp(-a*abs(v_r));
dx(1,1)=x(2,1);
dx(2,1)=-2 * zeta * x(2,1) - x(1,1) - k_i * x(3,1) - k_r * x(1,1) + k_r * x(4,1) - 2 * kappa * x(2,1) + 2 * kappa * x(5,1);
dx(3,1)=x(1,1);
dx(4,1)=x(5,1);
dx(5,1)=-2 * kappa * alpha * (x(5,1) - x(2,1)) - k_r * alpha * (x(4,1) - x(1,1)) - alpha * (sigma_0* x(6,1) + sigma_1*(v_r  -abs(v_r)*sigma_0 * x(6,1) / g) + sigma_2 * v_r);
dx(6,1)=  (v_r  -abs(v_r) * sigma_0 * x(6,1) / g);
end
