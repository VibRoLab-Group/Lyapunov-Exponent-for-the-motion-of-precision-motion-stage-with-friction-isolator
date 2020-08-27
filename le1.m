%This code determines the L.E based on the nonlinear floquet multiplier for
%the periodic solution of the system
clear all;
clc;
tic;
global zeta sigma_0 sigma_1 sigma_2 mu_s mu_k a mu v_rv k_r alpha kappa
sigma_0=110;
sigma_1=1.37;
sigma_2=0.0823;
mu_s=0.44;
mu_k=0.35;
a=2.5;
zeta=0.02;
kappa=0.001;
v_rv=1;
k_r=0.5;
alpha=2;
epsilon=0.0001;%This signifies perturbation to initial conditions corresponding to steady periodic solution
noi=150;%This defines the number of points for which L.E has to be determined
LE=zeros(noi,1);
floquet_multipliers1=zeros(6,noi);
F1=zeros(6,6);
%reult.mat file contains the steady periodic solutions for different values
%of integral gain k_i. In this file, x is 7timesN array, in which first row
%signifies bifurcation parameter (in this case it is k_i), and other  rows
%defines state space of the system, with the exception of 3rd row. 3rd rwo
%deines the time-period of solution. 
load result
x1=x;
for jjj=1:noi
    jj=jjj;
         mu=x1(1,jj);
         tic
         x0=x1(2:end,jj);%initial conditions for steady periodic solution, whose stability has to be determined
         x0(2)=0;
         tspan=linspace(0,x1(3,jj),25000);
         op = odeset('RelTol',1e-13,'AbsTol',1e-13);
         [t2,x2]=ode45('eqn',tspan,x0,op);
       for k=1:6
            ics=x0;
            delta=epsilon*((ics(k)));
            ics(k)=ics(k)+delta;%In this step we determine the solution with perturbation in initial conditions
            [t3,x3]=ode45('eqn',tspan,ics,op);
            F1(:,k)=(x3(end,:)'-x2(end,:)')/(delta);
       end
        toc
        k=2;
            ics=x0;
            delta=epsilon;
            ics(k)=ics(k)+delta;
            [t3,x3]=ode45('eqn',tspan,ics,op);
            F1(:,k)=(x3(end,:)'-x2(end,:)')/(delta);
       
floquet_multipliers1(:,1)=eig(F1);
kkk=0;
for kk=1:6
tt=fun(floquet_multipliers1(kk,1));
if tt==0
pp(1+kkk)=floquet_multipliers1(kk,1);
kkk=kkk+1;
end
end
%For this we use the concept that real part of Floquet exponent represent
%L.E.
LE(jjj,1)=max(real(log((pp))/x1(3,jj)));
      save lyapyanov floquet_multipliers1 LE
end
