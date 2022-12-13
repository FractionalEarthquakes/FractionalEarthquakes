clear all

close all



k=0.37;

alpha=0.79;

c=0.0078;

p=1.149; theta = p-1;

mc=4;

m1 = mc:0.05:8;

dm=m1-mc;



syms x

eqn = (k .* exp(alpha.*(x-mc)) - 1) ./ (k .* c^theta .* gamma(1-theta) .* exp(alpha.*(x-mc))) == 0; 

s = vpasolve(eqn,x);

s

m1_var = mc:0.05:double(s);

nu = (k .* exp(alpha.*(m1_var-mc)) - 1) ./ (k .* c^theta .* gamma(1-theta) .* exp(alpha.*(m1_var-mc)));



figure

plot(m1_var,nu,'k','Linewidth',2)

grid on

xlabel('m_1','FontSize',14)

ylabel('\nu','FontSize',18)

tstar = (-1./nu) .^ (1/theta);



figure

semilogy(nu,tstar,'k','Linewidth',2)

grid on

ylabel('t^*','FontSize',14)

xlabel('\nu','FontSize',18)



alphavar = 0.5:0.05:5;

kvar = 0.3:0.05:1;

[Xa, Yk]=meshgrid(alphavar,kvar);

eps = 0.5;

figure

subplot(2,1,1)

epstilde = (eps/c)^(theta);

f0 = Yk .* (gamma(1-theta) + epstilde);

f = (1./Xa) .* log(epstilde ./ f0);

surf(Xa,Yk,f)

zlim([0, ceil(max(max(f)))])

xlabel('\alpha','FontSize',12)

ylabel('K_0','FontSize',12)

zlabel('RHS','FontSize',12)

colorbar

text(2.7,0.9,1.5,['\theta = ', num2str(round(theta,2)), ',   t_0 = ', num2str(round(c,2))],'FontWeight','bold','FontSize',12)



pvar = 1:0.05:2;

thetavar = pvar-1;

cvar = 0.005:0.005:1;

[Xtheta, Yc]=meshgrid(thetavar,cvar);

subplot(2,1,2)

epstilde = (eps./Yc).^(Xtheta);

f0 = k .* (gamma(1-Xtheta) + epstilde);

f = (1./alpha) .* log(epstilde ./ f0);  

surf(Xtheta,Yc,f)

zlim([0, ceil(max(max(f)))])

xlabel('\theta','FontSize',12)

ylabel('t_0','FontSize',12)

zlabel('RHS','FontSize',12)

text(0.6,0.9,1.5,['\alpha = ', num2str(round(alpha,2)), ',   K_0 = ', num2str(round(k,2))],'FontWeight','bold','FontSize',12)

colorbar



