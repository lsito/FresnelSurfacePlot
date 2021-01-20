%% Fresnel surfaces

%% Cleaning 
clear;
clc;
close all; 

%% Parameters
% 
% eta1 = 1.618;
% eta2 = 1.620;
% eta3 = 1.627;
% 
% mu0 = 4*pi*10^(-7);
% eps0 = 8.854 * 10^(-12); 
% 
% c = 1/(sqrt(mu0*eps0));

% tau1 = 1/eta1;
% tau2 = 1/eta2;
% tau3 = 1/eta3;

tau1 = 0.8;
tau2 = 0.5;
tau3 = 0.2;



%% Main

%La prima superficie la plotto tutto, la seconda a metà
[teta, phi] = meshgrid(linspace(0,pi,50), linspace(0,2*pi,50));
[teta0, phi0] = meshgrid(linspace(0,pi,50), linspace(0,pi,50));

% rho = 1;

%Coordinate "intere"
a = sin(teta).*cos(phi);
b = sin(teta).*sin(phi);
c = cos(teta);

%Coordinate "a metà"
a0 = sin(teta0).*cos(phi0);
b0 = sin(teta0).*sin(phi0);
c0 = cos(teta0);


phi2 = PHI2(a ,b ,c , tau1, tau2, tau3);
psi4 = PSI4(a ,b ,c , tau1, tau2, tau3);

phi20 = PHI2(a0 ,b0 ,c0 , tau1, tau2, tau3);
psi40 = PSI4(a0 ,b0 ,c0 , tau1, tau2, tau3);

%Prima superficie intera
n1 = sqrt(1./(phi2 + sqrt((phi2).^2 - psi4)));
%Seconda superficie a metà
n2 = sqrt(1./(phi20 - sqrt((phi20).^2 - psi40)));

f1 = a .* n1;
f2 = b .* n1;
f3 = c .* n1;

g1 = a0 .* n2;
g2 = b0 .* n2;
g3 = c0 .* n2;

figure('Name','Surf');
plot1 = surf(f1,f2,f3)
hold
plot2 = surf(g1,g2,g3)
alpha(plot2,0.1)


figure('Name','Mesh');
plot3 = mesh(f1,f2,f3)
hold
plot4 = mesh(g1,g2,g3)
%alpha(plot2,0.3)

%% Definition of functions
function phi2 = PHI2(p ,q ,r , tau1, tau2, tau3)
    phi2 = p.^2./2.*(1./tau2 + 1./tau3) + q.^2./2.*(1./tau1 + 1./tau3) + r.^2./2.*(1./tau2 + 1./tau1); 
end

function psi4 = PSI4(p ,q ,r , tau1, tau2, tau3)
    psi4 = (p.^2 + q.^2 + r.^2).*(p.^2./(tau2.*tau3) + q.^2./(tau1.*tau3) + r.^2./(tau2.*tau1)); 
end