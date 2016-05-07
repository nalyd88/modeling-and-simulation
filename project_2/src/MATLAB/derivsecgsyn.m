function dxdt = derivsecgsyn(t,x,flag,rr,sfint,ti,ai,bi)

%% Work based on paper by McSharry "A Dynamical Model for Generating Synthetic 
% Electrocardiogram Signals"

% This file provides dxdt = F(t,x) taking input paramters: 
% rr: rr process 
% sfint: Internal sampling frequency [Hertz]
% Order of extrema: [P Q R S T]
% ti = angles of extrema [radians] 
% ai = z-position of extrema 
% bi = Gaussian width of peaks 


xi = cos(ti);
yi = sin(ti);
ta = atan2(x(2),x(1));
a0 = 1.0 - sqrt(x(1)^2 + x(2)^2);
ip = floor(t*sfint); 
w0 = 2*pi/rr(ip); %Equation (4) in paper 


fresp = 0.25;
zbase = 0.00015*sin(2*pi*fresp*t);

dx1dt = a0*x(1) - w0*x(2);
dx2dt = a0*x(2) + w0*x(1);

dti = rem(ta - ti, 2*pi);
dx3dt = - sum(ai.*dti.*exp(-0.5*(dti./bi).^2)) - 1.0*(x(3) - zbase); %Equation (1,2) in paper 

dxdt = [dx1dt; dx2dt; dx3dt];

