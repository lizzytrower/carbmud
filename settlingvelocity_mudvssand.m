%%Carbonate mud settling velocity and advection length scale calculations
%Elizabeth Trower, University of Colorado Boulder, April 2018
%This code was designed with Matlab R2017b

clear

%calculate mud settling time
Dmud = 100*10^-9; %particle intermediate diameter in m
Amud = pi()*(Dmud/2)^2; %cross sectional area in m^2
Lengthmud = 1*10^-6; %particle length in m
Vmud = Amud*Lengthmud;

rho_s = 2850; %density of aragonite [kg/m^3]
rho_f = 1025; %density of seawater in [kg/m^3]
R = (rho_s - rho_f)/rho_f; %submerged specific density [unitless]

nu = 1.3*10^-6; %kinematic viscosity of water [m^2/s]

g = 9.81; %gravitational acceleration [m/s^2]

wsmud = 1/12*R*g/nu*Vmud*Dmud/Amud;
settledistance = 1:.1:5; %[m]
settletime_m = settledistance./wsmud./60./60;

%calculate mud floc settling time
Dfloc = 50*10^-6;
rho_floc = rho_f + 8375*(Dfloc*10^6)^-1.084; %from Fennessy et al 1994
R_floc = (rho_floc - rho_f)/rho_f;

CSF = 0.8;  %Corey Shape Factor; 1 is for spheres, 0.8 is for natural
PS = 3.5;  %Powers Roundness; 6 is for spheres, 3.5 is for natural

Dstar_floc = (R_floc.*g.*Dfloc.^3)./(nu.^2);
Xf = log10(Dstar_floc);
R1f = -3.76715+1.92944.*Xf - 0.09815.*(Xf.^2) - 0.00575.*(Xf.^3) + ...
    0.00056.*(Xf.^4);
R2f = log10(1-((1-CSF)./0.85))-(((1-CSF).^2.3).*tanh(Xf-4.6)) + ...
    0.3.*(0.5-CSF).*((1-CSF).^2).*(Xf-4.6);
R3f = (0.65-((CSF./2.83).*tanh(Xf-4.6))).^(1+((3.5-PS)./2.5));
Wstarf = R3f.*10.^(R2f+R1f);

wsfloc = (R_floc.*g.*nu.*Wstarf).^(1./3);
settletime_f = settledistance./wsfloc./60./60;

%calculate sand settling time
Ds = 500; %range of grain sizes [um]
Ds = Ds*10^-6; %convert grain size to [m]
CSF = 0.8;  %Corey Shape Factor; 1 is for spheres, 0.8 is for natural
PS = 3.5;  %Powers Roundness; 6 is for spheres, 3.5 is for natural

Dstar = (R.*g.*Ds.^3)./(nu.^2);
X = log10(Dstar);
R1 = -3.76715+1.92944.*X - 0.09815.*(X.^2) - 0.00575.*(X.^3) + ...
    0.00056.*(X.^4);
R2 = log10(1-((1-CSF)./0.85))-(((1-CSF).^2.3).*tanh(X-4.6)) + ...
    0.3.*(0.5-CSF).*((1-CSF).^2).*(X-4.6);
R3 = (0.65-((CSF./2.83).*tanh(X-4.6))).^(1+((3.5-PS)./2.5));
Wstar = R3.*10.^(R2+R1);
wssand = (R.*g.*nu.*Wstar).^(1./3);

settletime_s = settledistance./wssand./60./60;

%set a flow velocity
U = .8; %[m/s]

la_mud = U.*settledistance./wsmud;
la_floc = U.*settledistance./wsfloc;
la_sand = U.*settledistance./wssand;


advectionratio_mudsand = wssand./wsmud;
advectionratio_mudfloc = wsfloc./wsmud;
advectionratio_flocsand = wssand./wsfloc;