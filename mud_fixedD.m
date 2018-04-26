%%Carbonate mud abrasion rate at fixed grain size
%Elizabeth Trower, University of Colorado Boulder, April 2018
%This code was designed with Matlab R2017b

clear

%set grain size
D1 = 420; %[um]
D1 = D1*10^-6; %grain size in [m]

ustar1 = 0.001:0.001:.2; %shear velocity in [m/s]

%physical parameters
rho_s = 2850; %particle density in [kg/m^3]
rho_f = 1025; %density of water in [kg/m^3]
R = (rho_s - rho_f)/rho_f; %submerged specific density [unitless]
kv = 25*10^6; %[unitless]
g = 9.81; %[m/s^2]
nu = 1.3*10^-6; %kinematic viscosity of water [m^2/s]
young = 144*10^9; %young's modulus [kg/m/s^2]
strength = 1*10^6; %tensile strength [kg/m/s^2]

tauc = 0.03; %Critical Shields number.  0.03 is good for sand.
gaurds2 = 1; %this sets limit to Ub if  = 1
Stc = 10;

%calculate settling velocity
CSF = 0.8;  %1 is for spheres, 0.8 is for natural
PS = 3.5;  %6 is for spheres, 3.5 is for natural
Dstar = (R.*g.*D1.^3)./(nu.^2);
X = log10(Dstar);
R1 = -3.76715+1.92944.*X - 0.09815.*(X.^2) - 0.00575.*(X.^3) + 0.00056.*...
    (X.^4);
R2 = log10(1-((1-CSF)./0.85))-(((1-CSF).^2.3).*tanh(X-4.6)) + 0.3.*...
    (0.5-CSF).*((1-CSF).^2).*(X-4.6);
R3 = (0.65-((CSF./2.83).*tanh(X-4.6))).^(1+((3.5-PS)./2.5));
Wstar = R3.*10.^(R2+R1);
ws1 = (R.*g.*nu.*Wstar).^(1./3);

%variables for impact rate eqn
Rep = (R*g.*D1).^(1/2).*D1./nu; %[unitless]
A_GP = 1.3*10^-7; %constant from Garcia and Parker
A1 = 0.36; %[unitless]
Vp = pi()/6.*D1.^3; %[m^3]

eps_v = kv*strength^2/(2*young); %kinetic energy per unit volume eroded [kg/m/s^2]

H = .5;  %Set Depth of water [m]

%pre-allocate space
c_b = zeros(length(D1),length(ustar1));
Ir = zeros(length(D1),length(ustar1));
Erate = zeros(length(D1),length(ustar1));

for nn = 1:length(ustar1)
    
    ustar = ustar1(nn); %[m/s]
    
    Z = ustar./ws1.*Rep.^0.6; %[unitless]
    c_b(:,nn) = A_GP.*Z.^5./(1+A_GP/0.3.*Z.^5); %[unitless]
    
    Ir(:,nn) = A1.*c_b(:,nn)./Vp'; %impact rate (without w_i) [1/m^3]

    c_b1 = c_b(:,nn);
    D = D1;
    ws = ws1;
    tau = ustar^2/(R*g*D1);
    tstage = tau/tauc;

    susp_abrasion_calculations_mud
    Ewi = E1_st*(g*D)^(3/2); %[m^3/s^3]

    
    V_i = 1/2.*Vp.*rho_s./eps_v; %volume eroded per impact (without w_i) [m*s^2]
    
    Erate(:,nn) = V_i.*Ir(:,nn).*Ewi; %[m^3/m^2/s]
    disp(nn)
    
end

Erate = real(Erate.*rho_s*1000*60*60*24*365); %convert to g/m^2/yr

figure
plot(ustar1,log10(Erate))

load('exptdata.mat','data_fixedD')
hold on
scatter(data_fixedD(:,2),data_fixedD(:,1))
