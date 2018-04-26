%%Carbonate mud abrasion rate
%Elizabeth Trower, University of Colorado Boulder, April 2018
%This code was designed with Matlab R2017b

clear

D1 = 100:5:1500;
D1 = D1*10^-6; %grain size in [m]

ustar1 = 0.001:0.001:0.15; %shear velocity in [m/s]

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

H = 0.5;  %Set Depth of water [m]

%pre-allocate space
c_b = zeros(length(D1),length(ustar1));
Ir = zeros(length(D1),length(ustar1));
Erate = zeros(length(D1),length(ustar1));
Ewi = zeros(length(D1),1);

for nn = 1:length(ustar1)
    
    ustar = ustar1(nn); %[m/s]
    
    Z = ustar1(nn)./ws1.*Rep.^0.6; %[unitless]
    c_b(:,nn) = A_GP.*Z.^5./(1+A_GP/0.3.*Z.^5); %[unitless]
    
    Ir(:,nn) = A1.*c_b(:,nn)./Vp'; %impact rate (without w_i) [1/m^3]
    
    for mm = 1:length(D1)
        
        c_b1 = c_b(mm,nn);
        D = D1(mm);
        ws = ws1(mm);
        tau = ustar^2/(R*g*D);
        tstage = tau/tauc;
        
        susp_abrasion_calculations_turb10_lizzy_mud
        Ewi(mm) = E1_st*(g*D)^(3/2); %[m^3/s^3]
        
    end
    
    V_i = 1/2.*Vp.*rho_s./eps_v; %volume eroded per impact (without w_i) [m*s^2]
    
    Erate(:,nn) = V_i'.*Ir(:,nn).*Ewi; %[m^3/m^2/s]
    disp(nn)
    
end

Erate = real(Erate.*rho_s*1000*60*60*24*365); %convert to g/m^2/yr

P = [0.8 1.2 2.5 7.5];
ustar_draw = zeros(length(ws1),length(P));
    
for pp = 1:length(P)
    
    betta = 2;

    ustar_draw(:,pp) = ws1./(P(pp).*.41.*betta);
    
end


figure 
p1 = contourf(ustar1,D1*10^6,log10(Erate),(0:.5:6),'LineStyle','none');
ylabel('grain diameter, D (\mum)')
xlabel('shear velocity, u_* (m/s)')
xlim([0 0.15])
ylim([0 1500])

hold on

p2 = plot(ustar_draw(:,1),D1*10^6,'-.k');
p3 = plot(ustar_draw(:,2),D1*10^6,'--k');
p4 = plot(ustar_draw(:,3),D1*10^6,'-k');

load('exptdata.mat','data_all')
p5 = scatter(data_all(1,3),data_all(1,1),100,data_all(1,2),'filled',...
    'MarkerEdgeColor','k');
p6 = scatter(data_all(9:12,3),data_all(9:12,1),100,data_all(9:12,2),'filled',...
    'MarkerEdgeColor','k');
p7 = scatter(data_all(2:8,3),data_all(2:8,1),100,data_all(2:8,2),'d','filled',...
    'MarkerEdgeColor','k');

colorbar
