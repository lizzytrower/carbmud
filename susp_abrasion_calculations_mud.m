%%Carbonate mud abrasion rate calculations
%Elizabeth Trower, University of Colorado Boulder, April 2018
%Michael Lamb, Caltech; Woodward Fischer, Caltech
%This code was designed with Matlab R2017b

%%%%%%%%%CALCULATIONS%%%%%%%%%%%%%%%%%%%%%

cdrag = (4/3).*(R.*g.*D)./(ws.^2); %[unitless]

%compute flow velocity
z0 = 3.*D./30; %This is a roughness coefficient that needs to be set.  Here set by grainsize
dz = (H-z0)./1000; %[m]
z=0;
z = [z0:dz:H]; %[m]
flow_vel = (ustar./0.41) .* log(z./z0); %[m/s]
flow_z = z; %[m]
Uf = sum((ustar./0.41) .* log(z./z0)).*dz./H; %[m/s]


%compute bed load height and velocity
hb = D.*1.44.*(tstage-1).^0.5; %height of the bed load layer [m]
Us = (R.*g.*D).^0.5.*1.56.*(tstage-1).^0.56; %bed load velocity [m/s]
Us1=Us;
%don't let particle velocity exceed fluid velocity - this is important
if gaurds2 == 1
Us(Us>Uf)=Uf;
end

%compute suspended load
if hb < H 
hb(hb<D)=D; %sets minimum ht of bed load layer [m]
b = hb; %bottom of the suspended load layer - same as height of bedload [m]
betta = 2; %Based on Scheingross et al. (2014) best fit
P = ws./(0.41.*ustar.*betta); %Rouse number [unitless]

%This Log scale cleans up the integration
di5 = 0;
i5=0;
res = 1000;
di5 = (log(H)-log(b))./res;
i5 = [log(b):di5:log(H)];
z=0;
dz=0;
z = exp(i5);
z(length(z))=H; %[m]
dz = diff(z); %[m]
dz = [dz(1) dz];
a1=sum((((1-(z(z>z0)./H))./(1-(b./H))).*(b./z(z>z0))).^P.*log(z(z>z0)./z0)...
    .*dz(z>z0))./ (Uf.*H) .* (ustar./0.41);
cb = c_b1;

%find concentration profile
c=0;
c(1) = cb;
c(2:length(z)+1) = cb.*(((1-(z./H))./(1-(b./H))).*(b./z)).^P ;
z=[0 z];
c(z==H)=0;

%calculate the fall distance
    gradc(1:length(c))=0;
    gradc(2:length(c)) = -diff(c);    
    Hfall = (1./cb).*sum(z.*gradc); %[m]
else
    hb = H;
    cb = 1./(Us.*hb);
    Hfall = hb;   
    a1=0;
end

if cb == 0
    Hfall = 0;
end

if P < 0.8
    Hfall = H./5;
end


%%%%%%Probability for fluctuations%%%%%%%%%%%%%%%%%%%%%%%
sig = ustar; %[m/s]
dx = sig./100; %the number of bins to subdivide [m/s]
X = [-6*sig:dx:6*sig]; %spread distribution for six sigma = w'/ws [m/s]
f = normpdf(X,0,sig); %centered at zero normal gausian distribution [m/s]
X = X./ws;  % to normalize as w'/ws same as psi [unitless]

%%%%%Calculate impact velocity due to gravity and turbulence
Scos = 1;  %cosine of the angle of the bed for impacts.  Assume flat for ocean

wfall = Scos.*((2.*(2./3).*D.*g./cdrag.*R)...
    .*(1-exp(-cdrag.*rho_f./rho_s.*(Hfall./Scos)./(2./3.*D)))).^0.5; %[m/s]
wfall(Hfall<=(0.5.*D))=0;

    psifall = wfall./ws; %[unitless]
    settlematrix = 0;
    settlematrix = psifall + X; %[unitless]
    settlematrix1=settlematrix;
    settlematrix(settlematrix<0) = 0; %no negative impacts
    psifall_turb = sum((settlematrix).*f).*dx; 
    psi_fall3 = sum((settlematrix.^3).*f).*dx;
    E1 = psi_fall3;  %erosion with turbulence


%Stokes number correction
wi_st = settlematrix;
wi_st((D.*wi_st.*ws.*rho_s./(18.*nu.*rho_f))<Stc) = 0; 
%critical stokes number of 70 based on Scheingross et al 2014
psi_fall3_st = sum((wi_st.^3).*f).*dx;
E1_st = psi_fall3_st;  %erosion with turbulence and stokes correction

if tstage <=1  %Set erosion to very small if particles stop moving
    E1_st = 0;
end

if Hfall == 0
    E1_st = 0;
end


