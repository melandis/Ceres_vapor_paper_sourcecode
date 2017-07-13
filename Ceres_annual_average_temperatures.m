%Ceres_annual_average_temperatures.m

%A code to calculate the Ceres annual average temperatures in linearly
%spaced intervals from -90 to 90N, assuming Ceres is a smooth billiard 
%ball. Can handle slopes but default is flat terrain.

%Written by M.E. Landis, used in M.E. Landis et al., 2017 JGR
%Created in MATLAB_2016a


clc
clear all


format long e
warning('off', 'MATLAB:circshift:ScalarShift')  %circshift can throw a warning in some versions of MatLab

f=0.5; %0/1 for fully explicit/implicit, 0.5 for Crank-Nicholson

%constants 
sigma=5.67e-8;           % Stefan-Boltzmann constant
f1au=1367;                              %solar flux at 1AU

%how many years for the model to be run

model_years=15; %total years to be run 
initial_years=5; %spin up time (temperatures will be average and reset to speed convergence after this number of years)

%degrees to radians conversions
d2r=pi/180;                             %convert degrees to radians
r2d=180/pi;                             %convert radians to degrees

%Set up the linearly spaced upper and lower limits of the latitude bands
Lat1=linspace(0,36,37)*5-90-2.5;   %lower limit of latitude bands
Lat1(1)=-90; 

Lat2=linspace(0,36,37)*5-90+2.5; %upper limit of latitude bands
Lat2(end)=90; 

aa=482640;   %a axis radius of the ellipsoid
cc=445570;   %c axis radius of the ellipsoid
rr = (aa*cc) ./ sqrt(aa^2*sin(0.5*(Lat1+Lat2).*(pi/180)).^2+cc^2.*cos(0.5.*(Lat1+Lat2).*(pi/180)).^2) ; 

area = 2.*pi*rr.^2.*(sin(Lat2*(pi/180))-sin(Lat1.*(pi/180))); 

avg_lat=(Lat1+Lat2)'./2; 

%Latitude and longitude
latd=avg_lat;         %north latitude, in degrees

lat_array=(latd./360)*2*pi; %lat and lon in radians, for the use in the rest of the code

%Slope dependence
slope=0     *d2r;                       %slope in degrees, will need to include 0 if running for the first time for a non-zero case 
az_slope=0  *d2r;                       %direction that the slope is pointing in degrees (0/360 is north)

%Model parameters
n_layers=200;     %number of layers in the model

TI_reg=15;        %thermal intertia of the regolith
rho_reg=1388.0;   %regolith density
c_reg=837.0;      %regolith specific heat

A_slope=0.09;            %Albedo for the surface being modeled 
eps_slope=0.95;          %Emissivity for the surface being modeled 

Q=0.002;                 %Upwelling flux at the bottom boundary condition  


%Set up arrays of thermal properties for the column of regolith
rho=rho_reg*ones(1,n_layers);        
c=c_reg*ones(1,n_layers); 
k=(TI_reg.^2/(rho(1)*c(1)))*ones(1,n_layers);


%Info on Ceres' orbit
a=2.76750591440571;                %Orbital semi-major axis
e=0.07582276595896797;             %Orbital eccentricity
lsp=302.11022    *d2r;             %Ls of perihelion in entered in degrees (converts automatically to radians)

au = 1.4959787061e11;
gc  = 6.67259e-11;               % Gravitational constant
sm  = 1.9891e30;                 % Solar mass
yr = 2*pi/sqrt(gc*sm/(a*au)^3);  % Year length 

obl=linspace(2,20,19) * d2r;         %Obliquity should be entered in degrees (converts in line to radians), this includes 2 to 20 degrees in one degree intervals

day_length=32667;


%Set up depths
dz_0=0.5*sqrt((k(2)*day_length)/(pi*rho(2)*c(2))); %initializing the first layer thickness based on the stability criterion
dz=dz_0*1.03.^(linspace(0,n_layers-1,n_layers));    %layer thicknesses are allowed to vary

%set reference temperatures
Tr=150;
Ts=150;

%Model times 
dt=600;

n_tsteps=round(yr/dt); 
time=linspace(0,yr,n_tsteps); 
ls=linspace(0,2*pi, n_tsteps); 

%Correct times so model years begin at the same HA/lst
HA_corr=mod((time/(day_length)), 1);
tsteps=yr;
tsteps=tsteps+((1-HA_corr(max(size(time)))*day_length)); 
time=linspace(0,yr,n_tsteps); 
n_tsteps=max(size(time)); 

%set reference temperatures
Tr=150;
Ts=150;

%Setting time steps
sf=0.1;  %Stability factor, the fraction of the maximum stable time step set (max. time step calculated from Courant criterion)
         %A jump will occur on highly inclined surfaces near sunrise/sunset
         %for ~1 time step because of the way SW is calculated. Be careful
         %on slopes greater than about 20 degrees because the time steps
         %must be smaller=>also running into issues with the sun rising
         %instantly

%Model times 

k3l=sqrt(gc*sm/((a*au)^3)); 
year = 2*pi/sqrt(gc*sm/(a*au)^3);                   %length of a Ceres year

%initialize temperature arrays (T_plot and T)
Ts_plot=zeros(1,n_tsteps);
T=200*ones(1,n_layers);

% Low-Res calculation of times at uniformly spaced true anomalies
n=10000; 

ta=(2*pi/n)*(linspace(0, n, n));    %create evenly spaced true anomalies
ea_inners= (e+cos(ta))./(1+e*cos(ta));  %calculate eccentric anomalies
ea=acos(ea_inners);
ma=ea-e*sin(ea);                        %calculate mean anomalies
t=ma/k3l;                               %time along Mars orbital path (will be irregularly spaced)

t(ta>pi)=yr-t(ta>pi);                      %correction of 2nd half of the orbit
 

%High-Res interpolation of true anomalies at uniformly spaced times

t2  = linspace(0,yr, n_tsteps); %Using uniformly spaced times to calculate eccentric anamolies
ta2= interp1(t,ta, t2); 
ea2_inners= (e+cos(ta2))./(1+e*cos(ta2)); %calculate the eccentric anomalies at these times 
ea2=acos(ea2_inners);  

% Solar distance/declination/longitude and local time calculations

sol_dist=a*(1-e*cos(ea2));              %solar distance
ls =mod(ta2+lsp, 2*pi);            %solar longitude


%calculate sine and cosine of solar declination for all obliquities
for n=1:max(size(obl))
sin_dec(n,:)=sin(obl(n)).*sin(ls);
cos_dec(n,:)=sqrt(1-sin_dec(n,:).^2); 
end

%initialize temperature arrays (T_plot and T)
Ts=zeros(1,n_tsteps);
T=Tr*ones(1,n_layers); 

for obl_v=1:max(size(obl))
for w=max(size(slope))
    slope=slope(w)
    az_slope=az_slope(w)
for q=1:max(size(lat_array))
    lat=lat_array(q)

HA=mod((time/(day_length)*2*pi), 2*pi)-pi;  %starts at dawn 

cosi=(sin(lat).*sin_dec(obl_v,:)+cos(lat).*cos_dec(obl_v,:).*cos(HA)); %sine of the solar incidence angle
sini=sqrt(1.0-cosi.^2);                                                 %cosine of the solar incidence angle

cosi(cosi<0)=0;                                                         %prevents cos(incidence angle) from being negative

az_star_start=((sin_dec(obl_v,:))-sin(lat)*cosi)./(cos(lat)*sini); %solar azimuth
az_star_start(az_star_start>1)=1;
az_star_start(az_star_start<-1)=-1;
az_star=acos(az_star_start);
az_star(HA'>0)=(2*pi)-az_star(HA'>0); %only flip this when hour angle is in the afternoon 

cosi_s=cosi.*cos(slope)+sini.*sin(slope).*cos(az_slope-az_star);  %force this to be between zero and one
cosi_s(cosi_s<0)=0;
cosi_s(cosi==0)=0; 

flux_flat=(f1au./(sol_dist.^2)); 
SW=       (f1au./(sol_dist.^2)).*cosi; 
SW_slope= (f1au./(sol_dist.^2)).*cosi_s;


%Calculate the sky view factor, E=1 for flat terrain, E=cos^2(slope/2) for
%a sloped surface on a flat plane

E=cos(slope/2)*cos(slope/2); 

%Constants for both explicit and implicit versions of the model
a2=(3.*dz(1).*0.5.*eps_slope.*sigma)./k(1); 
b_constants=(2.*eps_slope.*sigma.*dz(1))./k(1); %same for implicit and explicit

alpha_constants=dt./(c.*rho.*dz); 
alpha_u=alpha_constants.*(2.*k.*circshift(k,[0, 1]))./(k.*circshift(dz,[0, 1])+circshift(k,[0, 1]).*dz);
alpha_d=alpha_constants.*(2.*k.*circshift(k,[0,-1]))./(k.*circshift(dz,[0,-1])+circshift(k,[0,-1]).*dz);

alpha_u(1)=0; 
alpha_d(n_layers)=0;

beta=(k(1)*dt)/(c(1)*rho(1)*dz(1)*dz(1));

%EXPLICIT MODEL: Thermal parameters 

ke1=((1-f).*alpha_u).*ones(1,n_layers); 
ke1(1)=0; 

ke2=(1-(1-f).*alpha_u-(1-f).*alpha_d).*ones(1,n_layers); 
ke2(n_layers)=(1-(1-f).*alpha_u(n_layers)); 

ke3=alpha_d.*(1-f).*ones(1,n_layers); 
ke3(n_layers)=0;

boundary=zeros(1,n_layers); 
boundary(n_layers)=(dt*Q)/(rho(n_layers)*c(n_layers)*dz(n_layers));

%IMPLICIT MODEL: thermal parameters 

ki1=-f.*alpha_u.*ones(1,n_layers); 
ki1(1)=[ ]; 

ki2=(1+f.*alpha_u+f.*alpha_d).*ones(1,n_layers); 
ki2(n_layers)=(1+f.*alpha_u(n_layers)); 

ki3=-f.*alpha_d.*ones(1,n_layers); 
ki3(n_layers)=[ ];

S=sparse(diag(ki3,1)+diag(ki2)+diag(ki1, -1));  


if slope > 0 
    load('tempvars.mat', 'Ts', 'A_slope', 'eps_slope')
    T_flat=Ts; 
else T_flat=zeros(1,max(size(time)));
end


%Calculate incoming solar flux, minus effects from terrain since that must
%be calculated based on flat floor temperatures
F_solar = SW_slope; 

F_VIS_terr = (A_slope.*SW_slope).*(1-E); 

F_IR_terr = eps_slope*sigma*T_flat.^4.*(1-E); 

f_vis_incoming =F_solar+ F_VIS_terr; 
f_ir_incoming= F_IR_terr; 

fi_vis_incoming=circshift(f_vis_incoming,[0,-1]);
fi_ir_incoming=circshift(f_ir_incoming,[0,-1]);


%Calculate temperatures for a sloping surface
for m=1:model_years
for n=1:n_tsteps
    
    b=1./(1+b_constants*Tr^3); 
    ke2(1)=1-(1-f)*(alpha_d(1)+(2-2*b)*beta);
    
    ae=(dz(1)/(2.*k(1))).*((1-A_slope)* f_vis_incoming(n)+eps_slope* f_ir_incoming(n)+3.*eps_slope.*sigma*Tr^4).*b;
    ai=(dz(1)/(2.*k(1))).*((1-A_slope)*fi_vis_incoming(n)+eps_slope* fi_ir_incoming(n)+3.*eps_slope.*sigma*Tr^4).*b;
    
    boundary(1)=2*beta*((1-f)*ae+f*ai); 
    
    T_present=(ke3.*circshift(T,[0,-1])+ke2.*T+ke1.*(circshift(T,[0,1])))+boundary;
    
    S(1,1)=1+f*(alpha_d(1)+(2-2*b)*beta);  
   
    T=S\(T_present');
    T=T';
 
    Tr=ai+b*T(1);
    Ts(n)=Tr;


    
end 

m

if m==initial_years
     T_average=mean(Ts);
     Tr=T_average;
     T=T_average*ones(1,n_layers); 
     Ts=zeros(1,max(size(n_tsteps)));
     
     
end



end

Ts_slope_min=min(Ts)
Ts_slope_max=max(Ts)
Ts_slope_avg=mean(Ts)

T_lat_annual_avg(q)=mean(Ts);

end


M=[Lat1' Lat2' area' T_lat_annual_avg']
filename = ['Ceres_latitudinal_annual_avg_temps' num2str(obl(obl_v)*r2d) 'obliquity' num2str(slope*r2d) 'slope' num2str(az_slope*r2d) 'azimuth'  '.csv' ];

csvwrite(filename, M)
 
%IF slope =0, then save temperatures, emissivities, and albedos of the flat
%terrain case 

if slope == 0
    save('tempvars.mat', 'Ts', 'A_slope', 'eps_slope')
end

end
end



