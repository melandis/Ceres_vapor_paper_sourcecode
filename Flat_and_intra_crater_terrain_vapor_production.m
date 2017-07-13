%Flat_and_intra_crater_terrain_vapor_production.m

%Code to calculate the vapor production from exposed surface ice on flat terrain or within an impact crater on
%Ceres.

%Written by M.E. Landis, used in M.E. Landis et al., 2017 JGR
%Created in MATLAB_2016a


clc
clear all


format long e
warning('off', 'MATLAB:circshift:ScalarShift')

f=0.5; %0/1 for fully explicit/implicit, 0.5 for Crank-Nicholson

%how many years for the model to be run
model_years=10; 
initial_years=5; 

%degrees to radians conversions
d2r=pi/180;                             %convert degrees to radians
r2d=180/pi;                             %convert radians to degrees

%Can run this code in two modes: flat terrain and Oxo crater. Setting d to
%zero will force E~1, which is the skyview factor for flat terrain. Comment
%out the block of code for the case that is NOT being run. Always run a
%flat terrain case in the directory being used first in order to set up the
%tempvars.m file that provides the ambient surrounding terrain temperatures
%for use in the re-radiation from terrain

%Flat terrain 
D=10e3; 
d=0; 
r=0; 

lat_array=[0] * d2r ; 

%Oxo crater parameters and crater geometry calculations
% D=10e3;                         %diameter in m
% d=1.5e3;                        %depth in m, assuming a parabolic shape for the crater
% lat_array=42.2 * d2r;           %latitude of the center of Oxo crater
% 
% r=1/3*D;            %[0, 100, 1e3, 1/3*D, 4e3];         %radii of interest within the crater, in m, !!!!make sure it is less than half the diameter please!!


%Slope dependence

theta=atan(r*(4*d)/(D^2)); 
gamma=tan(theta)/(2*(d/D)); 

gamma=(2*r)/D; 


slopes=atan(2.*tan(theta)); 
az_slopes=[0]    *d2r; %[0, 90, 180, 270] *d2r;  %0 is north, etc. 

phi_skyview=linspace(0,2*pi, 1000);                   %direction(s) of interest

for x=1:max(size(r))
      
height=atan(((2.*d/D).*(1-gamma(x).^2))./(sqrt(1-gamma(x).^2.*sin(phi_skyview).^2)+gamma(x).*cos(phi_skyview))); 
E(x)=(0.5/pi) * (sum(cos(height).^2)) * (phi_skyview(2)-phi_skyview(1)); 
 

end


%Model parameters
n_layers=200;     %number of layers in the model

TI_ice=2100;        %thermal intertia of the ice

A_flat=0.09;     %Albedo of surrounding terrain
A_slope=0.135;   %Albedo of the water ice patch        
eps_slope=0.95; 
sigma=5.67e-8;  
Q=0.002; 
 
T0=150;


%Set up arrays of thermal properties for the column, including any
%variation
rho=925.0*ones(1,n_layers);        
rho(1)=925.0; 
c=1615.0*ones(1,n_layers); 
c(1)=1615.0;
k=(TI_ice.^2/(rho(2)*c(2)))*ones(1,n_layers);
k(1)=(TI_ice.^2/(rho(1)*c(1))); 
inv_k=1/k(1); 

f1au=1367;                              %solar flux at 1AU


%Info on planet's orbit (doing this one for the MOON) (this would change for a previous epoch)

a=2.76750591440571;                %Orbital semi-major axis
e=0.07582276595896797;             %Orbital eccentricity
lsp=302.11022    *d2r;                    %Ls of perihelion in DEGREES
obl=4.03        * d2r;                            %Obliquity in DEGREES
day_length=32667;
                   
au = 1.4959787061e11;
gc  = 6.67259e-11;               % Gravitational constant
sm  = 1.9891e30;                 % Solar mass
yr = 2*pi/sqrt(gc*sm/(a*au)^3);


%Set up depths
                            %dz=0.05*ones(1,n_layers); %constant dz

dz_0=0.5*sqrt((k(2)*day_length)/(pi*rho(2)*c(2))); %70.79e-6; %varying dz
dz=dz_0*1.03.^(linspace(0,n_layers-1,n_layers)); 


%set reference temperatures
Tr=150;
Ts=150;

%Model times 

sf=0.1;  %Stability factor, the fraction of the maximum stable time step set (max. time step calculated from Courant criterion)
         %A jump will occur on highly inclined surfaces near sunrise/sunset
         %for ~1 time step because of the way SW is calculated. Be careful
         %on slopes greater than about 20 degrees because the time steps
         %must be smaller=>also running into issues with the sun rising
         %instantly
dt=600; %min(sf.*(rho.*c.*dz.*dz.*0.5.*inv_k));  

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

%more times!

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

sin_dec=sin(obl)*sin(ls);              %sine of the solar declination
cos_dec=sqrt(1-sin_dec.^2); 
 
%initialize temperature arrays (T_plot and T)
Ts=zeros(1,n_tsteps);
T=Tr*ones(1,n_layers); 

HA=mod((time/(day_length)*2*pi), 2*pi)-pi;  %starts at dawn 



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

for w=1:max(size(slopes))
    slope=slopes(w)
    az_slope=az_slopes(w)
for q=1:max(size(lat_array))
    lat=lat_array(q)
    for p=1:max(size(r))
    slope=slopes(p); 

%Calculate incidence angles 
cosi=(sin(lat).*sin_dec+cos(lat).*cos_dec.*cos(HA)); 
sini=sqrt(1.0-cosi.^2);


az_star_start=((sin_dec)-sin(lat)*cosi)./(cos(lat)*sini); %solar azimuth
az_star_start(az_star_start>1)=1;
az_star_start(az_star_start<-1)=-1;
az_star=acos(az_star_start);
az_star(HA'>0)=(2*pi)-az_star(HA'>0); %only flip this when hour angle is in the afternoon 

cosi_s=cosi.*cos(slope)+sini.*sin(slope).*cos(az_slope-az_star);  %force this to be between zero and one
cosi_s(cosi_s<0)=0;
cosi_s(cosi==0)=0; 
cosi(cosi<0)=0; 

cosi(cosi<0)=0; 

%Now, calculate a factor to figure out where the sun is below the taller
%horizon due to the crater wall

phi_comp=pi-az_slope+az_star; %phi for comparison with the 90-incidence angle (90-i)
height_comp=atan(((2.*d/D).*(1-gamma(p).^2))./(sqrt(1-gamma(p).^2.*sin(phi_comp).^2)+gamma(p).*cos(phi_comp)));
incidence_angle=acos(cosi); 

solar_flux=(f1au./(sol_dist.^2)); 
solar_flux(((pi/2)-incidence_angle)<height_comp)=0;  %pi/2-incidence angle needs to be in parenths 


SW=solar_flux.*cosi; 

SW_slope= (solar_flux).*cosi_s;
flux_flat=(f1au./(sol_dist.^2));

%Calculate maximum noontime flux
cosi_noon=(sin(lat).*sin_dec+cos(lat).*cos_dec);
SW_noon=(f1au./(sol_dist.^2)).*cosi_noon; 

if slope > 0 
    load('tempvars.mat', 'Ts', 'A_slope', 'eps_slope')
    T_flat=Ts; 
else T_flat=0;
end


%Calculate incoming solar flux, minus effects from terrain since that must
%be calculated based on flat floor temperatures
F_solar = SW_slope; 

F_VIS_terr = (A_flat.*SW).*(1-E); 
 
F_IR_terr = eps_slope*sigma*T_flat.^4.*(1-E); 

f_vis_incoming =F_solar+ F_VIS_terr;
f_ir_incoming=  F_IR_terr;  

fi_vis_incoming=circshift(f_vis_incoming,[0,-1]);
fi_ir_incoming=circshift(f_ir_incoming,[0,-1]);


ae_guts=(dz(1)/(2.*k(1))).*((1-A_slope)* f_vis_incoming+eps_slope* f_ir_incoming);  
ai_guts=(dz(1)/(2.*k(1))).*((1-A_slope)*fi_vis_incoming+eps_slope*fi_ir_incoming); 
a_temp_prefix=(dz(1)/(2.*k(1))).*3.*eps_slope.*sigma; 


%Calculate temperatures for a sloping surface
for m=1:model_years
for n=1:n_tsteps
    
    b=1./(1+b_constants*Tr^3); 
    ke2(1)=1-(1-f)*(alpha_d(1)+(2-2*b)*beta);
    
    ae=(ae_guts(n)+a_temp_prefix*Tr^4).*b;
    ai=(ai_guts(n)+a_temp_prefix*Tr^4).*b;
    
    boundary(1)=2*beta*((1-f)*ae+f*ai); 
    
    T_present=(ke3.*circshift(T,[0,-1])+ke2.*T+ke1.*(circshift(T,[0,1])))+boundary;
    
    S(1,1)=1+f*(alpha_d(1)+(2-2*b)*beta);  %#ok<SPRIX> 
   
    T=S\(T_present');
    T=T';
 
    Tr=ai+b*T(1);
    Ts(n)=Tr;
    
    
end 
if m==initial_years
     T_average=mean(Ts);
     Tr=T_average;
     T=T_average*ones(1,n_layers); 
     Ts=zeros(1,max(size(n_tsteps)));
     
     
end

model_year=m 

end

Ts_slope_min=min(Ts)
Ts_slope_max=max(Ts)
Ts_slope_avg=mean(Ts)

T_lat_annual_avg(q)=mean(Ts);

M=[Ts]; 
filename = ['Ceres_latitudinal_annual_temps_for_' num2str(lat_array(q)*r2d) 'N_latitude_' num2str(A_slope) '_albedo' '.csv' ];

csvwrite(filename, M)

end
end


%IF slope =0, then save temperatures, emissivities, and albedos of the flat
%terrain case 

if slope == 0
    save('tempvars.mat', 'Ts', 'A_slope', 'eps_slope')
end

end

Ts_slope_min=min(Ts)
Ts_slope_max=max(Ts)
Ts_slope_avg=mean(Ts)


%Vapor loss: these equations do not include a diffusive barrier to escaping
%water vapor because this is for exposed surface ice

molec_m=2.99151e-26; %molecular mass of water in kg

kb=1.38065e-23;      %Boltzmann's constant in Jules per Kelvin difference
gas_constant=1/(2*pi*kb); 
Po=611;             %reference pressure in Pa
Lf=51058.;
Tref=273.16;
R=8.31;        %universal gas constant Jules per mol per Kelvin
rho_ice=925; 

%hardcode in the Kuppers et al. 2014 observation for reference
Kuppers=6*ones(5e5, 1);

for q=1:max(size(lat_array))
   
    M= ['Ceres_latitudinal_annual_temps_for_' num2str(lat_array(q)*r2d) 'N_latitude_' num2str(A_slope) '_albedo' '.csv']
    Ts=csvread(M); 
    
    T=Ts;
     
    P_vap=Po.*exp((-Lf/R).*((1./T)-(1./Tref))); 

    dmdt(q,:)=sqrt((molec_m.*gas_constant)./T).*P_vap; 

    dhdt(q,:)=dmdt(q,:)/rho_ice; 

    h(q,:)=cumsum(dhdt(q,:))*dt; 
    
    q
    
    name2=['Ceres_flat_terrain_vapor_production' num2str(lat_array(q)*r2d) 'N' 'for' num2str(model_years) 'Cyr_' num2str(A_slope) '_albedo' '.mat'];
    save(name2, 'T', 'P_vap', 'dmdt', 'dhdt', 'h') 

end



Ls_plot=circshift(ls, [1, -(find(ls==min(ls))-1)])*r2d;
dmdt_plot=circshift(dmdt, [1, -(find(ls==min(ls))-1)])*1e6; 
Ts_plot=circshift(Ts, [1, -(find(ls==min(ls))-1)]); 


figure(2)
plot(Ls_plot, dmdt_plot, 'k')
xlabel('Ls' , 'FontSize', 20)
ylabel('Vapor production (kg km^{-2} s^{-1})', 'FontSize', 20)
axis([0 360 0 inf])
set(gca,'FontSize',20)
