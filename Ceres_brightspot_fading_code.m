%Ceres_brightspot_fading.m--a script to calculate the fading lifetimes of
%bright spots on Ceres. Modern Ceres orbital parameters are used, but can
%be changed. First calculates a column of ice with a many thermal skin
%depts thick layer of regolith on the top, and then in a second time loop
%takes this away and allows that bright area to sublimate. Time,
%temperatures, vapor pressure, change in mass and linear dimension of an
%ice table are also included.

%Written by M.E. Landis, used in M.E. Landis et al., 2017 JGR
%Created in MATLAB_2016a


clc
clear all

format long e
warning('off', 'MATLAB:circshift:ScalarShift')


f=0.5; %0/1 for fully explicit/implicit, 0.5 for Crank-Nicholson

%degrees to radians conversions
d2r=pi/180;                             %convert degrees to radians
r2d=180/pi;                             %convert radians to degrees

%Latitude and longitude
lat=0          *pi/180;         %north latitude, in degrees

dust_fraction=0.05;             %volume fraction of regolith (so 1-c is the volume fraction of ice...)

r=50e-6;    %set regolith particle radius

%Slope dependence
slope=0         *d2r;
az_slope=0      *d2r;

E=cos(slope/2)*cos(slope/2); %Calculate the sky view factor, E=1 for flat terrain, E=cos^2(slope/2) for a sloped surface on a flat plane, see Shane's notes for the general
%formula. There is a seperate formula for when this ice patch is put in a
%crater. 

%Set up layers and assign thermal properties/constants to variables 

n_layers=300;     %number of layers in the model

TI_reg=15;        %thermal intertia of the regolith, assuming standard asteriod regolith properties 
rho_reg=1388.0;
c_reg=837.0;

TI_ice=2100;
A_slope=0.135;            %Albedo for the bright spot in Oxo crater as measured by Combe et al (2016)
eps_slope=0.95;          %emissivity for ice, which is what Fanale and Salvail (1989) used  
sigma=5.67e-8;  
rho_ice=925.0; 
c_ice= 1615.0; 

Q=0.002;                    %internal heat flux of Ceres 
 
T0=150;                     %temperature to initialize the initial run of the thermal model  

f1au=1367;                              %solar flux at 1AU


%Info on CERES orbit in the current epoch (this would change for a previous epoch)

a=2.76750591440571;                %Orbital semi-major axis
e=0.07582276595896797;             %Orbital eccentricity
lsp=302.11022    *d2r;                    %Ls of perihelion in DEGREES
obl=4.03        * d2r;                            %Obliquity in DEGREES
day_length=32667;
                   
au = 1.4959787061e11;
gc  = 6.67259e-11;               % Gravitational constant
sm  = 1.9891e30;                 % Solar mass
yr = 2*pi/sqrt(gc*sm/(a*au)^3);

dt=600; %min(sf.*(rho.*c.*dz.*dz.*0.5.*inv_k));  %time step, can be arbitrary or set according to the Courant criterion (which is in the comments but have to change these parameters to the new names they have above)

n_tsteps=round(yr/dt); 
time=linspace(0,yr,n_tsteps); 
ls=linspace(0,2*pi, n_tsteps); 

HA_corr=mod((time/(day_length)), 1);  %Correct times so model years begin at the same HA/lst
tsteps=yr;
tsteps=tsteps+((1-HA_corr(max(size(time)))*day_length)); 
time=linspace(0,yr,n_tsteps); 
n_tsteps=max(size(time)); 


%Set up parameters of the initial model run (this is the one with the
%regoltih layer over top of the ice column)


model_years=15; %how many years for the model to be run for the initial determination of temperatuers within the ice column 


%Set up arrays of thermal properties for the column which in this case is
%dust on top and ice all the way down 
rho=rho_ice*ones(1,n_layers);        
rho(1)=rho_reg; 
c=c_ice*ones(1,n_layers); 
c(1)=c_reg;
k=(TI_ice.^2/(rho(2)*c(2)))*ones(1,n_layers);
k(1)=(TI_reg.^2/(rho(1)*c(1))); 


%set reference temperatures and initialize temperature arrays
Tr=150;
Ts=150;

Ts_plot=zeros(1,n_tsteps);
T=Tr*ones(1,n_layers);

%Set up depths

dz_0=0.5*sqrt((k(2)*day_length)/(pi*rho(2)*c(2))); %70.79e-6; %varying dz
dz=dz_0*1.03.^(linspace(0,n_layers-1,n_layers)); 


%Setting time steps using the calculation of a year from where the
%planetary body is in the solar system, complete with true anamoly and
%orbital mechanics things 

sf=0.1;  %Stability factor, the fraction of the maximum stable time step set (max. time step calculated from Courant criterion)
         %A jump will occur on highly inclined surfaces near sunrise/sunset
         %for ~1 time step because of the way SW is calculated. Be careful
         %on slopes greater than about 20 degrees because the time steps
         %must be smaller=>also running into issues with the sun rising
         %instantly. For now this okay for large scale models, but be
         %careful!! If weird stuff appears this might be why. 

k3l=sqrt(gc*sm/((a*au)^3)); 
year = 2*pi/sqrt(gc*sm/(a*au)^3);                   %length of a year in seconds

% Low-Res calculation of times at uniformly spaced true anomalies
n=10000; 

ta=(2*pi/n)*(linspace(0, n, n));    %create evenly spaced true anomalies
ea_inners= (e+cos(ta))./(1+e*cos(ta));  %calculate eccentric anomalies
ea=acos(ea_inners);
ma=ea-e*sin(ea);                        %calculate mean anomalies
t=ma/k3l;                               %time along orbital path (will be irregularly spaced)

t(ta>pi)=yr-t(ta>pi);                      %correction of 2nd half of the orbit
 

%High-Res interpolation of true anomalies at uniformly spaced times

t2  = linspace(0,yr, n_tsteps); %Using uniformly spaced times to calculate eccentric anamolies
ta2= interp1(t,ta, t2); 
ea2_inners= (e+cos(ta2))./(1+e*cos(ta2)); %calculate the eccentric anomalies at these times 
ea2=acos(ea2_inners);  

% Solar distance/declination/longitude and local time calculations

sol_dist=a*(1-e*cos(ea2));              %solar distance
ls =mod(ta2+lsp, 2*pi);                 %solar longitude

sin_dec=sin(obl)*sin(ls);              %sine of the solar declination
cos_dec=sqrt(1-sin_dec.^2); 
 
%initialize temperature arrays (T, Ts)
Ts=zeros(1,n_tsteps);
T=Tr*ones(1,n_layers);

%incidence angles and fluxes/solar radiance parameters (with slopes)

HA=mod((time/(day_length)*2*pi), 2*pi)-pi;  %starts at dawn 

cosi=(sin(lat).*sin_dec+cos(lat).*cos_dec.*cos(HA)); 
sini=sqrt(1.0-cosi.^2);

cosi(cosi<0)=0; 

az_star_start=((sin_dec)-sin(lat)*cosi)./(cos(lat)*sini); %solar azimuth
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


cosi_noon=(sin(lat).*sin_dec+cos(lat).*cos_dec); %Calculate maximum noontime flux for those basic atmospheric models--this will be turned off in the Ceres model but is retained because who knows when I'll have to cannibalize this for Mars 
SW_noon=(f1au./(sol_dist.^2)).*cosi_noon; 


%Constants for both explicit and implicit versions of the model for the
%dust covered ice stability temperature determination 
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
    
             ae=(dz(1)/(2.*k(1))).*((1-A_slope)*f_vis_incoming(n)+eps_slope*f_ir_incoming(n)+3.*eps_slope.*sigma*Tr^4).*b;
             ai=(dz(1)/(2.*k(1))).*((1-A_slope)*fi_vis_incoming(n)+eps_slope*fi_ir_incoming(n)+3.*eps_slope.*sigma*Tr^4).*b;
    
            boundary(1)=2*beta*((1-f)*ae+f*ai); 
    
            T_present=(ke3.*circshift(T,[0,-1])+ke2.*T+ke1.*(circshift(T,[0,1])))+boundary;
    
            S(1,1)=1+f*(alpha_d(1)+(2-2*b)*beta);  %#ok<SPRIX> 
   
            T=S\(T_present');
            T=T';
 
            Tr=ai+b*T(1);
            Ts(n)=Tr; 
        end
        m
        
end

T_average=mean(Ts);

name1=['An_avg_temp_' num2str(lat*r2d) 'N' '.mat']; 

save(name1, 'T_average')


load(name1, 'T_average')

Tr=T_average;
T=T_average*ones(1,n_layers); 
Ts=zeros(1,max(size(n_tsteps)));

%And now, by some geological miracle, the top layer of regolith is
%instantly removed! 


%Set up parameters for a totally icy column of material, in order to
%calculate fading lifetimes, etc. Due to the changing albedo, etc, more
%detail must be calculated in each case. 

model_years=1; %how many years for the model to be run for the initial determination of temperatuers within the ice column 

A_fade=0.09*ones(1,max(size(time))); 
A_fade(1)=A_slope; 
h=zeros(1,max(size(time)+1)); 

%Constants needed for calculation of vapor production rates and ice
%sublimation 

molec_m=2.99151e-26; %molecular mass of water in kg
kb=1.38065e-23;      %Boltzmann's constant in Jules per Kelvin difference
gas_constant=1/(2*pi*kb); 
Po=611;             %reference pressure in Pa
Lf=51058.;
Tref=273.16;
R=8.31;        %universal gas constant Jules per mol per Kelvin
rho_ice=925; 



X=[0 2*r];    %starting & ending point point of amount of sublimation lag in microns, assuming a monolayer of regolith particles
Y=[0.135 0.09]; %starting & ending points of albedo for the bright spot and then background albedo of Ceres


iceloss_for_dark_dust=(2*r/dust_fraction); %how much dust is in the ice 

Albedo_change_constants=(A_slope-0.09)/iceloss_for_dark_dust; 

eyr=3.1558149e7;      


%Set up arrays of thermal properties for the column which in this case is
%dust on top and ice all the way down 
rho=rho_ice*ones(1,n_layers);        
rho(1)=rho_ice; 
c=c_ice*ones(1,n_layers); 
c(1)=c_ice;
k=(TI_ice.^2/(rho(2)*c(2)))*ones(1,n_layers);
k(1)=(TI_ice.^2/(rho(1)*c(1))); 


%reference temperatures have already been rest to the annual average temp,
%which is the most likely temp for all the ice layers 

%Set up depths

dz_0=0.5*sqrt((k(2)*day_length)/(pi*rho(2)*c(2))); %70.79e-6; %varying dz   %actually just cut off the first couple layers from the previous column
dz=dz_0*1.03.^(linspace(0,n_layers-1,n_layers)); 


%Constants for both explicit and implicit versions of the model for the
%dust covered ice stability temperature determination 
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


%This loop will calculate the fading and albedo change for one year. This
%first year will be used to estimate the total fading time used in the loop
%after this. 

 for n=1:n_tsteps    
     
        F_VIS_terr = (A_fade(n).*SW_slope).*(1-E); 
        f_vis_incoming =F_solar+F_VIS_terr;
        fi_vis_incoming=circshift(f_vis_incoming,[0,-1]);
    
        b=1./(1+b_constants*Tr^3); 
        ke2(1)=1-(1-f)*(alpha_d(1)+(2-2*b)*beta);
    
        ae=(dz(1)/(2.*k(1))).*((1-A_fade(n))*f_vis_incoming(n)+eps_slope*f_ir_incoming(n)+3.*eps_slope.*sigma*Tr^4).*b;
        ai=(dz(1)/(2.*k(1))).*((1-A_fade(n))*fi_vis_incoming(n)+eps_slope*fi_ir_incoming(n)+3.*eps_slope.*sigma*Tr^4).*b;
    
        boundary(1)=2*beta*((1-f)*ae+f*ai); 
    
        T_present=(ke3.*circshift(T,[0,-1])+ke2.*T+ke1.*(circshift(T,[0,1])))+boundary;
    
        S(1,1)=1+f*(alpha_d(1)+(2-2*b)*beta);  %#ok<SPRIX> 
   
        T=S\(T_present');
        T=T';
        T_check(:,n)=T; 
 
        Tr=ai+b*T(1);
       
        P_vap(:,n)=Po.*exp((-Lf/R).*((1./Tr)-(1./Tref))); 

        dmdt(:,n)=sqrt((molec_m.*gas_constant)./Tr).*P_vap(:,n);     %this needs to be multiplied by 10^6 to get the actual mass in kg 
        dhdt(:, n)=dmdt(n)/rho_ice; 
        h(:,n+1)=sum(dhdt*dt); 
        
        T_fade(:,n)=Tr;
       
        if h(n+1)<iceloss_for_dark_dust
            A_fade(n+1)=A_slope-Albedo_change_constants*h(:,n+1); 
        
            
        end
        
 end

 
name2=['Fading_vars_for_' num2str(lat*r2d) 'N' '.mat'];  

save(name2, 'T_fade', 'A_fade','h', 'dmdt', 'dhdt'); 

load(name2, 'T_fade', 'A_fade','h', 'dmdt', 'dhdt'); 


if (max(dmdt)*eyr*1e9)/rho_ice <=1     %this is the 1m/Gyr ice loss cut off, and will stop the code for running for forever if ice is stable at that latitude
 disp(' ICE IS STABLE!! stop the code')
 save_name=['ICE_STABLE_FOR' num2str(lat*r2d) 'N' '.mat']; 
 save(save_name)
end



if h(end)>iceloss_for_dark_dust    %if one year is enough to sublimate enough dust to fade the spot, this loop will be triggered
   fading_time_seconds=time(min(find(h>iceloss_for_dark_dust)));  %#ok<NASGU>
   disp ('Fading time (in seconds)')
   disp (fading_time_seconds)   
end




if h(end)<iceloss_for_dark_dust   %if the spot isn't dark after one year, then run the code for m years to get it to fade. There is a "break" statement for when the spot fades, though, so the code doesn't run forever
 
h_end=h(end); 

round(iceloss_for_dark_dust/h_end)+1
 
    for m=1:round(iceloss_for_dark_dust/h_end)+1

        A_start=A_fade(end); 
        A_fade=0.09*ones(1,max(size(time))); 
        A_fade(1)=A_start; 
 
        h_start=h(end);    
        h=zeros(1,max(size(time)+1)); 
        h(1)=h_start; 
        
        for n=1:n_tsteps    


            F_VIS_terr = (A_fade(n).*SW_slope).*(1-E); 
            f_vis_incoming =F_solar+ F_VIS_terr; 
            fi_vis_incoming=circshift(f_vis_incoming,[0,-1]);

            b=1./(1+b_constants*Tr^3); 
            ke2(1)=1-(1-f)*(alpha_d(1)+(2-2*b)*beta);

            ae=(dz(1)/(2.*k(1))).*((1-A_fade(n))*f_vis_incoming(n)+eps_slope*f_ir_incoming(n)+3.*eps_slope.*sigma*Tr^4).*b;
            ai=(dz(1)/(2.*k(1))).*((1-A_fade(n))*fi_vis_incoming(n)+eps_slope*fi_ir_incoming(n)+3.*eps_slope.*sigma*Tr^4).*b;

            boundary(1)=2*beta*((1-f)*ae+f*ai); 

            T_present=(ke3.*circshift(T,[0,-1])+ke2.*T+ke1.*(circshift(T,[0,1])))+boundary;

            S(1,1)=1+f*(alpha_d(1)+(2-2*b)*beta);  %#ok<SPRIX> 

            T=S\(T_present');
            T=T';
            T_check(:,n)=T; 

            Tr=ai+b*T(1);

            P_vap(:,n)=Po.*exp((-Lf/R).*((1./Tr)-(1./Tref))); 

            dmdt(:,n)=sqrt((molec_m.*gas_constant)./Tr).*P_vap(:,n);  
            dhdt(:, n)=dmdt(n)/rho_ice; 
            h(:,n+1)=sum(dhdt*dt); 

            T_fade(:,n)=Tr;

            if h(n+1)<iceloss_for_dark_dust
                A_fade(n+1)=A_slope-Albedo_change_constants*h(:,n+1); 
            else A_fade(n+1)=0.09; 

            end

        end
      m  
      h(end)
      if h(end)>iceloss_for_dark_dust
          t_fade=m*time(end)+time(min(find(h>iceloss_for_dark_dust)))
          break
      end
    end
    
    save_name=['Fade_lifetime_for' num2str(lat*r2d) 'N' '.mat']; 
    save(save_name, 't_fade')
end



