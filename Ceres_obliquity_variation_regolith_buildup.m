%Ceres_obliquity_variation_ice_retreat.m 

%A script to give the sine wave for A. Ermkov et al. (2017)'s solution to Ceres' past obliqity
%to put into the ice retreat code. 

%Written by M.E. Landis, used in M.E. Landis et al., 2017 JGR
%Created in MATLAB_2016a


%Obliquity is basically a sine wave with a period of 25kyr and max 20
%degrees and min 2 degrees. In this chunk of code, that sine wave equation
%is calculated and then each fraction of 25kyr is spent at each obliquity.
%This is going to provide the detailed modeling for each obliquity cycle
%for Ceres (since obliquity is moving quickly compared to the lifetime of
%the solar system)

clc 
clear all

spy=3.1558149e7;               %seconds per terrestrial year 
period_of_obliquity_cycle=25e3;  %how many terrestrial years does it take for one obliquity cycle to occur?

a=(20-2)/2;                       %amplitude of the change in obliquity
f=1/(period_of_obliquity_cycle*spy);                     %frequency 

n_steps=1e4; 

time=linspace(0, 25e3*spy, n_steps); 
dt=time(2)-time(1); 

obl=a*sin(2*pi*f*time+(3*pi/2))+12;         %this starts the obliquity cycle at the 4 (modern Ceres) obliquity 

obl_step_f=round(obl, 0); 


for n=4:20
    number=sum(obl_step_f==n); 
    array(n-3)=number; 
    
end

time_at_different_obliquities=array*dt; %how much time is spent at each obliquity. Most of
                                         %the time is spent at either 4 or
                                         %20 degrees obliquity
                                        

% This part of the code reassembles a series of CSV files with annual
% average temperatures at different latitudes from a previous model run so
% that there is a menu of obliquity/temperature combinations 


Vars=csvread('Ceres_latitudinal_annual_avg_temps2obliquity0slope0azimuth.csv');
Lat1=Vars(:,1); %min latitude of band
Lat2=Vars(:,2); %max latitude of band
A_msqr=Vars(:,3); %area of latitude band

for n=2:20
    filename=['Ceres_latitudinal_annual_avg_temps' num2str(n) 'obliquity0slope0azimuth.csv']; 
    Vars=csvread(filename); 
    T(:,n-1)=Vars(:,4);   
end

                                       
%Set up variables that will be used in the vapor diffusion equations based
%on Schorghofer, 2008

%constants
h_start=3e-2*ones(37,1);              %initial regolith thickness  
phi=0.5;             %Porosity
inv_phi=1./phi; 
tau=2;               %Tortuosity 
inv_tau=1./tau; 
r=50e-6;             %Pore size
molec_m=2.99151e-26; %molecular mass of water in kg
c=[0.5];             %volume fraction of regolith (so 1-c is the volume fraction of ice...)
                        % c=0 is a pure water ice table, >0.5c>0 is excess ice, >0.5 is pore
                        % filling ice
kb=1.38065e-23;      %Boltzmann's constant in Jules per Kelvin difference
gas_constant=1/(2*pi*kb); 
Po=611;             %reference pressure in Pa
Lf=51058.;
Tref=273.16;
inv_Tref=1./Tref; 
R=8.31;        %universal gas constant Jules per mol per Kelvin
rho_ice=952; 
inv_rhoice=1/rho_ice; 

model_years=[500e6 1e9 2e9 3e9 4e9]; %how long to run the regolith generation/vapor production model

%hardcode in the Kuppers et al. 2014 observation for reference
Kuppers=6*ones(5e5, 1);


%general code that will do either massive or pore filling ice cases for
%multiple Ceres obliquity cycles 

%h is the thickness of regolith above the ice 

inv_T=1./T;
P_vap=Po*exp((-Lf/R)*(inv_T-inv_Tref));
J_constants=((4*pi)/(8+pi)).*(phi/(1-phi)).*inv_tau.*r.*sqrt(molec_m.*gas_constant*inv_T).*P_vap;


model_year=max(model_years); 
    
for obl_cycles=1:model_year/period_of_obliquity_cycle; 
    
for t=1:max(size(time))

obl_tstep=obl_step_f(t);   

h=sqrt(h_start.^2+2*(J_constants(:,obl_tstep-2) *inv_rhoice*(1/(1-phi))*(c/(1-c))*dt));

h_start=h; 

end

h_save(:, obl_cycles)=h; 
h_end(:, obl_cycles+1)=h; 
h_start=h_end(:,obl_cycles+1); 

if mod(obl_cycles, 100)==0
    obl_cycles
end


if obl_cycles==(model_years(1)/period_of_obliquity_cycle)
    filename=['Ceres_full_obliquity_var_run_for_50_micron_particles' num2str(model_years(1)/period_of_obliquity_cycle) '_obliquity_cycles.mat'] 
    save(filename, 'h_save');     
end

if obl_cycles==(model_years(2)/period_of_obliquity_cycle)
    filename=['Ceres_full_obliquity_var_run_for_50_micron_particles' num2str(model_years(2)/period_of_obliquity_cycle) '_obliquity_cycles.mat'] 
    save(filename, 'h_save');     
end

if obl_cycles==(model_years(3)/period_of_obliquity_cycle)
    filename=['Ceres_full_obliquity_var_run_for_50_micron_particles' num2str(model_years(3)/period_of_obliquity_cycle) '_obliquity_cycles.mat'] 
    save(filename, 'h_save');     
end

if obl_cycles==(model_years(4)/period_of_obliquity_cycle)
    filename=['Ceres_full_obliquity_var_run_for_50_micron_particles' num2str(model_years(4)/period_of_obliquity_cycle) '_obliquity_cycles.mat'] 
    save(filename, 'h_save');     
end

if obl_cycles==(model_years(5)/period_of_obliquity_cycle)
    filename=['Ceres_full_obliquity_var_run_for_50_micron_particles' num2str(model_years(5)/period_of_obliquity_cycle) '_obliquity_cycles.mat'] 
    save(filename, 'h_save');     
end
 
end


    