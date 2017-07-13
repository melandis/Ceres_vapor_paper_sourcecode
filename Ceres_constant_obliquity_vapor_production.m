%Ceres_constant_obliquity_vapor_production.m a script to calculate the growth of a sublimation lag
%over ice on Ceres, either pore-filling or massive depending on the value
%of "C", regolith content by volume

%Written by M.E. Landis, used in M.E. Landis et al., 2017 JGR
%Created in MATLAB_2016a

clc
clear all

%import temperature array from previous code
Vars=csvread('Ceres_latitudinal_annual_avg_temps4obliquity0slope0azimuth.csv'); 
Lat1=Vars(:,1); %min latitude of band
Lat2=Vars(:,2); %max latitude of band
A_msqr=Vars(:,3); %area of latitude band
T=Vars(:,4);       %annual average temperature for the latitude band
inv_T=1./T;
inv_A=1./A_msqr;

%Model parameters
model_years=4700000000;     %how many EARTH YEARS the model will run for 
h0=3e-2;              %initial regolith thickness


%constants
cday=32680; %length of a Ceres day in seconds
cyr=1.451521976718465e+08; %length of a Ceres year in seconds
eyr=3.1558149e7;            %length of an Earth year in seconds        
t_s=eyr*model_years;       %amount of time to run model in seconds
phi=0.5;             %Porosity
inv_phi=1./phi; 
tau=2;               %Tortuosity 
inv_tau=1./tau; 
r=50e-6;            %grain radius
A=0.09;              %Albedo
TI=15;               %Thermal inertia
molec_m=2.99151e-26; %molecular mass of water in kg
c=[0, 7e-5, 0.02 0.25, 0.5];             %volume fraction of regolith (so 1-c is the volume fraction of ice...)
                        % c=0 is a pure water ice table,, >0.5c>0 is excess ice, >0.5 is pore
                        % filling ice
kb=1.38065e-23;      %Boltzmann's constant in Jules per Kelvin difference
gas_constant=1/(2*pi*kb); 
Po=611;              %reference pressure in Pa
Lf=51058.;
inv_Lf=1/Lf;
Tref=273.16;
inv_Tref=1./Tref; 
R=8.31;        %universal gas constant Jules per mol per Kelvin
rho_ice=952; 
inv_rhoice=1/rho_ice; 

time=linspace(100, t_s, 5e5);

%hardcode in the Kuppers et al. 2014 observation for reference
Kuppers=6*ones(5e5, 1);

%general code that will do either massive or pore filling ice cases
for m=1:max(size(c))

P_vap=Po*exp((-Lf/R)*(inv_T-inv_Tref)); 

J_constants=((4*pi)/(8+pi)).*(phi/(1-phi)).*inv_tau.*r.*sqrt(molec_m.*gas_constant*inv_T).*P_vap; 

h=sqrt(h0^2+2*(J_constants*inv_rhoice*(1/(1-phi))*(c(m)/(1-c(m)))*time));

for n=1:max(size(h))
   J(:,n)=(J_constants)./h(:,n);
   J_area(:,n)=A_msqr.*J(:,n);
   tot_kg_lost(:,n)=sum(J_area(:,n));
   
end
    J_save(:,:,m)=J; 
    h_save(:,:,m)=h; 
    J_area_save(:,:,m)=J_area;
    J_tot_ceres(:,m)=sum(J_area);
end

time_plot=time/eyr; 
kg_lost=sum(J_area_save(:,1:1:95745,:), 2);
kg_lost_per_area=sum(J_save(:,1:1:95745,:), 2);


figure(105)
loglog(time_plot, J_tot_ceres(:,1), 'r', 'LineWidth', 3)
hold on 
loglog(time_plot, J_tot_ceres(:,2), 'b', 'LineWidth', 3)
loglog(time_plot, J_tot_ceres(:,3), 'Color',1/169.6879*[23,128,109], 'LineWidth', 3)
loglog(time_plot, J_tot_ceres(:,5), 'k', 'LineWidth', 3)
loglog(time_plot, Kuppers, 'k-.', 'LineWidth', 3)
ylabel('Vapor output (whole Ceres), kg*s^{-1}', 'FontSize', 24)
xlabel('time (years)', 'FontSize', 24)
set(gca,'fontsize',20)
leg1=sprintf('c=%6.5f' ,c(1));
leg2=sprintf('c=%6.5f',c(2));
leg3=sprintf('c=%4.2f',c(3));
%leg4=sprintf('c=%4.2f',c(4));
leg5=sprintf('c=%4.2f',c(5));
%leg6=sprintf('c=%4.2f',c(6));
legend(leg1, leg2, leg3, leg5, 'Kuppers et al observation', 'Location', 'best') %, leg3, leg4, leg5, leg6,
axis([1e6 7e9 0.03 100])

figure(5)
hold on
plot((Lat1+Lat2)/2, h_save(:,478724,1), 'r', 'LineWidth', 3)
plot((Lat1+Lat2)/2, h_save(:,478724,2), 'b', 'LineWidth', 3)
plot((Lat1+Lat2)/2, h_save(:,478724,3), 'LineWidth', 3, 'Color', 1/169.6879*[23,128,109])
plot((Lat1+Lat2)/2, h_save(:,95745,5), 'k', 'LineWidth', 3)
set(gca,'Ydir','reverse')
set(gca,'Xdir','reverse')
set(gca,'fontsize',20)
xlabel('N Latitude', 'FontSize', 24)
ylabel('Regolith thickness (m)', 'FontSize', 24)
axis([-90 90 -inf 3])

