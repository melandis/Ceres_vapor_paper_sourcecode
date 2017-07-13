%finite_area_particle_lofting.m 

%Code to plot the regolith particle size that can be lofted with height above the
%surface of a Ceres water ice patch

%Written by M.E. Landis, used in M.E. Landis et al., 2017 JGR
%Created in MATLAB_2016a

clc
clear all

%constants and presets 
kb=1.38065e-23;      %Boltzmann's constant in Jules per Kelvin difference
gas_constant=1/(2*pi*kb); 
Po=611;             %reference pressure in Pa
Lf=51058.;
Tref=273.16;
R=8.31;
density_carb_chondrite=2.5e3; 
molec_m=2.99151e-26; %molecular mass of water in kg

gc  = 6.67259e-11;
r_Ceres= 	473e3; 
m_Ceres= 8.958e20 ; 
g_Ceres= (gc*m_Ceres) / (r_Ceres^2); 

hot=180;
cold=155; 
T=linspace(cold, hot, 1000); 

v_thermal=sqrt( (3 * kb * T) / molec_m); 


%Start calculating parameters
%vapor pressure and initial diameter

P_vap=Po.*exp((-Lf/R).*((1./T)-(1./Tref))); 

d_o=(1.5*P_vap)./(density_carb_chondrite*g_Ceres); 


%finite ice sheet with radius r, theta max =/= pi/2, arctan(r/h)
%also decrease due to reduction in vapor velocity in height 


% r=5e3;            %set up for a ~5km radius circular water ice patch
% area=pi*r.^2/1e6; %area of ice patch in square km 

area=1e6;               %set up for a 1km^2 circular water ice patch
r=sqrt(area/pi); 


h=linspace(0, 10e3, 100); 

for n=1:max(size(h))
    d_finite_ice_sheet_redux(:,n)=[1-cos(atan(r./h(n)))]; 
    v_height_redux(:,n)=sqrt(1-((2.*g_Ceres.*h(n))./v_thermal.^2)); 

end

%calculate reduction in diameter due to these two effects 

for n=1:max(size(h))
    
    d_final(:,n)=d_o'.*(d_finite_ice_sheet_redux(:,n).*v_height_redux(:,n)); 

end

Contours=linspace(min(min(d_final)),max(max(d_final)),1e4); 

h_km=h/1e3; 

figure(1)
hold on 
imagesc(T, h_km, log10((d_final'*1e6)))
c=colorbar
c.Label.String='Particle size ({\mu}m)'
c.YTick=[-4 -3 -2 -1 0 1]
c.TickLabels={'10^{-4}','10^{-3}','10^{-2}', '10^{-1}', '10^0', '10'}    
colormap jet
xlabel ('Surface temperature (K)', 'FontSize', 30)
ylabel('Height above surface (km)', 'FontSize', 30)

set(gca, 'Ydir', 'normal')
set(gca,'fontsize',20)
axis([cold hot 0 5])

[C, h]=contour(T, h_km, log10((d_final'*1e6)), [-1, 0, 1], 'w', 'LineWidth', 3)




