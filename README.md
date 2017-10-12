# Ceres_vapor_paper_sourcecode

This is the source code used in <i>Landis et al.,</i> [2017] on the vapor production from various configurations of water ice on Ceres (DOI:10.1002/2017JE005335)


<b>Buried ice table sublimation</b>

<i>"Ceres_annual_average_temperatures.m"</i> should be run first. It calculates the annual average temperatures with latitude over Ceres, assuming it is a smooth sphere. This must be run first to give temperature files used later. The default is to calculate annual average temperatures for 2-20 obliquity. 

<i>"Ceres_obliquity_variation_regolith_buildup.m"</i> calculates depth of sublimation lag for a pore-filling water ice table using the varying obliquity described in the literature (<i>Ermakov et al.,</i> [2017]). 

<i>"Ceres_constant_obliquity_vapor_production.m"</i> loads one of the annual average temperature files generated from "Ceres_annual_average_temperatures.m" and calculates the vapor and depth of sublimation lag produced on Ceres over solar system history. Can handle pore-filling and excess water ice cases. 

<b>Exposed surface ice sublimation</b>

<i>"Flat_and_intra_crater_terrain_vapor_production.m" </i> calculates sublimation for exposed surface ice (e.g. ice without a dessicated regolith overlayer). The code can be set to model flat terrain or terrain within an impact crater (impact crater parameters currently set to those of Oxo crater, <i>Combe et al.,</i> [2016])

<i>"Ceres_brightspot_fading_code.m"</i> determines the fading lifetime of water ice patch on Ceres (depending on latitude and assumed dust content) by having a monolayer of regolith particles build up over the surface. 

<b>Particle lofting</b>

<i>"finite_area_particle_lofting.m"</i> plots the particle diameter lofted at certain heights at a given temperature above a water ice patch on Ceres. 


