clear; close all; clc;

EjCarbon = 1928; 		%J/kg 	specific energy of carbon fiber legs
carbonMass = 0.0125; 	%kg  		mass of carbon fiber legs
payloadMass = 0.007; 	%kg		payload mass (assumed to be at top of body)
motorMass = 0.0103; 	%kg		motor mass (assumed to be at top of body)
footMass = 0.001; 		%kg 		approximate mass of items at base of legs (hinges)
totalMass = carbonMass+payloadMass+motorMass+footMass;
g = 9.8; 			%m/s^2 	gravity
A = 0.018^2 ;			%m^2  	frontal area
Cd = 1;			%		coefficient of drag
rho = 1.225;			% kg/m^3	density of air


eff = (.63+((payloadMass+motorMass)/totalMass-.5)*.53)*.888*.94 *.92 
	%efficiency; this is an approximation based off a simulation) =>Use this for bow jumpers
%eff = (payloadMass+motorMass+carbonMass/2)/(totalMass) 
	%approximation for simple jumpers=> use this for linear jumpers, like a pogo stick

energy = EjCarbon*carbonMass; 
	%J		total energy
heightNoEffNoDrag = energy/totalMass /g; 
	%m		nominal jump height with no efficiency loss during jump and no air drag
heightNoDrag = heightNoEffNoDrag*eff	 
	%m		jump height with efficiency loss, but no drag loss
takeoffV = sqrt(g*heightNoDrag*2) 
	%m/s		takeoff velocity
Ds = Cd*rho*A*takeoffV^2/(2*totalMass*g); 
	% drag term
height = heightNoDrag .* 1./Ds.*log(1+Ds/1) 
	%m		jump height with efficiency loss and drag loss
