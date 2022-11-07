%this script estimates the spring motor mass ratio for extended fig. 6

us=1950; %Effective energy density of springs with structure
ps=1000;%density of springs
g=9.81; %gravitational acceleration
hE=32.9; %experimetnal jump height

%the robot weighs 30g
mr=.03; %[kg] total mass of the robot
ms=.0124; %weight of the spring on the robot
mm=.01033; %mass of the motor gear assembly
md=mr-ms-mm; %amount of dead mass on the robot
l=.3; %[m] length of the robot in meters
A=.018^2; %estimated frontal area of the robot 

%estimated drag parameters
Cd=1;%estimated drag parameters
rhoAir=1.225;%kg/m^3
rho1=mr/(A*l);%1000;%kg/m^3 
vtoE=real(sqrt(3*((ms/mr-2).*g.*l+2*ms/mr*us)./(3-2*ms/mr)));%theoretical top mass velocity based on experimental parameters
vCOME=vtoE./l.*(l/2.*ms/mr+l.*(1-ms/mr));%theoritical COM velocity
eta=fzero(@(eta) hE-rho1*l/(Cd*rhoAir)*log(Cd*rhoAir/(rho1*l*g)*.5*vCOME^2*eta+1),.9);%fitted efficiency parameter


%vary mu of the robot to show the effect of spring mass jumper
mmVec=[2*ms ms/1.2 ms/10 0];
m=mmVec+ms+md; 

rho=m/(.018^2*.3);%1000;%kg/m^3 
rhoMat=ones(1000,1)*rho;

mu=ms./m; %for ms+mm=.14;

d=logspace(-7.5,0,1000); 
dMat=d'*(2*us*mu./(g*(2-mu))); %scale
muMat=ones(length(d),1)*mu;

tto=sqrt((3-2*muMat)./(6*us*muMat)).*real(acos((muMat-2)*g.*dMat./((muMat-2).*g.*dMat+4*us.*muMat))).*dMat; %for some reason it gives imaginary numbers
vto=real(sqrt(3*((muMat-2).*g.*dMat+2*muMat*us)./(3-2*muMat)));
w=us.*muMat;
vCOM=vto./dMat.*(dMat/2.*muMat+dMat.*(1-muMat));%COM velocity
Tcom=1/2*vCOM.^2;%com ke
h=rhoMat.*dMat./(Cd*rhoAir).*log(eta*Cd*rhoAir*vCOM.^2./(2*rhoMat*g.*dMat)+1);


%COM KE
figure
colororder([0 255 255;0 191 255;0 128 255;0 64 255]/255) 

loglog(dMat,Tcom,'LineWidth',2)
hold on
loglog([1e-6 1e4],[us us],'k--','LineWidth',2)
xlabel('Scale (m)')
ylabel('Specific Energy (J/kg)')
xlim([1e-4 1e1])


title('COM Kinetic Energy')

grid on

figure
colororder([0 255 255;0 191 255;0 128 255;0 64 255]/255) 
loglog(dMat,w,'LineWidth',2)
hold on
loglog([1e-6 1e4],[us us],'k--','LineWidth',2)
xlabel('Scale (m)')
ylabel('Specific Energy (J/kg)')
xlim([1e-4 1e1])
title('Produced Energy')
% legend('\mu='+string(mu),'Location','southeast')
% hleg=legend(string(mu/.14),'Location','southeast'); %extended data 4 legend
hleg=legend(string(ms./mmVec),'Location','southeast'); %extended data 4 legend for ms+mm=.14

htitle=get(hleg,'Title');
set(htitle,'String','Spring-motor mass ratio')
hleg.NumColumns=2;
grid on
grid on


figure
colororder([0 255 255;0 191 255;0 128 255;0 64 255]/255) 
loglog(dMat,h,'LineWidth',2)
hold on
loglog([1e-6 1e4],[us/g us/g],'k--','LineWidth',2)
plot(.3,32.9,'p','Color',[0 0 0])
xlabel('Scale (m)')
ylabel('Height (m)')
xlim([1e-3 1e3])
ylim([hE max(max(h))])
title('Estimated Jump Height')
% hleg=legend(string(mu/(5/30)),'Location','southeast'); %motor mass ratio preserving, increaseing spring mass
hleg=legend(string(ms./mmVec),'Location','southeast'); %dead mass preserving legend

htitle=get(hleg,'Title');
set(htitle,'String','Spring-motor mass ratio')
hleg.NumColumns=2;
grid off
xlim([1e-3 1e2])
ylim([1 100])