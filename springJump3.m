%calculates the jump height with effective spring energy density
us=1922; %energy density of springs with structure
ps=1000;%density of springs
g=9.81; %gravitational acceleration


%estimated drag parameters
Cd=1;%estimated drag parameters
rhoAir=1.225;%kg/m^3

ms=[1 1 1 1 1]; %relative spring mass 
mm=[1000 100 10 1 .1];%relative motor mass to spring
mu=ms./(ms+mm); %for ms+mm=.14;

d=logspace(-7.5,0,1000); 
dMat=d'*(2*us*mu./(g*(2-mu))); %scale
muMat=ones(length(d),1)*mu;

tto=sqrt((3-2*muMat)./(6*us*muMat)).*real(acos((muMat-2)*g.*dMat./((muMat-2).*g.*dMat+4*us.*muMat))).*dMat; %for some reason it gives imaginary numbers
vto=real(sqrt(3*((muMat-2).*g.*dMat+2*muMat*us)./(3-2*muMat)));
w=us.*muMat;
vCOM=vto./dMat.*(dMat/2.*muMat+dMat.*(1-muMat));%COM velocity
Tcom=1/2*vCOM.^2;%com ke


%COM KE
figure
colororder([0 255 255;0 191 255;0 128 255;0 64 255 ;0 32 122]/255) 

loglog(dMat,Tcom,'LineWidth',2)
hold on
loglog([1e-6 1e4],[us us],'k--','LineWidth',2)
xlabel('Scale (m)')
ylabel('Specific Energy (J/kg)')
xlim([1e-4 1e3])
title('COM Kinetic Energy')

figure
colororder([0 255 255;0 191 255;0 128 255;0 64 255 ;0 32 122]/255) 
loglog(dMat,w,'LineWidth',2)
hold on
loglog([1e-6 1e4],[us us],'k--','LineWidth',2)
xlabel('Scale (m)')
ylabel('Specific Energy (J/kg)')
xlim([1e-4 1e3])
title('Produced Energy')
% legend('\mu='+string(mu),'Location','southeast')
% hleg=legend(string(mu/.14),'Location','southeast'); %extended data 4 legend
hleg=legend(string(ms./mm),'Location','southeast'); %extended data 4 legend for ms+mm=.14

htitle=get(hleg,'Title');
set(htitle,'String','Spring-motor mass ratio')
hleg.NumColumns=2;


figure
colororder([0 255 255;0 191 255;0 128 255;0 64 255 ;0 32 122]/255) 
loglog(dMat,tto,'LineWidth',2)
xlabel('Scale (m)')
ylabel('Time (s)')
xlim([1e-4 1e3])
title('Acceleration Time')

% legend('\mu='+string(mu),'Location','southeast')
