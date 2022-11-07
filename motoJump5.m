%this does the motor jump for model that parameterizes specific energy and
%specific power. Also includes the effect of linkage mass. Assume no
%payload

g=9.81;
u=200; %[J/kg] specific energy bio
p=200;%200*.14; %[W/kg] specific power bio
ulb=1130;%[J/kg] energy transfer ability of biological linkages

u2=u*ulb/(ulb+u);
p2=p*ulb/(ulb+u);


pm=2000;%200*.14; %[W/kg] specific power for non animated things
ulm=2650; %[J/kg] linkage energy transfer ability for motorized systems 


%Drag parameters
Cd=1;%estimated drag parameters
rhoAir=1.225;%kg/m^3
rho=1000;%kg/m^3 

G=[.01 .1 1]; %gear ratios
Gs=G; %G*
Gs(Gs<1)=1;

d=logspace(-6,0,1000);
dMat=d'*G*u2/g; %scale for the gear ratios
GMat=ones(length(d),1)*G;
GsMat=ones(length(d),1)*Gs;



z=(16*p2^2*dMat.^2.*GsMat.*(dMat*g-GMat*u2)-GMat.^4*u2^4)./(16*dMat.^2*p2^2.*GsMat.*(dMat*g-GMat*u2)); %this is the exponent


tto=(16*z.*dMat.^2*p2^2+16*dMat.^2*p2^2.*lambertw(-exp(-z)))./(4*GMat.^2*u2^2*p2);
vto=(4*dMat*p2.*(GMat*u2-dMat*g).*(1+lambertw(-exp(-z))))./(GMat.^2*u2^2);
w=1/2*vto.^2+g*dMat./GsMat;
h=rho*dMat/(Cd*rhoAir).*log(Cd*rhoAir.*vto.^2./(2*rho*g*dMat)+1);


%this does optimal continually variable gear ratio for bio jumper.
%Nonvarying linkage mass
% tms=u2/(2*p2); %stroke length constraint for bio jumpers
% t=logspace(-6,log(tms)/log(10),1000);
% xO=-p2*(-2*g^2*t+p2*(1+lambertw(-exp(-1-g^2*t/p2)).*(2+lambertw(-exp(-1-g^2*t/p2)))))/(2*g^3);
% % xO2=[xO,logspace(log(xO(end))/log(10),log(xO(end)*1000)/log(10),1000)];%vector of scales
% % xO3=[xO,xO(end)*ones(1,1000)];%vector of displacements
% vO=p2/g*(1+lambertw(-exp(-1-g^2*t/p2)));%calculate sppeds
% % vO2=[vO,vO(end)*ones(1,1000)];
% wO=1/2*vO.^2+g*xO;
% hO=rho*xO/(Cd*rhoAir).*log(Cd*rhoAir.*vO.^2./(2*rho*g*xO)+1);

%max power bio jumper with varying linkage mass;
tms=u2/(2*p2); %stroke length constraint for bio jumpers. same value as varying linkage
t=logspace(-6,log(tms)/log(10),1000);
p3=p*ulb./(ulb+p*t);
xO=-p2*(-2*g^2*t+p2*(1+lambertw(-exp(-1-g^2*t/p2)).*(2+lambertw(-exp(-1-g^2*t/p2)))))/(2*g^3);
vO=p2/g*(1+lambertw(-exp(-1-g^2*t/p2)));%calculate sppeds
% vO2=[vO,vO(end)*ones(1,1000)];
wO=1/2*vO.^2+g*xO;
hO=rho*xO/(Cd*rhoAir).*log(Cd*rhoAir.*vO.^2./(2*rho*g*xO)+1);


%for the ideal motorized jumper
tim=logspace(-5,5,1000);
 xim=pm*(2*g^2*tim-pm*(1+lambertw(-exp(-1-g^2*tim/pm)).*(2+lambertw(-exp(-1-g^2*tim/pm)))))/(2*g^3);
% xim=-pm*(-2*g^2*tim+pm*(1+lambertw(-exp(-1-g^2*tim/pm)).*(2+lambertw(-exp(-1-g^2*tim/pm)))))/(2*g^3);
% vim=pm/g*(1+lambertw(-exp(-1-g^2*tim/pm)));
vim=pm/g*(1+lambertw(-exp(-1-g^2*tim/pm)));
wim=1/2*vim.^2+g*xim;
him=rho*xim/(Cd*rhoAir).*log(Cd*rhoAir.*vim.^2./(2*rho*g*xim)+1);

%this calculates the idealized motorized jumper as linkages scale
pm2=pm*ulm./(ulm+pm*tim);
xim2=pm2.*(2*g^2*tim-pm2.*(1+lambertw(-exp(-1-g^2*tim./pm2)).*(2+lambertw(-exp(-1-g^2*tim./pm2)))))/(2*g^3);

vim2=pm2./g.*(1+lambertw(-exp(-1-g^2*tim./pm2)));
wim2=1/2*vim2.^2+g*xim2;
him2=rho*xim2/(Cd*rhoAir).*log(Cd*rhoAir.*vim2.^2./(2*rho*g*xim2)+1);

%%
%make the plots
%take off time plot
figure
colororder([238 64 39;238 64 39;76 70 252;106 250 253; 0 115 252; 72 67 252]/255)
loglog(xim,tim,'LineWidth',2)
hold on
loglog(xim2,tim,'--','LineWidth',2)
loglog(xO,t,'LineWidth',2)

loglog(dMat,tto,':','LineWidth',2)

legend({'Eng., max power no linkage','Eng., max power opt. linkage','Bio., max power opt. linkage', 'Bio., G=0.01','Bio., G=0.1','Bio., G=1'},'Location','southeast')
xlim([1e-4,1e4])
xlabel('Scale (m)')
ylabel('Time (s)')
title('Acceleration Time')
% grid on


%take off velocity plot
figure
colororder([238 64 39;238 64 39;76 70 252;106 250 253; 0 115 252; 72 67 252]/255)
loglog(xim,vim,'LineWidth',2)
hold on
loglog(xim2,vim2,'--','LineWidth',2)
loglog(xO,vO,'LineWidth',2)

loglog(dMat,vto,':','LineWidth',2)

legend({'Eng., max power no linkage','Eng., max power opt. linkage','Bio., max power opt. linkage', 'Bio., G=0.01','Bio., G=0.1','Bio., G=1'},'Location','southeast')
xlim([1e-4,1e4])
xlabel('Scale (m)')
ylabel('Velocity (m/s)')
title('Takeoff Velocity')
% grid on

%take off KE
figure
colororder([238 64 39;238 64 39;76 70 252;106 250 253; 0 115 252; 72 67 252]/255)
loglog(xim,.5*vim.^2,'LineWidth',2)
hold on
loglog(xim2,.5*vim2.^2,'--','LineWidth',2)
loglog(xO,.5*vO.^2,'LineWidth',2)

loglog(dMat,.5*vto.^2,':','LineWidth',2)
loglog([1e-4 1e4],[u2 u2],'k--','LineWidth',2)

legend({'Eng., max power no linkage','Eng., max power opt. linkage','Bio., max power opt. linkage', 'Bio., G=0.01','Bio., G=0.1','Bio., G=1'},'Location','southeast')
xlim([1e-4,1e4])
xlabel('Scale (m)')
ylabel('Specific Energy(J/kg)')
title('COM Kinetic Energy')
% grid on

%plot of work
figure
colororder([238 64 39;238 64 39;76 70 252;106 250 253; 0 115 252; 72 67 252]/255)
loglog(xim,wim,'LineWidth',2)
hold on
loglog(xim2,wim2,'--','LineWidth',2)
hold on
loglog(xO,wO,'LineWidth',2)

loglog(dMat,w,':','LineWidth',2)
loglog([1e-4 1e4],[u2 u2],'k--','LineWidth',2)

legend({'Eng., max power no linkage','Eng., max power opt. linkage','Bio., max power opt. linkage', 'Bio., G=0.01','Bio., G=0.1','Bio., G=1'},'Location','southeast')
xlim([1e-4,1e4])
xlabel('Scale (m)')
ylabel('Specific Energy (J/Kg)')
title('Produced Energy')
% grid on

%Jump height
figure
colororder([238 64 39;238 64 39;76 70 252;106 250 253; 0 115 252; 72 67 252]/255)
loglog(xim,him,'LineWidth',2)
hold on
loglog(xim2,him2,'--','LineWidth',2)

loglog(xO,hO,'LineWidth',2)

loglog(dMat,h,':','LineWidth',2)



legend({'Eng., max power no linkage','Eng., max power opt. linkage','Bio., max power opt. linkage', 'Bio., G=0.01','Bio., G=0.1','Bio., G=1'},'Location','southeast')
xlim([1e-4,1e4])
xlabel('Scale (m)')
ylabel('Height(m)')
title('Estimated Jump Height')
% grid on



