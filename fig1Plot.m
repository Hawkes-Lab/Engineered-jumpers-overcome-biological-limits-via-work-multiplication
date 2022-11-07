mm=.005;%[kg] mass of motor
mp=.013;%mass of dead weight
ms=.012;%total mass of spring
m1=mm+mp;
m2=ms/2;
g=9.81;
l0=.3;%legnth of the legs
l=l0/2;%[m] length of legs in model
I=m2*(.05^2+l^2)/12;
thf=pi/2; %Rest angle of legs
th0=pi/6;%inital leg angle
CdA=.018^2;%6.4e-4; %m^2 combined coefficient of drag
rho=1.225; %kg/m^3 
pm=2000; %[w/kg] estimated speific power of motor
p=pm*mm;%[W] power of the motor
us=2555; %J/kg effective specific energy of legs
k=ms*us/(2*(thf-th0)^2); %Effective spring constant Nm/rad


%%
%Spring loading process
tspan0=[0 5];
opts0 = odeset('Events', @(t,th) springStop(t,th,th0),'RelTol',1e-3,'AbsTol',1e-4);
[t0,th_0,te0,the0,ie0]=ode45(@(t,th) Spring(t,th,m1,m2,I,k,g,l,thf,p),tspan0,[thf;-0.01],opts0); %inital velocity to prevent blow up
Us0=2*k*(thf-th_0(:,1)).^2;

%%
%Takeoff process
tspan=[0 0.05];
opts = odeset('Events', @(t,th) liftOff(t,th,m1,m2,I,k,g,l,thf),'RelTol',1e-6,'AbsTol',1e-6);

[t,th,te,the,ie]=ode45(@(t,th) odeFun1(t,th,m1,m2,I,k,g,l,thf),tspan,[th0;0],opts);

%Calculate energy and momentum
y1=2*l*sin(th(:,1)); vy1=2*l*th(:,2).*cos(th(:,1));
x2=-l/2*cos(th(:,1)); vx2=l/2*th(:,2).*sin(th(:,1));
y2=3/2*l*sin(th(:,1));vy2=3/2*l*th(:,2).*cos(th(:,1));
y3=1/2*l*sin(th(:,1)); vy3=1/2*l*th(:,2).*cos(th(:,1));

Us=2*k*(thf-th(:,1)).^2;
Ug=g*(m1+m2)*y1;%potential gravity energy
Tx=m2*vx2.^2;
Ty=1/2*m1*vy1.^2+1/2*m2*vy2.^2+1/2*m2*vy3.^2;
Tw=I*th(:,2).^2;
Py=m1*vy1+m2*vy2+m2*vy3;

Wp=Tx+Ty+Tw+Ug-Ug(1); %produced energy
Wy=Ty+Ug-Ug(1);%utilized energy
WyCOM=1/2*(m1+2*m2)*(Py/(m1+2*m2)).^2+Ug-Ug(1);%utilized energy from COM perspective

vto=Py(end)/(m1+2*m2);

tspan2=[0 1];
opts2 = odeset('Events', @(t,y) Apex(t,y),'RelTol',1e-6,'AbsTol',1e-6);


%%
%flight process 
[t2,y,te2,ye,ie2]=ode45(@(t,y) Flight(t,y,m1,m2,g,CdA,rho),tspan2,[0;vto],opts2);

KE2=1/2*(m1+2*m2)*y(:,2).^2;
Ug2=(m1+2*m2)*g*(y(:,1)-y(1,1))+Ug(end)-Ug(1);
Wy2=KE2+Ug2;%prduced energy 2


%%
%plots 
c=1000;%time stretch factor for launch process
c2=10;%time strech factor for the flight

plot(t0-t0(end),Us0,'--','LineWidth',2,'Color',[242 133 54]/255)
hold on
plot(t*c,Us,'--','LineWidth',2,'Color',[242 133 54]/255)
plot(t*c,Wp,'LineWidth',2,'Color',[80 142 252]/255);
plot([t; t(end)]*c,[Wy; Wy2(1)],'LineWidth',2,'Color',[237 59 46]/255)
plot([t(end)*c;t(end)*c+t2*c2],[1/2*(Wy(end)+Wy2(1));Wy2],'LineWidth',2,'Color',[237 59 46]/255)
xlabel('Time')
ylabel('Energy')
set(gca,'YTick',[])
set(gca,'XTick',[])

%for COM based KE
figure
plot(t0-t0(end),Us0,'--','LineWidth',2,'Color',[242 133 54]/255)
hold on
plot(t*c,Us,'--','LineWidth',2,'Color',[242 133 54]/255)
plot(t*c,Wp,'LineWidth',2,'Color',[80 142 252]/255);
plot([[t; t(end)]*c;t(end)*c+t2*c2],[WyCOM; Wy2(1);Wy2],'LineWidth',2,'Color',[237 59 46]/255)

xlabel('Time')
ylabel('Energy')
set(gca,'YTick',[])
set(gca,'XTick',[])

% xlim([0 t(end)*2]) % for the detail about the jump section
%%
%diffeq that describes launch
function dthdt=odeFun1(t,th,m1,m2,I,k,g,l,thf)

    dthdt=[th(2);
            (8*k*(thf-th(1))-4*(m1+m2)*g*l*cos(th(1))+2*l^2*(2*m1+m2)*th(2)^2*sin(2*th(1)))/...
            (4*I+l^2*(4*m1+3*m2)+2*l^2*(2*m1+m2)*cos(2*th(1)))];
end
function [position,isterminal,direction] = liftOff(t,th,m1,m2,I,k,g,l,thf)

  position = (m1+2*m2)*g-2*l*(m1+m2)*(th(2)^2*sin(th(1))-cos(th(1))*(8*k*(thf-th(1))-4*(m1+m2)*g*l*cos(th(1))+2*l^2*(2*m1+m2)*th(2)^2*sin(2*th(1)))/...
            (4*I+l^2*(4*m1+3*m2)+2*l^2*(2*m1+m2)*cos(2*th(1)))); % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end

%diffeq that describes flight
function dydt=Flight(t,y,m1,m2,g,CdA,rho)
    dydt=[y(2);
        -1/(m1+2*m2)*(1/2*CdA*rho*y(2)^2+(m1+2*m2)*g)];
end
function [position,isterminal,direction] = Apex(t,y)

  position = y(2);
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end

%diffeq that describes spring loading
function dthdt=Spring(t,th,m1,m2,I,k,g,l,thf,p)

    dthdt=[th(2);
            (2*p+(8*k*(thf-th(1))-4*(m1+m2)*g*l*cos(th(1))+2*l^2*(2*m1+m2)*th(2)^2*sin(2*th(1)))*th(2))/...
            (th(2)*(4*I+l^2*(4*m1+3*m2)+2*l^2*(2*m1+m2)*cos(2*th(1))))];
end
function [position,isterminal,direction] = springStop(t,th,th0)

  position = th(1)-th0;
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end