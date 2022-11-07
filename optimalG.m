%this function finds the optimal gear ratio for a given scale
u=200; %[J/kg] specific energy bio
p=200;%200*.14; %[W/kg] specific power bio
ulb=1130;%[J/kg] energy transfer ability of biological linkages

u2=u*ulb/(ulb+u);
p2=p*ulb/(ulb+u);

g=9.81;
n=100; %number of indices 
x=logspace(-4,4,n);

syms G d
vto1=4*d*p2*(G*u2-g*d)*(1+lambertw(-exp(-(16*p2^2*d^2*(d*g-G*u2)-G^4*u2^4)/(16*d^2*p2^2*(d*g-G*u2)))))./(G^2*u2^2); %velocity function for G<=1
vto2=4*d*p2*(G*u2-g*d)*(1+lambertw(-exp(-(16*G*p2^2*d^2*(d*g-G*u2)-G^4*u2^4)/(16*G*d^2*p2^2*(d*g-G*u2)))))./(G^2*u2^2); %velocity function for G>1;
dvDG1=matlabFunction(diff(vto1,G));
dvDG2=matlabFunction(diff(vto2,G));
Gopt1=zeros(1,n);
Gopt2=zeros(1,n);

for i=1:n
    if x(i)<=u2/g
        G_o1=fzero(@(G) dvDG1(G,x(i)),.8546*x(i)^.6913); %the guess is emprically determined
%         G_o1=fzero(@(G) dvDG1(G,x(i)),.5*x(i)); 
        Gopt1(i)=G_o1;
    end
    G_o2=fzero(@(G) dvDG2(G,x(i)),x(i)*2.5*g/u2);
    Gopt2(i)=G_o2;
end

%this makes the plot for the optimal gear ratio
Gopt1(Gopt1>1)=1; %for bio jumper. Comment out for non bio jumper
mask1=Gopt1==0 | Gopt1>1;
% mask2=mask1 & Gopt2>1;
% Gopt=[Gopt1(~mask1) Gopt2(mask2)];

% loglog([x(~mask1) x(mask2)],Gopt,'LineWidth',2)
loglog(x(~mask1), Gopt1(~mask1),'Color',[85 76 252]/255,'LineWidth',2)
xlabel('Scale (m)')
ylabel('G_{opt}')
title('Optimal Fixed Gear Ratio (G\leq1)')
% grid on
xlim([1e-4 1e4])
ylim([1e-3 1e1])