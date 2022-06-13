% This numerical scheme is to analyse the system of equations (Eq 1a, 1b
% and 1c: as written down in the manuscript)

% This can be compiled and run using MATLAB. This generates a movie in .avi
% format on completion of the run time.

clc
clear all
close all

%advection terms are computed using upwinding schemes
%diffusion terms are computed using central difference schemes






%Parameters

Gamma=1;
B=4;
c=1;
D=1;
k1=1;
k2=0;k5=0;
k3=1;
zeta1=2.35;
zeta2=-2;
k4=3;
chi=1;
chip=0.01;
chipp=0.001;
chippp=0.001;
rho_s=(k1-(k2*k5))/((k1*k1)-(k2*k2));
phi_s=((k1*k5)-k2)/((k1*k1)-(k2*k2));
B1 = B - (c*c) - (2*c*zeta1*rho_s) - (2*c*zeta2*phi_s);


%dx=0.08;
%dx=0.2;
%dx=0.22;

%mesh details

dx=0.22;
l=6*pi;
x = (0:dx:l);
N=length(x);


%initial conditons

% u=u+(0.05*sin(x*pi)+0.05*sin(2*x*pi)+0.05*sin(3*x*pi)+0.05*sin(x*4*pi))/100;
% u=u+(0.05*cos(x*pi)+0.05*cos(x*2*pi)+0.05*cos(3*x*pi)+0.05*cos(x*4*pi))/100;
% 
% rho= rho_s*ones(1,N);
% rho=rho+(0.05*sin(x*pi)+0.05*sin(x*2*pi)+0.05*sin(3*x*pi)+0.05*sin(4*x*pi)) ;
% rho=rho+(0.05*cos(x*pi)+0.05*cos(x*2*pi)+0.05*cos(3*x*pi)+0.05*cos(x*4*pi));


sigma=0.4;

epsilon=zeros(1,N);
%u(5:end-4) = (normpdf(x(5:end-4),mean(x),sigma))/100;
noise=rand(1,N)-0.5;
%epsilon(5:end-4)=epsilon(5:end-4)+(noise(5:end-4))/200;
%epsilon=epsilon+(0.005*cos(8*x)) ;

%%noise=0;
%u=u-(0.005*sin(10*x*pi));

rho= rho_s*ones(1,N);
%rho(5:end-4) = rho_s+rho_s*(normpdf(x(5:end-4),mean(x),sigma))/100;
%rho(5:end-4)=rho(5:end-4)+(noise(5:end-4))/100;

phi= phi_s*zeros(1,N);
%phi=phi+(0.005*sin(10*x));
%phi=phi+(0.5*sin(8*x)) ;
phi(5:end-4)=phi(5:end-4)+(noise(5:end-4))/6;
%phi=phi+0.5*(1+tanh(x-40)-tanh(x-5));
%phi=phi+0.1*(1+tanh(x-30)-tanh(x-25));


%phi=phi+ 0.5 - exp(-(x-25).*(x-25)/100);

%phi(5:end-4) = phi_s+0.1*(normpdf(x(5:end-4),mean(x),sigma))/100;
%phi(5:end-4)=0.01*sin(8*x((5:end-4)));


epsilon0=epsilon;
rho0=rho;
phi0=phi;

%strain=du_dx_central(u,dx);
%strain0=strain;
dt=min(0.1*dx^2/B, 0.1*dx^2/D)/200 
%dt=0.000001;

% eps=10^(-6);
T=0;


folder_name=sprintf('coexistence');
folder=fullfile(folder_name);
 mkdir(folder);
video_name=sprintf('trav_stand');
video_file=fullfile(folder,video_name);
vd = VideoWriter(video_file);  
 
%plotting initial conditons 
 
 open(vd);
 xp0=10;
yp0=10;
width=1000;
height=750;
set(gcf,'position',[xp0,yp0,width,height])
 plot_size=10;
 subplot(3,1,1); plot(x,phi,'g','LineWidth',6);  title({[' T= ',num2str(T)]});set(gca,'FontSize',24,'linewidth',6);ylabel('\phi \rightarrow','FontSize',24);
  xlim([0 l])
 subplot(3,1,2); plot(x,epsilon,'g','LineWidth',6); set(gca,'FontSize',24,'linewidth',6);ylabel('\epsilon \rightarrow','FontSize',24);
  xlim([0 l])
subplot(3,1,3); plot(x,rho,'g','LineWidth',6); set(gca,'FontSize',24,'linewidth',6);ylabel('\rho \rightarrow','FontSize',24);
 xlim([0 l])
 xlabel('x \rightarrow','FontSize',24)
 
 %Main loop
 
 Fr = getframe(gcf);
 writeVideo(vd,Fr);j=1; aj=100;ajj=1;
 tic
for i=1:1:42000

% evaluation equation 1

% term 1:= B u_xx central difference scheme
depsilon=du_dx_higher(epsilon,dx);ddepsilon=du_dx_higher(depsilon,dx);
drho=du_dx_higher(rho,dx);ddrho=du_dx_higher(drho,dx);
dphi=du_dx_higher(phi,dx);ddphi=du_dx_higher(dphi,dx);



t1= chi * 2.* ((zeta1.*drho)+(zeta2.*dphi));

t2= ((B-(c*c)-(2*chip*c).* ((zeta1.*rho)+(zeta2.*phi))).* depsilon) + ((-(2*chip*c).*((zeta1.*drho)+(zeta2.*dphi))).* epsilon);

t3= (((chipp*c*c).* (zeta1.*drho)+(zeta2.*dphi)) .* epsilon .* epsilon) + (((chipp*c*c).* (zeta1.*rho)+(zeta2.*phi)).* 2 .* epsilon .* depsilon);

t4= ((-chippp*c*c*c*(1/3)).* ((zeta1.*drho)+(zeta2.*dphi)) .* epsilon .* epsilon .*epsilon) + (((-chippp*c*c*c*(1/3)).* ((zeta1.*rho)+(zeta2.*phi))).* 3 .* epsilon .* epsilon .* depsilon);


v=t1+t2+t3+t4;


%tt1= chi * 2.* ((zeta1.*ddrho)+(zeta2.*ddphi));

%tt2= ((B-(c*c)-(2*chip*c).* ((zeta1.*rho)+(zeta2.*phi))).* ddepsilon)+ ((-(2*chip*c).*((zeta1.*drho)+(zeta2.*dphi))).* depsilon) +  ((-(2*chip*c).*((zeta1.*drho)+(zeta2.*dphi))).* depsilon) +  ((-(2*chip*c).*((zeta1.*ddrho)+(zeta2.*ddphi))).* epsilon);

%tt31= (  ((chipp*c*c).* (zeta1.*drho)+(zeta2.*dphi)).* 2 .* epsilon .* depsilon  ) + (  ((chipp*c*c).* (zeta1.*rho)+(zeta2.*phi)).* 2 .* depsilon .* depsilon  ) + (  ((chipp*c*c).* (zeta1.*rho)+(zeta2.*phi)).* 2 .* epsilon .* ddepsilon  );

%tt32= (  ((chipp*c*c).* (zeta1.*drho)+(zeta2.*dphi)).* 2 .* epsilon .* depsilon  ) + (  ((chipp*c*c).* (zeta1.*ddrho)+(zeta2.*ddphi)).* epsilon .* epsilon  );

%tt41= (   (-chippp*c*c*c*(1/3)).* ((zeta1.*drho)+(zeta2.*dphi)) .* 3 .* epsilon .* epsilon .* depsilon  ) +  (   (-chippp*c*c*c*(1/3)).* ((zeta1.*rho)+(zeta2.*phi)) .* 6 .* epsilon .* depsilon .* depsilon  ) + (   (-chippp*c*c*c*(1/3)).* ((zeta1.*rho)+(zeta2.*phi)) .* 3 .* epsilon .* epsilon .* ddepsilon  );

%tt42= (   (-chippp*c*c*c*(1/3)).* ((zeta1.*ddrho)+(zeta2.*ddphi)) .* epsilon .* epsilon .* epsilon  ) + (   (-chippp*c*c*c*(1/3)).* ((zeta1.*drho)+(zeta2.*dphi)) .* 3 .* epsilon .* epsilon .* depsilon  );


vhigh= du_dx_higher (v,dx);

%vhigh= tt1+tt2+tt31+tt32+tt41+tt42;

 
epsilon_new= epsilon + vhigh*dt/Gamma;
 
% evaluation equation 2


 T1= D.*ddrho;
 T2= 1-(c.* epsilon);
 T3= -(k1+(k3.* epsilon)).*rho;
 T4= -(k2+(k4.* epsilon)).*phi;
 T5=du_dx_upwind_flux(rho,v,dx,dt);
 
 rho_new= rho+dt.*(T1+T2+T3+T4-T5); 
 
 
 %evaluation equation 3
 
 TT1= D.*ddphi;
 TT2= k5.*(1-(c.* epsilon));
 TT3= -(k1+(k3.* epsilon)).*phi;
 TT4= -(k2+(k4.* epsilon)).*rho;
 TT5=du_dx_upwind_flux(phi,v,dx,dt);
 phi_new= phi+dt.*(TT1+TT2+TT3+TT4-TT5); 

 
 T=T+dt;
 epsilon=epsilon_new;
 rho=rho_new;
 phi=phi_new;
  if rem(i,ajj)==0 
         plot_size=10;
        subplot(3,1,1); plot(x,phi,'b','LineWidth',6);title({['T= ',num2str(T)]}); ylabel('\phi','fontsize',24);set(gca,'FontSize',24,'linewidth',6);
         ylim([-0.5 0.9]);xlim([0 l]);
        subplot(3,1,2); plot(x,epsilon,'r','LineWidth',6);ylabel('\epsilon \rightarrow','fontsize',24);set(gca,'FontSize',24,'linewidth',6);
         xlim([0 l]);ylim([-2 2]);
        subplot(3,1,3); plot(x,rho,'g','LineWidth',6); set(gca,'FontSize',24,'linewidth',6);ylabel('\rho \rightarrow','FontSize',24);
        ylim([0.2 2.2]);xlim([0 l]);
        
        xlabel('x \rightarrow','fontsize',24)
         if rem(i,aj)==0&& j<=100
             dt=dt*2;j=j+1; aj=aj*2;ajj=ajj*j*2; ajj=min(ajj,50);
         end
          
            
            Fr = getframe(gcf);
            writeVideo(vd,Fr);
   end

end
toc
close(vd);

%functions

function res=du_dx_central(y,dx)

yf=[y(2:end),y(2)];                                  % 1step forward node in y
yb=[y(end-1),y(1:end-1)];                            % 1step backward node in y

% xf=[x(2:end),x(end)+x(2)-x(1)];                    % 1step forward node in x
% xb=[x(1)-x(end)+x(end-1),x(1:end-1)];
res=(yf-yb)./(2*dx);

end

function res=d2u_dx2_central(y,dx)
yf=[y(2:end),y(2)];                                 % 1step forward node in y
yb=[y(end-1),y(1:end-1)];                           % 1step backward node in y

%xf=[x(2:end),x(end)+x(2)-x(1)];                    % 1step forward node in x
%xb=[x(1)-x(end)+x(end-1),x(1:end-1)];

res=(yf-2*y+yb)/(dx^2);

end

function res=du_dx_upwind(y,dx,v)
yf=[y(2:end),y(2)];                                     % 1step forward node in y
yb=[y(end-1),y(1:end-1)];                               % 1step backward node in y
yff=[yf(2:end),yf(2)];
ybb=[yb(end-1),yb(1:end-1)];

% xf=[x(2:end),x(end)+x(2)-x(1)];                    % 1step forward node in x
% xb=[x(1)-x(end)+x(end-1),x(1:end-1)];

res1=(3*y-4*yb+ybb)./(2*dx);
res2=(-3*y+4*yf-yff)./(2*dx);
if v >0
    res =res1;
else
    res=res2;
end



end


function res=du_dx_higher(y,dx)

yf=[y(2:end),y(2)];                                     % 1step forward node in y
yb=[y(end-1),y(1:end-1)];                               % 1step backward node in y
yff=[yf(2:end),yf(2)];   
yfff=[yff(2:end),yff(2)];   
ybb=[yb(end-1),yb(1:end-1)];     
ybbb=[ybb(end-1),ybb(1:end-1)];   


% xf=[x(2:end),x(end)+x(2)-x(1)];                    % 1step forward node in x
% xb=[x(1)-x(end)+x(end-1),x(1:end-1)];
res=(yfff-9*yff+45*yf-45*yb+9*ybb-ybbb)/(60*dx);

end



function res=du_dx_upwind_flux(rho,v,dx,dt)
rhof=[rho(2:end),rho(2)];                            % 1step forward node in y
rhob=[rho(end-1),rho(1:end-1)];                      % 1step backward node in y
rhoff=[rhof(2:end),rhof(2)];
rhobb=[rhob(end-1),rhob(1:end-1)];
vf=[v(2:end),v(2)];  
% xf=[x(2:end),x(end)+x(2)-x(1)];                    % 1step forward node in x
% xb=[x(1)-x(end)+x(end-1),x(1:end-1)];
v_int=0.5*(vf+v);
v_intb=[v_int(end-1),v_int(1:end-1)]; 

res=zeros(1,length(rho));
for i=1:1:length(rho)
    if v_int(i) >= 0
        theta=1;
        r=(rho(i)-rhob(i)+10^-6)./(rhof(i)-rho(i)+10^-6);
    else
        theta=-1;
        r=(rhoff(i)-rhof(i)+10^-6)./(rhof(i)-rho(i)+10^-6);
    end
    
    limiter=(r+abs(r))/(1+abs(r));
    
	j1 = 0.5*v_int(i)*((1+theta)*rho(i)+(1-theta)*rhof(i));
	j2 = 0.5*abs(v_int(i))*(1-abs(v_int(i)*dt/dx))*limiter*(rhof(i)-rho(i));
	j_plus = j1 + j2;
    if v_intb(i) >=0

        theta = 1;
        r = (rhob(i)-rhobb(i)+10^-6)/(rho(i)-rhob(i)+10^-6);

    else

        theta = -1;
        r= (rhof(i)-rho(i)+10^-6)/(rho(i)-rhob(i)+10^-6);	
    end
        limiter=(r+abs(r))/(1+abs(r));
	
	j1 = 0.5*v_intb(i)*((1+theta)*rhob(i)+(1-theta)*rho(i));
	j2 = 0.5*abs(v_intb(i))*(1-abs(v_intb(i)*dt/dx))*limiter*(rho(i)-rhob(i));
	j_mins = j1 + j2;

	res(i) = (j_plus - j_mins)/dx;
end
    
end




















