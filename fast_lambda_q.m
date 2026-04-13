close all;
clear all;
p = load("./test_hit_div");
hit_t=load("test_hit_div_t");
Z=load("Z.dat");
ini_p=load("initial_posi");
psi_a = 0;
number_of_beams = 150;
dt=1500e-10;
mi=2*1.6726e-27;
e=1.6e-19;
T=200;
vth=sqrt(T*e/mi);
%%%%%%%%%%%%%% two points define the divertor plate%%%%%%%%%%%%%%%%%%%%%%%5
r1 = 1.2;  %1.23
z1 = 0.25773+Z(1); %-1.4
r2 = 1.5;    %1.48
z2 = 0.25773+Z(1);  %-1.2
m = (z2-z1)/(r2-r1);
b = z1-m*r1;
R_div = linspace(r1,r2,400);
Z_div = linspace(z1,z2,400);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = load("R.dat");
Z = load("Z.dat");
psi = load("psi_p.dat");
%energy= load("./energeez");
%u2 = load("initialv_par");
phi=load("phisolps.dat");
psip= psi';
phip=phi';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phip(:,:)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hit_num(1:max(hit_t))=0;%zeros(length(hit_t));
for i=1:length(hit_t)
    if (hit_t(i)~=0)
    hit_num(hit_t(i))=hit_num(hit_t(i))+1;
    end
end
for i=2:length(hit_num)
    hit_num(i)=hit_num(i-1)+hit_num(i);
end
plot(hit_num);
figure;

[Mid_z,zmid_id] = min(abs(Z)); 

Rmid_id = floor(length(R)/2);
sgin_of_psi = (psip(length(R),zmid_id)-psip(Rmid_id,zmid_id))/abs(psip(length(R),zmid_id)-psip(Rmid_id,zmid_id));

dR=R(2)-R(1);
dZ=Z(2)-Z(1);
[dpdR,dpdZ] = gradient(psip,dZ,dR);

[xg,zg] = ndgrid(R,Z);
%test=[1,2,3;5,6,7;9,10,11];
%[d1,d2] = gradient(test,1,1); 
gradpsiabs = sqrt(dpdR.^2+dpdZ.^2);

[M,Z_x_id] = min(min(gradpsiabs(2:end-2,5:zmid_id-5)));
[M,R_x_id] = min(min(gradpsiabs(2:end-2,5:zmid_id-5)'));
R_x_id = R_x_id+1;
Z_x_id = Z_x_id+4;
R_x = R(R_x_id);
Z_x = Z(Z_x_id);


%number_of_particles = i-1;
%t_max = length(p)/number_of_particles;
f = griddedInterpolant(xg,zg,psip);
g = griddedInterpolant(xg,zg,gradpsiabs);
phis = griddedInterpolant(xg,zg,phip);

%g(0,2.1197)

for i = floor(length(R)/2):length(R)
    if (psip(i-1,zmid_id)-psi_a)*(psip(i,zmid_id)-psi_a)<=0
        R_m_a=(R(i-1)*abs(psip(i,zmid_id)-psi_a)+R(i)*abs(psip(i-1,zmid_id)-psi_a))/abs(psip(i,zmid_id)-psip(i-1,zmid_id));
        gradpsim = (gradpsiabs(i-1,zmid_id)*abs(psip(i,zmid_id)-psi_a)+gradpsiabs(i,zmid_id)*abs(psip(i-1,zmid_id)-psi_a))/abs(psip(i,zmid_id)-psip(i-1,zmid_id));
        %gradpsim = g(R_m_a,Mid_z);
        idx_R_a = i;
        break;
    end
end
R_M_a = 2.1177;
[M,idx_R_div]= min(abs(f(R_div,Z_div)));
R_inter = R_div(idx_R_div);
Z_inter = Z_div(idx_R_div);
%f(R_div(idx_R_div),Z_div(idx_R_div))
gradpsidiv= g(R_inter,Z_inter);

Ep=0;
dpsi=0;
Rim=0;
R0_loss=0;
Z0_loss=0;
udiv=0;
Rim_w_e=0;
number_of_particles=0
for i = 1:length(p)
    if p(i,2) ~= 0%(i-1,2)
        number_of_particles=number_of_particles+1;
        dpsi(number_of_particles) = p(i,1)-psi_a;
        Rim(number_of_particles)   = sgin_of_psi*dpsi(number_of_particles)/gradpsim;
        Rim_w_e(number_of_particles) = Rim(number_of_particles)*p(i,3);
        Ep(number_of_particles) = p(i,3);%(energy(i)+e*(phis(xtrg(1),ztrg(1))-phis(xtrg(k),ztrg(k))))/(R_m_a+Rim(np));%energy(i)/(R_m_a+Rim(np));
        udiv(number_of_particles) = abs(p(i,2));

    end
end




p_inner_divertor=0;
plot([r2,r1],[z2,z1]);
hold on
contour(R,Z,psi);
plot(ini_p(:,1),ini_p(:,2),LineStyle="none",Marker=".");
axis equal;


%figure
%contour(R,Z,psi,50);
%hold on;
%scatter(R0_loss,Z0_loss);
%axis equal;
%title("Initial positions of the particles hit the divertor plate​")
%scatter(p(1:number_of_particles,3),p(1:number_of_particles,4))
%set(gca,'FontSize',24)

%lost_t=sort(lost_t);
%lost_p=0;



Ridiv = dpsi/gradpsidiv;
%Rim   = dpsi/gradpsim;

lambdadiv = sum(abs(Ridiv))/length(dpsi)
lambdamid = sum(abs(Rim))/length(dpsi)
lambdamid_w_e=sum(abs(Rim_w_e))/sum(Ep)


R_midgrid = linspace(min(Rim)-0.0001,max(Rim)+0.0001,number_of_beams+1);
dRmid = (R_midgrid(end)-R_midgrid(1))/number_of_beams;
R_beam=0;
for i = 2:length(R_midgrid)
    R_beam(i-1) = 0.5*(R_midgrid(i-1)+ R_midgrid(i));
end
q_mid=zeros(1,length(R_beam));
%q_mid=[];
%min_i=min(floor(((Rim)-min(R_midgrid))/dRmid));
%max_i=max(floor(((Rim)-min(R_midgrid))/dRmid));
for i = 1:length(dpsi)
    pp=floor((Rim(i)-min(R_midgrid))/dRmid)+1;
    q_mid(pp)=q_mid(pp)+Ep(i)/dRmid;
end

%x = R_beam*1000;
%y = q_mid;
%lambdamid=lambdamid*1000;
%testfittype = fittype(@(a,b,c,x) a*exp(c^2/4/lambdamid^2-x/lambdamid)+b);%.*erfc(-x/b)+c);
%eich = fittype(@(q0,qbg,S,s0,x) q0/2*exp((S/(2*lambdamid))^2-(x-s0)/lambdamid).*erfc(S/(2*lambdamid)-(x-s0)/S)+qbg)
%coeffnames(eich)
%[fitted,gov]=fit(x',y',eich)
q0_guess =max(q_mid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lambdamid=lambdamid*1000;
%R_beam=R_beam*1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i=1:10
%    R_beam(number_of_beams+i)=R_beam(number_of_beams)+i*(R_midgrid(2)-R_midgrid(1));
%    q_mid(number_of_beams+i)=0;
%end


start_guess = [2.4*q0_guess,0,0.9*lambdamid,-0.000];
%s0f=0.001;
%[fitted,gov]=fit (R_beam',q_mid',eich,'Start',start_guess)

%[testfitted,gottest] = fit (R_beam',q_mid',testfittype);
%coeffvals = coeffvalues(fitted)

xx=R_beam(1):0.0001:R_beam(end);
%y=fitted(x);
%x=x';
%c=coeffvalues(fitted);
c=start_guess;
q0f=c(1);
qbgf=c(2);
Sf=c(3);
s0f=c(4);
yy=q0f/2*exp((abs(Sf)/(2*lambdamid))^2-(xx-s0f)/lambdamid).*erfc(abs(Sf)/(2*lambdamid)-(xx-s0f)/abs(Sf))+qbgf;


figure
plot(R_beam,q_mid);
hold on
plot(xx,yy)
xlabel('s-s_0 (m)');
ylabel('q (a.u.)')
set(gca,'FontSize',24)


figure
hist(udiv/vth,100);
xlabel('v_{//div}/v_{i th}');
ylabel('Number of particles hit the outer divertor plate')
set(gca,'FontSize',24)
udiv_avg=sum(udiv)/length(udiv)/vth
lambdamid
lambdamid_w_e
figure


plot([r2,r1],[z2,z1]);
hold on
contour(R,Z,psi);
axis equal;

lambda_s=0.0005;
lambda_d=0.2;
nuj0=3.8e-2;%2e5;
for i=1:length(R_beam)
    if R_beam(i)<0
        nujp(i)=nuj0*exp(-(R_beam(i))/lambda_s)-1;
    else
        nujp(i)=0;
    end
end


nuj=nujp.*exp(-0.2/lambda_d);
factor=1;%./(1+nuj*0.2/vth);
lambda_avg=sum(q_mid.*R_beam.*factor)/sum(q_mid)

%plot(factor);
figure
plot(R_beam,q_mid)%.*factor);
hold on
plot(xx,yy)
xlabel('s-s_0 (m)');
ylabel('q (a.u.)')
set(gca,'FontSize',24)