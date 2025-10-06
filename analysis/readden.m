clear;
close all;
dt=0.3;
omegaA=1.1447515591120262E-002;%2.8040572020799567E-002;%2.9933892787182888E-002;
phirphi=load('testphi_r_phi');
%ttt=importdata('testAtPhitrhotjt');
R=load('R.dat');
Z=load('Z.dat');
%xforw=load('testxforw');
%zforw=load('testzforw');
BR=load('BR.dat');
BZ=load('Bz.dat');
Bt=load('Bt.dat');
%testB=load('testB');
%testBR=load('testBR');
%testBZ=load('testBZ');
%testBt=load('testnBt');
mask2=load('mask2');
mask3=load('mask3');
mask4=load('mask4');
%testphi=load('testphi');
%testphis=load('testphis');
masktest=load('mask');
%psitest=load('psi_test');
%nitest=load('testni0');
%apar0=load('testapar0');
ne1=load('testden');
ne2=load('testden2');
dden=load('testdiffden');
upar=load('testupar');
%apar1=load('testapars');
%apar2=load('testapar');
%jpar1=load('testjpars');
%jpar2=load('testjpar');
%jpar0=load('testj0');
%ne0zeta=load('testne0_zeta');
%j0zeta=load('testjpar0_zeta');
c2_over_vA2=load('test_1_over_vA2');
%mask=load('debug.dat');

p=0:2*pi/33:2*pi*(1-1/33);
%phi=zeros(434,450);
for i=1:434
    for j=1:450

        ii=floor(((i-1)*450+j-1)/3)+1;
        jj=(i-1)*450+j-3*(floor(((i-1)*450+j-1)/3));

        %xf(i,j)=xforw(ii,jj);
        %zf(i,j)=zforw(ii,jj);
%        phi1(i,j)=testphis(ii,jj);
%        phi2(i,j)=testphi(ii,jj);
%        a0(i,j)=apar0(ii,jj);
        m1(i,j)=masktest(ii,jj);
%        psi(i,j)=psitest(ii,jj);
        dn(i,j)=dden(ii,jj);
        n1(i,j)=ne1(ii,jj);
        n2(i,j)=ne2(ii,jj);
        up(i,j)=upar(ii,jj);
         %ni(i,j)=nitest(ii,jj);
%        a1(i,j)=apar1(ii,jj);
%        a2(i,j)=apar2(ii,jj);
%        j0(i,j)=jpar0(ii,jj);
%        j1(i,j)=jpar1(ii,jj);
%        j2(i,j)=jpar2(ii,jj);
        m2(i,j)=mask2(ii,jj);
        m3(i,j)=mask3(ii,jj);
        m4(i,j)=mask4(ii,jj);
        vA2(i,j)=c2_over_vA2(ii,jj);
        %tBR(i,j)=testBR(ii,jj);
        %tBZ(i,j)=testBZ(ii,jj);
        %tBt(i,j)=testBt(ii,jj);
        %tB(i,j)=testB(ii,jj);
    end
end

%for i=1:16
%    for j=1:450
%        ii=floor(((i-1)*450+j-1)/3)+1;
%        jj=(i-1)*450+j-3*(floor(((i-1)*450+j-1)/3));
%        prphi(i,j)=phirphi(ii,jj);
%    end
%end

figure
[c,h]=contourf(R,Z,n1,100);
axis equal
set(h,'LineColor','none')
title('n_{i0}','Fontsize',24)
%Fontsize(gcf,24,"points")


figure
[c,h]=contourf(R,Z,n2,100);
axis equal
set(h,'LineColor','none')
title('n_{iend}','Fontsize',24)
%fontsize(gcf,24,"points")



%figure
%[c,h]=contourf(R,Z,ni,100);
%axis equal
%set(h,'LineColor','none')
%title('ni0','Fontsize',24)

figure
%[c,h]=contourf(R,Z,dn,100);
[c,h]=contourf(R,Z,n2-n1,100);
axis equal
set(h,'LineColor','none')
title('dn_i','Fontsize',24)

figure
[c,h]=contourf(R,Z,up,100);
axis equal
set(h,'LineColor','none')
title('n_iu_{||}','Fontsize',24)

figure
[c,h]=contourf(R,Z,vA2,100);
set(h,'LineColor','none')
xlabel('R(m)','Fontsize',24)
title('1/v_A^2','Fontsize',24)
axis equal






%valuea=(m1-1).*a2;
%valuej=(m2-2).*j2;
%valuen=(m3-3).*n2;
%valuephi=(m1-1).*phi2;

%t=floor(length(ttt)/2);
%tt=t*2;
%omega=[0:t-1]/t/omegaA;
%figure
%plot(ttt(1:2:tt,2))
%ppp=fft(ttt(1:2:tt,2));
%pp=abs(ppp);
%for i=2:length(pp)
%    pp(i)=2*pp(i)/t;
%end
%figure
%plot(omega,pp)
