clear all
clc
close all
lineWidth = 1.5;
load('SaTrat')  % Acelerogramas Tratados y recortados
load('EspResp') % Espectros de Respuesta
n=length(SaTrat);
 
%% ----------------------- ESPECTRO DE DISEÑO ------------------------------
Fa=1.2;
Fd=1.19; 
Fs=1.28; 
r=1;             % Todo tipo de suelo excepto tipo E 
eta=2.48;        % Provincias de la sierra 
dT=0.01;
Tmax= 3;
zona_sismica=0.4*981; %m/s2
[SaNEC,T]=EspectroNEC2020A(Fa,Fd,Fs,r,eta,dT,Tmax,zona_sismica);
 
%% ---------------------- Graficos Espectros sin SS ------------------------
% Espectros componente N
figure(1)
plot(T,SaNEC,'LineWidth',lineWidth)
hold on
for i=1:n
    plot(T,EspResp{1,i})
end
grid on
grid minor
title('Espectros componente N-S')
legend('NEC15','s368','s540','s741','s752','s753','s953','s960','s1004','s1080','s1048','s4886')
xlabel('Periodo [s]')
ylabel('Psa N-S [cm/s^2]')
 
figure(2)
plot(T,SaNEC,'LineWidth',lineWidth)
hold on
for i=1:n
    plot(T,EspResp{2,i})
end
grid on
grid minor
title('Espectros componente E-O')
legend('NEC15','s368','s540','s741','s752','s753','s953','s960','s1004','s1080','s1048','s4886')
xlabel('Periodo [s]')
ylabel('Psa E-O [cm/s^2]')
    
%% -------------------------------------- Escalamiento Sismico ------------------------------------
% ------------------------------------------ Espectro SRSS ---------------------------------------
figure(3)
plot(T,SaNEC,'LineWidth',lineWidth)
hold on
SRSS = {0};
for i=1:n
    SRSS{i} = ((EspResp{1,i}).^2+(EspResp{2,i}).^2).^0.5;
     plot(T,SRSS{i})
end
grid on
grid minor
title('Espectros SRSS')
legend('NEC15','s368','s540','s741','s752','s753','s953','s960','s1004','s1080','s1048','s4886')
xlabel('Periodo [s]')
ylabel('Aceleración [cm/s^2]')
 
%% ------------------------------------------- Periodos ETABS -------------------------------------------------------
%Período en x
Tx=0.78; %sec
%Período en y
Ty=0.78; %sec
% Promedio
Tavg=(Tx+Ty)/2;
%% Rango de periodos
% 0.2T a 1.5T
Tinf=0.2*Tavg;
Tsup=1.5*Tavg;
 
% Truncar periodos
Tavg=fix(Tavg*100)/100;
display(Tavg)
Tinf=fix(Tinf*100)/100;
display(Tinf)
Tsup=fix(Tsup*100)/100;
display(Tsup)
 
% Valor de NEC en Tavg
index=find(T==Tavg);
index_inf=find(T==Tinf);
display(index_inf)
index_sup=find(T==Tsup);
display(index_sup)
NECavg=SaNEC(index);
disp('El valor de la aceleración NEC en Tavg es')
display(NECavg)
%% Factor S1
S1=zeros(n,1);
for i=1:n
    S1(i)=NECavg/SRSS{i}(index);
end
disp('S1=')
disp(S1)
%% SRSS*S1
figure(4)
plot(T,SaNEC,'LineWidth',lineWidth)
hold on
for i=1:n
    plot(T,SRSS{i}*S1(i))
end
grid on
grid minor
title('Espectros SRSS')
legend('NEC15','s368','s540','s741','s752','s753','s953','s960','s1004','s1080','s1048','s4886')
xlabel('Periodo $[s]$')
ylabel('PSa [cm/s^2]')
%% SRSS*S1 Promedio
figure(5)
plot(T,SaNEC,'LineWidth',lineWidth)
hold on
SRSS_S1={0};
SRSSavg=zeros(length(SaNEC),1);
for i=1:n
    SRSS_S1{i}=S1(i)*SRSS{i};
    SRSSavg = SRSSavg + SRSS_S1{i};
end
SRSSavg = SRSSavg/n;
plot(T,SRSSavg)
grid on
grid minor
title('Espectros Promedio')
legend('NEC15','SRSS\cdot S1_{avg}')
xlabel('Periodo [s]')
ylabel('PSa [cm/s^2]')
 
%% Relacion espectral
ratio=SaNEC./SRSSavg;
figure(6)
plot(Tinf:1/100:Tsup,ratio(index_inf:index_sup),'LineWidth',1)
xlim([-inf inf])
ylim([-inf inf]) 
grid on
grid minor
title('Relacion Espectral entre 0.2T y 1.5T')
xlabel('Periodo [s]')
ylabel('Relaciones Espectrales [cm/s^2]')
 
%% Factor S2
S2=max(ratio(index_inf:index_sup));
display(S2)
%% Factor SS
SS=S1*S2;
display(SS)
%% SRSS_S1 promedio * SS 
figure(7)
plot(T,SaNEC,'LineWidth',lineWidth)
hold on
SRSS_SS={0};
SRSSavg_SS=zeros(length(SaNEC),1);
for i=1:n
    SRSS_SS{i}=SS(i)*SRSS{i};
    SRSSavg_SS = SRSSavg_SS + SRSS_SS{i};
end
SRSSavg_SS = SRSSavg_SS/n;
plot(T,SRSSavg_SS)
grid on
grid minor
title('Espectro SRSS promedio Escalado con SS')
legend('NEC15','$SRSS\cdot SS_{avg}$')
xlabel('Periodo [s]')
ylabel('PSa [cm/s^2]')
 
%% Escalamiento de registros sismicos
SaEsc = SaTrat;
for i=1:n
    SaEsc{1,i} = SaTrat{1,i} * SS(i); % Escalar componente NS
    SaEsc{2,i} = SaTrat{2,i} * SS(i); % Escalar componente EW
end
 
dlmwrite('SaEsc11.txt',SaEsc{1,1}./100);
dlmwrite('SaEsc21.txt',SaEsc{1,2}./100);
dlmwrite('SaEsc31.txt',SaEsc{1,3}./100);
dlmwrite('SaEsc41.txt',SaEsc{1,4}./100);
dlmwrite('SaEsc51.txt',SaEsc{1,5}./100);
dlmwrite('SaEsc61.txt',SaEsc{1,6}./100);
dlmwrite('SaEsc71.txt',SaEsc{1,7}./100);
dlmwrite('SaEsc81.txt',SaEsc{1,8}./100);
dlmwrite('SaEsc91.txt',SaEsc{1,9}./100);
dlmwrite('SaEsc10_1.txt',SaEsc{1,10}./100);
dlmwrite('SaEsc11_1.txt',SaEsc{1,11}./100);
 
 
dlmwrite('SaEsc12.txt',SaEsc{2,1}./100);
dlmwrite('SaEsc22.txt',SaEsc{2,2}./100);
dlmwrite('SaEsc32.txt',SaEsc{2,3}./100);
dlmwrite('SaEsc42.txt',SaEsc{2,4}./100);
dlmwrite('SaEsc52.txt',SaEsc{2,5}./100);
dlmwrite('SaEsc62.txt',SaEsc{2,6}./100);
dlmwrite('SaEsc72.txt',SaEsc{2,7}./100);
dlmwrite('SaEsc82.txt',SaEsc{2,8}./100);
dlmwrite('SaEsc92.txt',SaEsc{2,9}./100);
dlmwrite('SaEsc10_2.txt',SaEsc{2,10}./100);
dlmwrite('SaEsc11_2.txt',SaEsc{2,11}./100);
