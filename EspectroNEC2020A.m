function [Sa,T]=EspectroNEC2020A(Fa,Fd,Fs,r,eta,dT,Tmax,zg)
% Fa,Fd,Fs:  Factores de amplificacione del suelo
% r,eta,z:   Factores que dependen de la ubicacion de la estructura
% Td:        Periodo maximo del expecto
% dT:        Intrvalo de periodos para el espectro de diseño
% Este programa devuelve un valor de Sa(i) para cada valor de T(i)
T=(0:dT:Tmax); 
n=length(T);
Sa=zeros(n,1);

T0=0.1*Fs*Fd/Fa;
Tc=0.55*Fs*Fd/Fa;
for i=1:n
    if  T(i)<T0
        Sa(i)=zg*Fa*(1+(eta-1)*T(i)/T0);
    elseif T(i)>Tc
        Sa(i)=eta*zg*Fa*(Tc/T(i))^r;
    else 
        Sa(i)=eta*zg*Fa;
    end
end
        

