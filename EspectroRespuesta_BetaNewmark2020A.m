function[Sa,Sv,Sd,T,Sdr,Saa]=EspectroRespuesta_BetaNewmark2020A(sg,dT,Tmax,z,beta,gamma,dt)
% Funsion para encontrar el esprectro de respuesta elastica de un sistema  % simple 
% PARAMETROS DE ENTRADA
% sg: aceleracion del suelo 
% dT: intervalo de periodo del espectro de respuesta 
% Tmax: Periodo maximo del espectro de respuesta 
% z: Factor de amortiguamiento viscoso 
% beta:  1/8<=beta<=1/4
% gamma: 0<=gamma<=1 (usualmente se toma 1/2)
% dt: intervalo de paso del tiempo del acelerograma
 
% PARÁMETROS DE SALIDA
% Saa: Pseudo Aceleracion
 
T=(0:dT:Tmax)';
p=length(T);
Sa=zeros(p,1);
Sv=zeros(p,1); 
Sd=zeros(p,1);
Sdr=zeros(p,1);
Saa=zeros(p,1);
n=length(sg);
t=(0:dt:dt*(n-1))';
 
PGA=max(abs(sg));
Sa(1)=PGA;
Saa(1)=PGA;
for j=1:(p-1)
    wn=2*pi/T(j+1);
    Ma_m=1+2*z*wn*gamma*dt+beta*wn^2*dt^2;
    d=zeros(n,1);
    v=zeros(n,1);
    a=zeros(n,1);
    at=zeros(n,1);
    for i=1:(n-1)
        incsg=sg(i+1)-sg(i);
        incfa_m=-incsg-wn^2*v(i)*dt-a(i)*(2*z*wn*dt+wn^2*dt^2/2);
           
        inca=incfa_m/Ma_m;
        incv=a(i)*dt+inca*gamma*dt;
        incd=v(i)*dt+a(i)*dt^2/2+beta*inca*dt^2;
 
        a(i+1)=a(i)+inca;
        v(i+1)=v(i)+incv;
        d(i+1)=d(i)+incd;
        at(i+1)=a(i+1)+sg(i+1);
    end
    Sv(j+1)=max(abs(v));
    Sa(j+1)=wn*Sv(j+1);
    Sd(j+1)=Sv(j+1)/wn;
    Sdr(j+1)=max(abs(d));
    Saa(j+1)=max(abs(at));
end