function [t,dc,vc,ac]=Corec_Filt_2020A(Sa,paso,ind)
% --------------------- CORRECCION A LA LIENA BASE ------------------------ 
 
a0=Sa;% Registro Original 
 
dt=paso;
n=length(a0);
t=(0:dt:(n-1)*dt)';
 
% -------INTEGRACIÓN DE LA SEÑAL (Método de los Trapecios) PARA OBTENER LA % VELOCIDAD y DESPLAZAMIENTO. 
v0=cumtrapz(a0)*dt;
d0=cumtrapz(v0)*dt;
% ------------------ QUITAR LA LÍNEA DE TENDENCIA ------------------------- 
adetrend=detrend(a0);
vdetrend=cumtrapz(adetrend)*dt;
ddetrend=cumtrapz(vdetrend)*dt; 
% ------------------- FILTRADO DE SEÑAL PASA BANDAS ----------------------- 
flc=0.1;%Hz
fhc=30; %Hz 
Fn=1/2/dt; %frecuencia Nyquist
[fb,fa]=butter(4,[flc/Fn,fhc/Fn]);
vc=filter(fb,fa,vdetrend);
dc=cumtrapz(vc)*dt;
ac=diff(vc)/dt;
ac(end+1)=0;

