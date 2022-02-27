function[t_recort,signal_recort]=Intensidad_Areas(signal,paso,ind)
% ia=intensidad areas
% tt tiempo recortado correspondiente a la señal recortada
 
% clear all
% clc 
% C18OE=importdata('2015_09_16_C18O_E.txt');
% signal=C18OE;
% paso=0.01;
 
l=length(signal);                   % Longitud de la senal
dt=paso;                            % Intervalo de tiempo
t=0:dt:(l-1)*paso;                  % Vector de tiempo
 
ac2=signal.^2;                      % Aceleracion al cuadrado
integrar=0;                         % Variable para intensidad de arias
acumulada=0;                        % Variable para Ia acumulada
 
for i=2:l
integrar=integrar+(ac2(i-1)+ac2(i))/2*dt;        % Integrar
acumulada(i-1)=pi/(2*981)*integrar;              % Ia acumulada
end
 
ia=pi/(2*981)*integrar;                          % Intensidad de arias
acumulada=acumulada./ia*100;                     % Porcentaje
acumulada=acumulada';
 
a=length(acumulada);
signal_recort=0;
for i=1:a
% recorta y almacena en una nueva   
if acumulada(i)>=5 && acumulada(i)<=95 
  	% variable los valore significativos de la seña 
        signal_recort(i)=signal(i);     
    end  
end   
 
signal_recort=signal_recort';              
% Elimina los valores = o del vector señal acumulada
signal_recort(signal_recort==0)=[];       
aa=length(signal_recort);
t_recort=(0:paso:(aa-1)*paso);

