
%% recuparation des données à partir du fichier csv
IEX = readtable('TB_Ferreira2020_FS_MLneuron_PLS_FoM.csv','Range','A3:A311');
 

Iex=table2array(IEX)% transformer le tableau en un vecteur

PRMS = readtable('TB_Ferreira2020_FS_MLneuron_PLS_FoM.csv','Range','B3:B311');

Prms=table2array(PRMS)%transformer le tableau en un vecteur

FSPIKE = readtable('TB_Ferreira2020_FS_MLneuron_PLS_FoM.csv','Range','C3:C311');

 Fspike=table2array(FSPIKE)%transformer le tableau en un vecteur
 %% interpolation 
    Iex_reduit=Iex.*10^9; % changement d'échelle
    Fspike_reduit=Fspike.*10^(-3);% changement d'échelle
    Prms_reduit=Prms.*10^(9)
    
P=polyfit(Iex_reduit,Fspike_reduit,15)%interpolation de Fspike(Iex)
P2= polyfit(Iex_reduit,Prms_reduit,15)%interpolation de Prms(Iex)




Fspike2= polyval(P,Iex_reduit).*10^(3)% calcul de nos valeurs de fspike 
Prms2=polyval(P2,Iex_reduit).*10^(-9)% calcul de nos valeurs de Prms 


%% Estimation de l'erreur RMS De notre interpolation de Prms(Iex)
% normalise les valeurs 
%  Prms2_normalise=Prms2.*10^12; 
%  Prms_normalise=Prms.*10^12;

N=size(Prms2,1)          %nombre d'élements dans nos deux vecteurs 


s = 0; %somme

for i=1:1:N
    s = s + (Prms2(i) - Prms(i))^2 ;
end
s= s./N
s=sqrt(s)
RMS_Prms = s 
%% Estimation de l'erreur RMSPE (RMS en pourcentage) De notre interpolation de Prms(Iex)
N=size(Prms2,1)          %nombre d'élements dans nos deux vecteurs 
% normalise les valeurs 
%  Prms2_normalise=Prms2.*10^12; 
%  Prms_normalise=Prms.*10^12;

s = 0; %somme

for i=1:1:N
    s = s + ((Prms(i) - Prms2(i))./(Prms(i)))^2 ;
end
s=100*s./N
s=sqrt(s)
RMSPE_Prms = s 
%% Estimation de l'erreur RMS De notre interpolation de Fspike(Iex)
% Fspike2_normalise=Fspike2.*10^(-3); 
% Fspike_normalise=Fspike.*10^(-3);

N=size(Fspike2,1)          %nombre d'élements dans nos deux vecteurs 

s = 0; %somme

for i=1:1:N
    s = s + ((Fspike2(i) - Fspike(i))^2) ;
end

s= 100*s./N
s=sqrt(s)

RMS_Fspike = s 


%% Estimation de l'erreur RMSPE (RMS en pourcentage) De notre interpolation de Fspike(Iex)
N=size(Prms2,1)          %nombre d'élements dans nos deux vecteurs 
% Fspike2_normalise=Fspike2.*10^(-3); 
% Fspike_normalise=Fspike.*10^(-3);

s = 0; %somme

for i=1:1:N
    s = s + ((Fspike(i) - Fspike2(i))./(Fspike(i)))^2 ;
end
s= s./N
s=sqrt(s)
RMSPE_Fspike = s 





%% BILAN des erreurs 

RMS_Fspike 
RMS_Prms
Vecteur_erreur_percentage_Fspike_Prms=[RMSPE_Fspike,RMSPE_Prms]

%% test de nos polynomes
% REMARQUE : L'on constate que notre interpolation atteint ses limites (
% pas forcément pertinente) 

%Test valeur présente dans le tableau 
Iex_test=(3E-12).*10^9

Fspike_test= polyval(P,Iex_test).*10^(3)% calcul de nos valeurs de fspike de test
Prms_test=polyval(P2,Iex_test).*10^(-9)

% Valeur non présente dans le tableau
Iex_test2=(2.67E-09).*10^9

Fspike_test2= polyval(P,Iex_test2).*10^(3)% calcul de nos valeurs de fspike de test
Prms_test2=polyval(P2,Iex_test2).*10^(-9)



%% figures DE COMPARAISON AVEC LE MODELE CADENCE

figure(1)
plot( Iex,Prms2)
title('LIF', 'fontsize', 14)

    xlabel('Iex (A)', 'fontsize', 12);

    ylabel('Prms (W)', 'fontsize', 12)
    hold on 
    
plot(Iex,Prms,'r')
 title('LIF', 'fontsize', 14)

    xlabel('Iex (A)', 'fontsize', 12);

    ylabel('Prms (W)', 'fontsize', 12)
    legend ('modele lif matlab','modele lif cadence')
hold off
% 
% Iex2=Iex(1:25)
% Prms2=Prms(1:25)
 

figure(2)
plot( Iex,Fspike2)
 
    title('LIF', 'fontsize', 14)

    xlabel('Iex (A)', 'fontsize', 12);

    ylabel('Fspike (Hz)', 'fontsize', 12)
    
    hold on

plot( Iex,Fspike)
  title('LIF', 'fontsize', 14)

    xlabel('Iex (A)', 'fontsize', 12);

    ylabel('Fspike (Hz)', 'fontsize', 12)
    legend ('modele lif matlab','modele lif cadence')
    
    
    
hold off

%% Calcul de notre potentiel d'action pour un IEX1
% calcul d'une certaine valeur de fspike en fonctio d'un certain Iex

Iex1= (3E-9).*10^7
Fspike_final1 = polyval(P,Iex1).*10^(3)



[Vtemps,AMPFS]=LIF(0.5,2000,Fspike_final1,90)

 
%% Affichage de notre potentiel d'action pour un IEX1

%Afin d'afficher une constante on doit passer par la fonction impulse qui
%utilise le nombre complexe P

p=tf('p');
Iex1=Iex1*10*10;

 figure (7)
    subplot(2,1,1);
    plot(Vtemps,AMPFS)
    set(gcf,'color','w');

    
    title('LIF', 'fontsize', 14)

    xlabel('Time (ms)', 'fontsize', 12);

    ylabel('Membrane Potential (mV)', 'fontsize', 12)

    
subplot(2,1,2); 
    impulse(Iex1/p)
    set(gcf,'color','w');
   

    
    title('LIF', 'fontsize', 14)

    xlabel('Time (ms)', 'fontsize', 12);

    ylabel('Iex (nA)', 'fontsize', 12)

    

% 4.4998e+04

% xlsread('suite1.csv','A12:A320');
% 
% 
% 
% Fspike_final
%% Calcul de notre potentiel d'action pour un IEX1
% calcul d'une certaine valeur de fspike en fonctio d'un certain Iex

Iex2= (1E-9).*10^7
Fspike_final2 = polyval(P,Iex2).*10^(3)



[Vtemps,AMPFS]=LIF(0.5,2000,Fspike_final2,90)
% Affichage du potentiel d'action pour un Iex2
%Afin d'afficher une constante on doit passer par la fonction impulse qui
%utilise le nombre complexe P

Iex2=Iex2*10*10;

 figure (8)
    subplot(2,1,1);
    plot(Vtemps,AMPFS)
    set(gcf,'color','w');

   

    
    title('LIF', 'fontsize', 14)

    xlabel('Time (ms)', 'fontsize', 12);

    ylabel('Membrane Potential (mV)', 'fontsize', 12)

    subplot(2,1,2); 
    impulse(Iex2/p)
    set(gcf,'color','w');

   

    
    title('LIF', 'fontsize', 14)

    xlabel('Time (ms)', 'fontsize', 12);

    ylabel('Iex (nA)', 'fontsize', 12)

    

% 4.4998e+04
% 
% xlsread('suite1.csv','A12:A320');

%% Comparaison avec le modele cadence 

Iex1= (30E-12).*10^9
Fspike_final1 = polyval(P,Iex1).*10^(3)

[Vtemps,AMPFS]= LIF(0.5,2000,Fspike_final1,90)

% transformer le tableau en un vecteur

Vtemps_cadence =readtable('Vout_LIF_30p','Range','A2:A12823');

Vtemps_cadence =table2array(Vtemps_cadence );% transformer le tableau en un vecteur


 figure (9)
    subplot(2,1,1);
    set(gcf,'color','w');

    plot(Vtemps,AMPFS,'b')
hold on
    plot(Vtemps_cadence*10^3,Vout_cadence*10^3,'r')
     title('LIF', 'fontsize', 14)

    xlabel('Time (ms)', 'fontsize', 12);

    ylabel('Spike (mV)', 'fontsize', 12)
    legend ('modele lif matlab','modele lif cadence')
    
    subplot(2,1,2); 
    impulse(Iex1/p)
    set(gcf,'color','w');
    title('LIF', 'fontsize', 14)

    xlabel('Time (ms)', 'fontsize', 12);

    ylabel('Iex (nA)', 'fontsize', 12)
   
%% ERREUR RMSPE entre les potentiels d'action simulés sous matlab et cadences
T1=0.318
T2=0.374
i1= T1*2000
i2=T2*2000
A=AMPFS(i1)
A=AMPFS(i2)
    
 A= Vout_cadence(i1)
  A=Vout_cadence(i2)
    
for i=i1:1:i2
    s = s + ((Vout_cadence(i)*10^(4) - AMPFS(i))./(Vout_cadence(i)*10^(4)))^2 ;
end
s= s./N
s=sqrt(s)
RMSPE_AMPFS = s     
    
 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




