function [Vtemps,AMPFS,Lambda1,Lambda2,T] = FS(Iex,DC,Fe,Fsp,Vpp)

 
%% Initialisation des paramètres 
 nb_spike=5; %nombre de pic

Te=1/Fe; %période d'Échantillonnag

T =1/(3.40*Fsp);% période du signal
T=T*10^4
 

Vtemps=[0:Te:1]%Vecteur de temps avec un pas d'échantillonnage de Te

Ton = DC*T;%Durée de l'exponentielle croissante

 

 



Lambda1= log(Vpp)/Ton;%Facteur de croissance 

 


Lambda2= log(Vpp)/(T-Ton);% Facteur d'atténuation

%% Optimisation
A1=0.01  % controler l'allure des exponentielles 
A2=100
Lambda1=log(1+(Vpp/A1))/Ton

Lambda2=log(1-(Vpp/A2))/(Ton-T)

B1=-A1
B2=-A2
 
%% Fonction principale
AMPFS= zeros(1,length(Vtemps)) %vecteur d'amplitudes de FS

 

 

    x = 0;

    y = Ton;

    z = T;

    a=1;

   

        

    for i=1:1:length(Vtemps) % Boucle pour passer par toutes les valeurs du temps

 

        if (Vtemps(i)>=y && AMPFS(i)==0)

            if(a==1)

                a=Vtemps(i)

            end

             x = x+a;%préparation de la condition qui suit notre spike actuelle<

             y = y+a;

             z = z+a;

            

             %T=T+T; % on augmente la valeur de notre période

         end

            

          if (Vtemps(i)>=(x) && Vtemps(i)<=y)   %Condition pour l'exponentielle croissante

              AMPFS(i) = A1*exp(Lambda1*((Vtemps(i))-x))+B1           % Le facteur B1 fait passer successivement l'exponentielle croissante à zero lorsque Vtemps=0

     

          elseif (Vtemps(i)>y && Vtemps(i)<=z)%Condition pour l'exponentielle décroissante   

              AMPFS(i) =A2*exp(-Lambda2*((Vtemps(i)-x) ))+B2         % Le facteur B1 fait passer successivement l'exponentielle decroissante à Vpp lorsque Vtemps=Ton

          end 

 

    end

   

    

 

 

end