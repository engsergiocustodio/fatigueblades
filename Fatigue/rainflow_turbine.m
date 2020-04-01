%--------------------------------------------------------------------
%%%% "rainflow_turbine.m" FUNCTION
%%%% IMPLEMENTED BY SERGIO CUSTODIO - engsergiocustodio@gmail.com
%%%% CONTRIBUTIONS OF PROF. DR. LEONARDO D. RODRIGUES
%--------------------------------------------------------------------

function [rf] = rainflow_turbine(tensao)

n_even=length(tensao);
ident=zeros(n_even,1);
%% Identificar se evento de inicio é pico ou vale

if tensao(1,:)<tensao(2,:)
    ident1=0;
else
    ident1=1;
end

%% Identificar quem é pico e quem é vale na história

for i=2:n_even
    k=i-1;
    if ident(k,:)==ident1 
      ident(i,:)=1; 
    else 
      ident(i,:)=0;
    end 
end

cont_ant=zeros(n_even,1);
t_max=zeros(n_even,1);
t_min=zeros(n_even,1);
t_maxi=zeros(n_even,1);
t_mini=zeros(n_even,1);

for i=1:(n_even-1)
   if ident(i,:)==1
      t_max(i,:)=tensao(i,:);
      t_maxi(i,:)=i;
   else
      t_min(i,:)=tensao(i,:);
      t_mini(i,:)=i;
   end

   for j=(i+1):n_even
     sai=0;
     if j==(i+1)
         max_min=tensao(j,:);
         mm_j=j; 
     end
     
     %verifica contagem anterior
     
     if ((ident(i,:)==0) & (ident(j,:)==1))
        if (cont_ant(j,:)~=isempty(cont_ant(j,:))) % não vazio
           temp=max_min;
           temp1=mm_j;
           max_min=cont_ant(j,:);
           mm_j=ant_j(j,:);
           sai=1;
        end
            if((j>(i+1)) & (tensao(j,:)>=max_min)) 
                if sai==1
                cont_ant(j,:)=temp;
                ant_j(j,:)=temp1;
                else 
                cont_ant(j,:)=max_min;
                ant_j(j,:)=mm_j;
                max_min=tensao(j,:);
                mm_j=j;
                end
            end
     end
     
    if((ident(i,:)==1) & (ident(j,:)==0))
        if (cont_ant(j,:)~=isempty(cont_ant(j,:))) % diferente a vazio
           temp=max_min; 
           temp1=mm_j; 
           max_min=cont_ant(j,:); 
           mm_j=ant_j(j,:); 
           sai=1; 
        end
        if((j>(i+1)) & (tensao(j,:)<=max_min)) 
          if sai==1
            cont_ant(j,:)=temp; 
            ant_j(j,:)=temp1;
          else 
            cont_ant(j,:)=max_min; 
            ant_j(j,:)=mm_j; 
            max_min=tensao(j,:); 
            mm_j=j;
          end
        end
    end
     
     % Se o pico ou vale é maior ou menor
     
     if((ident(i,:)==1) & (ident(j,:)==1) & (tensao(j,:)>=t_max(i,:)))
         sai=1;
     end
     
     if((ident(i,:)==0) & (ident(j,:)==0) & (tensao(j,:)<=t_min(i,:)))
         sai=1;
     end
     
     % Saida
     
     if sai==1 
         break
     end
   end
   
   if ident(i,:)==0
       t_max(i,:)=max_min; 
       t_maxi(i,:)=mm_j; 
   else 
       t_min(i,:)=max_min; 
       t_mini(i,:)=mm_j; 
   end
   
end
 
tensao_a=abs(t_max-t_min)./2;
tensao_m=(t_max+t_min)./2;

rf=[tensao_a';tensao_m';ones(1,size(tensao_m,1))*0.5];