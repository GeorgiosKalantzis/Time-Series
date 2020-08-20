
fileID = fopen('C:\Users\law_g1\Desktop\dat1.dat'); % Specify full path dat1.dat
line = fgetl(fileID);
i = 0;
V = zeros(3521,1);

while ischar(line)
    V(i+1,1) = str2double(line);
    line = fgetl(fileID);
    i = i + 1; 
end
fclose(fileID);

ed = 1;
st = 1;
K = zeros(14,250);

%---------------------------------
% Passing the values of dat1 file
%---------------------------------
for i = 1:14
   K(i,:) = V(st:250*ed);
   st = 250*ed + 1;
   ed = ed + 1;
   if K(i,250) == 0
       newsize = i - 1;
      break;
   end
    
end

sz = size(K);

while i < sz(1,1)
    sz = size(K);
    K(i,:) = []; 
end



%--------------------------------
% Remove trend with n-differences
%--------------------------------

Kn = zeros(newsize,249);

h=zeros(newsize,1);


for i=1:newsize
    for j=1:249
      Kn(i,j) = K(i,j+1) - K(i,j);
    end
    h(i,1)=adftest(Kn(i,:));
     if h(i,1)==0   %stationarity check
         for j=1:248
            Kn(i,j)=K(i,j+2)-2*K(i,j+1)+K(i,j);
        end
       
     end
     h(i,1)=adftest(Kn(i,:));
    if h(i,1)==0  %stationarity check
       
        for j=1:247
            Kn(i,j)=K(i,j+3)-3*K(i,j+2)+3*K(i,j+1)-K(i,j);
        end
    end
end




%---------------------------------------
%Computing and figureing autocorrelation
%---------------------------------------

for i = 1:newsize 
    %figure(i);
   %[hV,pV,QV,xautV]=portmanteauLB(Kn(i,:),248,0.124,'gg');
end

%-------------------------------------------
% Check for the best model and compute NRMSE
%-------------------------------------------


h=zeros(newsize,1);
for i=1:newsize
         acm= autocorrelation(Kn(i,:),30);
          pacf=acf2pacf(acm(:,2),1);
          ylim([-0.125 0.125])
          autocorr(Kn(i,:),'NumLags',30,'NumSTD',0.124)
          figure(i);
          parcorr(Kn(i,:),'NumLags',15,'NumSTD',0.127);
           
end


nrmseV=zeros(10,1);
thetaV=zeros(newsize,2);
aicS=zeros(10,1);
phiV=zeros(5,12);
nrmseV(1,1),phiV,thetaV(1,[1]),SDz,aicS(1,1),fpeS,armamodel]=fitARMA(Kn(11,:),0,1,1);


%---------------------------
%Prediction with ar(5) model
%---------------------------

nrmseV=zeros(newsize,1);
phiV=zeros(6,newsize);
preM=zeros(100,newsize);
preV=zeros(newsize,1);


for i=1:newsize
   [nrmseV(i,1),preM(:,i),phiV(:,i),thetaV]=predictARMAnrmse(Kn(i,:),5,0,1,100);
   preV(i,1) = predictARMAmultistep(Kn(i,:),249,5,0,1,'1');
end


%---------------
%Second section
%---------------

Z=zeros(1,459);
%creating our new time series
    for j=1:1:459
         Z(1,j)=V(6*j);
    end
    
   
    
%--------------------------------------
% Checking the best way to remove trend
%--------------------------------------
    
Xn=zeros(1,458);
for j=1:458
   Xn(1,j) = log(Z(1,j+1)) - log(Z(1,j));  
end
h=adftest(Xn(1,:));
portmanteauLB(Xn,457,0.05,'g');
%{
for j=1:458
        Zn(1,j) = (Z(1,j+1) - Z(1,j))/(Z(1,j+1));
       
end
h(2,1)=adftest(Zn(1,:));
for j=1:458
    Yn(1,j)=log(Z(1,j+1)) - log(Z(1,j));
    
end
h(3,1)=adftest(Yn(1,:));
%}  

%------------------------
%PART A
%-------------------------

[hV,pV,QV,xautV]=portmanteauLB(Xn(1,:),457,0.05,'gg');
autocorr(Xn(1,:),'NumLags',457,'NumSTD',1.98);
arcorr(Xn(1,:),'NumLags',100,'NumSTD',1.96);
autocorr(Xn(1,:),'NumLags',100,'NumSTD',1.96);
[nrmseV,phiV,thetaV,SDz,aicS,fpeS,armamodel]=fitARMA(Xn(1,:),5,0,1);
pre= predictARMAmultistep(Xn(1,:),458,5,0,1,'1');
[nrmseV,preM,phiV,thetaV] = predictARMAnrmse(Xn(1,:),0,0,1,183,'1');



    
    
