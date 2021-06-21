clc;

clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%
%mean and standarddev are taken from core circuit
miR200_standarddev =5.14;
ZeB1_standarddev = 9.7;
PDL1_standarddev = 2.96;
CDH1_standarddev = 6.56;
SLUG_standarddev = 5.34;


miR200_mean = 1.5;
ZeB1_mean = -3.51;
PDL1_mean = 4.03;
CDH1_mean =-0.38;
SLUG_mean = 1.53;
%%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('NonNegative',[1 2 3 4 5]);
var_cond1=100*rand(1,1000)
var_cond2=100*rand(1,1000)
var_cond3=100*rand(1,1000)
var_cond4=100*rand(1,1000)
var_cond5=100*rand(1,1000)

[t,y] = ode45(@immune_eqn,[0 200],[var_cond1(k) var_cond2(k) var_cond3(k) var_cond4(k) var_cond5(k)],options);

%z score normalisation
Z_score_miR200 = (log2(y(:,2)) - miR200_mean)./(miR200_standarddev );
Z_score_ZeB1 = (log2(y(:,1)) - ZeB1_mean)./(ZeB1_standarddev );
Z_score_PDL1 = (log2(y(:,3)) - PDL1_mean)./(PDL1_standarddev );
Z_score_CDH1 = (log2(y(:,4)) - CDH1_mean)./(CDH1_standarddev );
Z_score_SLUG = (log2(y(:,5)) -SLUG_mean)./(SLUG_standarddev );


%--------------------------------------------------
%save('EMTscore.txt','collate_zscore_emtscore','-ascii','-tabs')
if (EMTscore(200,1)>0.7784)
     count1=count1+1;
     plot(t,Z_score_PDL1,'Color','#FF7F7F','Linewidth',1);
     %ylim([-1.2 1])
   
     hold on;
elseif (EMTscore(200,1)>0) && (EMTscore(200,1)<0.2)
     count2=count2+1;
     plot(t,Z_score_PDL1,'Color','#FFD17F','Linewidth',3.5);
     %ylim([-1.2 1])
     hold on;
elseif (EMTscore(200,1)<-0.2)
     count3=count3+1;
     plot(t,Z_score_PDL1,'Color','#7F7F7F','Linewidth',1);
     %ylim([-1.2 1])
     hold on;
%     count_3(k)=count3;
else
     count4=count4+1;

%plotting
plot(t,Z_score_PDL1);