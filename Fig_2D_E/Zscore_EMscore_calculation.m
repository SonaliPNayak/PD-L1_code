clear all;
close all;

result = load('allcollate_solution.dat');%loading all collated steady state solution 
%node representation
%1	ZEB1
%2	miR200
%3	PD-L1
%4	CDH1
%5	SLUG

%%%%%%%%%%%%%%%%%
%since each row contain one node like row 1 is one gene, row 2 is another
miR200 = result(:,2);
ZeB1 = result(:,1);
PDL1 = result(:,3);
CDH1 = result(:,4);
SLUG = result(:,5);
%GRHL2 = result(:,6);
%SNAIL = result(:,6);
%NRF2 = result(:,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% VERY VERY IMPORTANT----->>> IN CASE OF OE/DE THE MEAN AND STDEVIATION OF THE "CORE" TAKEN TO FIND ZSCORE OF OE/DE GENES 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%z score normalisation
Z_score_miR200 = (miR200 - miR200_mean)./(miR200_standarddev );
Z_score_ZeB1 = (ZeB1 - ZeB1_mean)./(ZeB1_standarddev );
Z_score_PDL1 = (PDL1 - PDL1_mean)./(PDL1_standarddev );
Z_score_CDH1 = (CDH1 - CDH1_mean)./(CDH1_standarddev );
Z_score_SLUG = (SLUG -SLUG_mean)./(SLUG_standarddev );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%storing data in matrix in txt file
Z_score = [Z_score_ZeB1 Z_score_miR200 Z_score_SLUG Z_score_CDH1 Z_score_PDL1];
save('z_score.txt','Z_score','-ascii','-tabs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--------------------------------------------------

%EMT score calculation
EMTscore1=(Z_score_ZeB1+Z_score_SLUG-(Z_score_miR200+Z_score_CDH1))./4; 
%--------------------------------------------------
collate_zscore_emtscore=[Z_score EMTscore1];
%storing data in matrix in txt file
save('col_zscore_emtscore.txt','collate_zscore_emtscore','-ascii','-tabs')
%------------------------------------------------------