clear all;
close all;

result = load('allcollate_solution.dat');%loadind all collated steady state solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%since each row contain one node like row 1 is one gene, row 2 is another
miR200 = result(:,2);
ZeB1 = result(:,1);
PDL1 = result(:,3);
CDH1 = result(:,4);
SLUG = result(:,5);
%GRHL2 = result(:,6);
%SNAIL = result(:,6);
%NRF2 = result(:,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------------------------------------------------------------

Z_score_miR200 = zscore(miR200);  %zscore( X )%doing Zscore normalization values will rescale and will varry between -3 to 3
Z_score_ZeB1 = zscore(ZeB1);
Z_score_PDL1 = zscore(PDL1);
Z_score_CDH1 = zscore(CDH1);
Z_score_SLUG = zscore(SLUG);
%Z_score_GRHL2 = zscore(GRHL2);
%Z_score_SNAIL = zscore(SNAIL);
%Z_score_NRF2 = zscore(NRF2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z_score = [Z_score_ZeB1 Z_score_miR200 Z_score_SLUG Z_score_CDH1 Z_score_PDL1];
%save('z_score.txt','Z_score','-ascii','-tabs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%--------------------------------------------------

%EMT score calculation
EMTscore1=(Z_score_ZeB1+Z_score_SLUG-(Z_score_miR200+Z_score_CDH1))./4; 
%EMTscore2=Z_score_ZeB1-Z_score_miR200;  %ZEB1-miR200
%EMTscore3=Z_score_ZeB1-(Z_score_miR200+Z_score_CDH1);

%--------------------------------------------------
collate_zscore_emtscore=[Z_score EMTscore1];
save('col_zscore_emtscore.txt','collate_zscore_emtscore','-ascii','-tabs')
%------------------------------------------------------


%=========================
%-------figure 1B---------
%=========================
%lets plot the bar graph distribution of EMT score with kernel fit.

figure()
h = histfit(EMTscore1,50,'kernel') %loading EMTscore1
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', (yt/numel(EMTscore1))*10);
h = findobj(gca,'Type','line');
x=get(h,'Xdata');
y=get(h,'Ydata');
Y = y/numel(EMTscore1);
xlabel('EMTscore=((ZeB1+SLUG)-(miR200+CDH1))/4')
%%for changing axis properties

ax=gca;
%set axis color to dark black
ax.YColor = 'k';
ax.XColor = 'k';
% Set x and y font sizes.
%%ax.XAxis.FontSize = 20;
%%ax.YAxis.FontSize = 20;
% The below would set everything: title, x axis, y axis, and tick mark label font sizes.
ax.FontSize = 14;
% Bold all labels.
%ax.FontWeight = 'bold';
%fig.PaperPositionMode = 'auto';
% print('EMTscore_distribution','-dpng','-r300')%300 is the resolution in dpi
% print('EMTscore_distribution','-depsc','-tiff','-r300')%300 is the resolution in dpi
% savefig('EMT_score_distribution.fig')
% 





%=======================================================
%-----------scatter plot figure 1C------------------
%=======================================================

sz=3;
figure()
scatter(EMTscore1,Z_score_PDL1,sz,'k','filled');
[RHO,PVAL] = corr(EMTscore1,Z_score_PDL1,'Type','Spearman');
ylabel('PD-L1','FontSize', 17);
xlabel('EM score','FontSize', 17);
%%for changing axis properties
ax=gca;
%set axis color to dark black
ax.YColor = 'k';
ax.XColor = 'k';
% Set x and y font sizes.
%%ax.XAxis.FontSize = 20;
%%ax.YAxis.FontSize = 20;
% The below would set everything: title, x axis, y axis, and tick mark label font sizes.
ax.FontSize = 17;
% Bold all labels.
ax.FontWeight = 'bold';
%savefig('scatterplot_emscore_pdl1.fig')


%=========================================
%PD-L1 postive conditional frequency plot
%=============================
%-----figure 1D---------------
%=============================
data2=readmatrix('col_zscore_emtscore.txt');
countt1=0;
countt2=0;
countt3=0;
H_PDL1_pos_count=0;
M_PDL1_pos_count=0;
E_PDL1_pos_count=0;
H_PDL1_neg_count=0;
M_PDL1_neg_count=0;
E_PDL1_neg_count=0;
for i=1:numel(data2(:,1))
    p=data2(i,6);
    if (p<-0.25)
        countt1=countt1+1;
        %for j=1:numel(data1(1,:))
         %   Emat(count1,j)=data1(i,j);
        %end
        E_PDL1=data2(i,5)
        if (E_PDL1>0)
            E_PDL1_pos_count=E_PDL1_pos_count+1
        else
            E_PDL1_neg_count=E_PDL1_neg_count+1
        
    elseif (p>0.5)
        countt2=countt2+1;
        %for j=1:numel(data1(1,:))
         %   Mmat(count2,j)=data1(i,j);
        %end
        M_PDL1=data2(i,5)
        if (M_PDL1>0)
            M_PDL1_pos_count=M_PDL1_pos_count+1
        else
            M_PDL1_neg_count=M_PDL1_neg_count+1
        
    else
        countt3=countt3+1;
        %for j=1:numel(data1(1,:))
         %   Hmat(count3,j)=data1(i,j);
        %end
        H_PDL1=data2(i,5)
        if (H_PDL1>0)
            H_PDL1_pos_count=H_PDL1_pos_count+1
        else
            H_PDL1_neg_count=H_PDL1_neg_count+1
        
    

total = countt1 + countt2 + countt3;
frequency_Ephenotype=countt1./total;
frequency_Mphenotype=countt2./total;
frequency_Hphenotype=countt3./total;
collfreq =[frequency_Ephenotype frequency_Hphenotype frequency_Mphenotype];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
frequency_E_PDL1_pos = E_PDL1_pos_count./countt1;
frequency_M_PDL1_pos = M_PDL1_pos_count./countt2;
frequency_H_PDL1_pos = H_PDL1_pos_count./countt3;

frequency_E_PDL1_neg = E_PDL1_neg_count./countt1;
frequency_M_PDL1_neg = M_PDL1_neg_count./countt2;
frequency_H_PDL1_neg = H_PDL1_neg_count./countt3;

collfreq_PDL1 = [frequency_E_PDL1_pos frequency_E_PDL1_neg; frequency_H_PDL1_pos frequency_H_PDL1_neg; frequency_M_PDL1_pos frequency_M_PDL1_neg];
%bar(collfreq_PDL1)
collfreq_PDL1_pos = [frequency_E_PDL1_pos; frequency_H_PDL1_pos; frequency_M_PDL1_pos];
%save('collfreq_PDL1_pos.txt','collfreq_PDL1_pos','-ascii','-tabs')
figure;
bar(collfreq_PDL1_pos)
xlabel('Phenotype')
ylabel('Conditional Probability')



