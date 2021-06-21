clc;
clear all;
close all;
%===========================================
%-----------figure 2D,E---------------------
%===========================================

%format in which matrix are stored
%Z_score = [Z_score_ZeB1 Z_score_miR200 Z_score_SLUG Z_score_CDH1 Z_score_PDL1];
%EMTscore1=(Z_score_ZeB1+Z_score_SLUG-(Z_score_miR200+Z_score_CDH1))./4; 
%collate_zscore_emtscore=[Z_score EMTscore1];
%-----------------------------------------------------------
%loading data
data1=readmatrix('col_zscore_emtscore.txt');
count1=0;
count2=0;
count3=0;
for i=1:numel(data1(:,1))
    p=data1(i,7);
    if (p<-0.25) % EM score<-0.25---> E phenotype
        count1=count1+1;
        for j=1:numel(data1(1,:))
            Emat(count1,j)=data1(i,j);
        
    elseif (p>0.5) % EM score<0.5---> M phenotype
        count2=count2+1;
        for j=1:numel(data1(1,:))
            Mmat(count2,j)=data1(i,j);
        
    else           % EM score between -0.25 to 0.5---> H phenotype
        count3=count3+1;
        for j=1:numel(data1(1,:))
            Hmat(count3,j)=data1(i,j);
        
    

total = count1 + count2 + count3;
frequency_Ephenotype=count1./total;
frequency_Mphenotype=count2./total;
frequency_Hphenotype=count3./total;
collfreq =[frequency_Ephenotype frequency_Hphenotype frequency_Mphenotype];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------------------------------------------------------------
figure()
%for ploting bar graph
bar(collfreq)
y = collfreq;
grid on
ylabel('Frequency')
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

%no. of X tick
ax.XTick = [1 2 3]; 
%label X tick
ax.XTickLabels = {'Epithelial','Hybrid','Mesenchymal'}
%rotate X label by 45 degree
ax.XTickLabelRotation = 45;
%for getting y axis value at the tip of bar graph
text(1:length(y),y,num2str(y'),'vert','bottom','horiz','center'); 
box off
% savefig('barplot.fig')
% print('barplot','-dpng','-r300')%300 is the resolution in dpi
% print('barplot','-depsc','-tiff','-r300')%300 is the resolution in dpi