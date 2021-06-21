clear all;
close all;

result = load('allcollate_solution.dat');%loadind all collated steady state solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------------------------------------------------

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================
%-----figure S1C----
%===================

figure()
h = histfit(Z_score_PDL1,50,'kernel') %loading EMTscore1
yt = get(gca, 'YTick');
set(gca, 'YTick', yt, 'YTickLabel', (yt/numel(Z_score_PDL1))*10);
h = findobj(gca,'Type','line');
x=get(h,'Xdata');
y=get(h,'Ydata');
Y = y/numel(Z_score_PDL1);
xlabel('PDL1 distribution','FontSize', 17)

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

%title('EMT SLUG PDL1 NRF2','FontSize', 17)
%saveas(gcf,'EMTscore_distribution.tif')%this is another alternative option but does not gie good resolution
%fig = gcf;
%fig.PaperPositionMode = 'auto';
% print('PDL1','-dpng','-r300')%300 is the resolution in dpi
% print('PDL1','-depsc','-tiff','-r300')%300 is the resolution in dpi
% 
% savefig('PDL1_distribution.fig')

%=======================================================
%-----------scatter plot figure S1E------------------
%=======================================================
%======
figure()
scatter(Z_score_CDH1,Z_score_PDL1,sz,'k','filled');
[RHO,PVAL] = corr(Z_score_CDH1,Z_score_PDL1,'Type','Spearman');
xlabel('CDH1 expression','FontSize', 17);
ylabel('PD-L1 expression','FontSize', 17);

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
%savefig('scatterplot_cdh1_pdl1.fig')

%=========================
figure()
scatter(Z_score_ZeB1,Z_score_PDL1,sz,'k','filled');
[RHO,PVAL] = corr(Z_score_ZeB1,Z_score_PDL1,'Type','Spearman');
xlabel('ZEB1 expression','FontSize', 17);
ylabel('PD-L1 expression','FontSize', 17);

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
%savefig('scatterplot_ZeB1_pdl1.fig')

%==================
figure()
scatter(Z_score_SLUG,Z_score_PDL1,sz,'k','filled');
[RHO,PVAL] = corr(Z_score_SLUG,Z_score_PDL1,'Type','Spearman');
xlabel('SLUG expression','FontSize', 17);
ylabel('PD-L1 expression','FontSize', 17);

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
%savefig('scatterplot_SLUG_pdl1.fig')

%====================================================
%-------PCA(ward method)--=>figure S1A,D-------------
%====================================================

[wcoeff,score,latent,tsquared,explained] = pca(Z_score);
PC1 = score(:,1);
PC2 = score(:,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save('score.txt','score','-ascii','-tabs')
%save('wcoeff.txt','wcoeff','-ascii','-tabs')
figure()
scatter(score(:,1),score(:,2),3,'k','filled')
xlabel('PC-1(84.4% var.)')
ylabel('PC-2(7.4% var.)')
%xlabel('PC-1(84.4% var.)')
%ylabel('PC-2(7.4% var.)')
%title('Ward method clustering ')
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

%savefig('pcascatter_wards.fig')
%print('pcascatter_wards','-dpng','-r300')%300 is the resolution in dpi
%print('pcascatter_wards','-depsc','-tiff','-r300')%300 is the resolution in dpi
%==============================
%-------figure S1A-------------
%==============================
figure()
sz=3;
scatter(PC1,PC2,sz,EMTscore1,'filled');
c = colorbar;
set(c, 'ylim', [-3 3])
%c.Label.String = 'PD-L1';
c.Label.String = 'EMT score';
c.Label.FontSize = 17;

xlabel('PC-1(84.44% var.)')
ylabel('PC-2(7.36% var.)')
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
%savefig('pca_emtscore_colorbar.fig')
%============================
%------figure S1D------------
%============================
figure()
sz=3;
scatter(PC1,PC2,sz,Z_score_PDL1,'filled');
c = colorbar;
set(c, 'ylim', [-3 3])
%c.Label.String = 'PD-L1';
c.Label.String = 'EMT score';
c.Label.FontSize = 17;

xlabel('PC-1(84.4% var.)')
ylabel('PC-2(7.4% var.)')
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



