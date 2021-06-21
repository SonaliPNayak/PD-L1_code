clear all;
close all;
clc;
result1=readmatrix('PDL1_paraset4582_bistable_kmir1.txt');
%result1=readmatrix('PDL1_paraset4582_bistable_kmirplus.txt');

%f.write(str(j)+"\t"+str(ZEB1[i])+"\t"+str(miR200[i])+"\t"+str(PDL1[i])+"\t"+str(CDH1[i])+"\t"+str(SLUG[i])+"\n")
k_ZeB1 = result1(:,1);
% miR200 = result1(:,3);
% ZeB1 = result1(:,2);
 PDL1 = result1(:,4);
%CDH1 = result1(:,5);
%SLUG = result1(:,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extraction of k PDL1 dataset

stable_upper_PDL1 = PDL1(1:443,:);
unstable_PDL1 = PDL1(443:518,:);
stable_lower_PDL1 = PDL1(518:900,:);
% %------------------------------------
% %extraction of k zeb dataset
stable_upper_k_ZeB1 = k_ZeB1(1:443,:);
%unstable_k_ZeB1 = k_ZeB1(443:546,:);
unstable_k_ZeB1 = k_ZeB1(443:518,:);
stable_lower_k_ZeB1 = k_ZeB1(518:900,:);
% %----------------------------------------
% %extraction of slug dataset
% stable_upper_SLUG = SLUG(1:443,:);
% unstable_SLUG = SLUG(443:546,:);
% stable_lower_SLUG = SLUG(546:900,:);
% %-------------------------------------------------
% %extraction of mir200 dataset
% stable_upper_miR200 = miR200(1:348,:);
% unstable_miR200 = miR200(348:546,:);
% stable_lower_miR200 = miR200(546:900,:);
%-------------------
%extraction of mir200 dataset
% stable_upper_CDH1 = CDH1(1:443,:);
% unstable_CDH1 = CDH1(443:546,:);
% stable_lower_CDH1 = CDH1(546:900,:);
%--------------------
% %extraction of zeb dataset
% %stable_upper
% stable_upper_ZeB1 = ZeB1(1:348,:);
% %unstable_middle
% unstable_ZeB1 = ZeB1(348:546,:);
% %stable_lower
% stable_lower_ZeB1 = ZeB1(546:900,:);
%===================================
% %ploting
% figure;
% plot(k_ZeB1,miR200,'k','Linewidth',3);
% xlim([-0.01 1])
% ylim([-5 80])
% %%%%%%%%%%%%%%%%%%
% figure;
% plot(k_ZeB1,ZeB1,'k','Linewidth',3);
% xlim([-0.01 1])
% ylim([-5 140])
% % figure;
% plot(stable_upper_k_ZeB1,stable_upper_ZeB1,'k','Linewidth',3);
% hold on
% plot(unstable_k_ZeB1,unstable_ZeB1,'color','#FF0A37','Linestyle','-.','Linewidth',3);
% hold on
% plot(stable_lower_k_ZeB1,stable_lower_ZeB1,'k','Linewidth',3);
% xlabel('k ZeB1');
% ylabel('ZeB1');
% xlim([-0.01 1])
% ylim([-5 140])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %EMTscore_calculation
% EMTscore = ((ZeB1+SLUG)-(miR200+CDH1));
% stable_upper_EMTscore = EMTscore(1:348,:);
% unstable_EMTscore = EMTscore(348:546,:);
% stable_lower_EMTscore = EMTscore(546:900,:);
% %-------------------------------
% figure%to plot emtscore bifurcation
% %plot
% plot(stable_upper_k_ZeB1,stable_upper_EMTscore,'k','Linewidth',3);
% hold on
% plot(unstable_k_ZeB1,unstable_EMTscore,'color','#FF0A37','Linestyle','-.','Linewidth',3);
% hold on
% plot(stable_lower_k_ZeB1,stable_lower_EMTscore,'k','Linewidth',3);
% xlabel('k ZeB1');
% ylabel('EMTscore');
% ax=gca;
% %set axis color to dark black
% ax.YColor = 'k';
% ax.XColor = 'k';
% % Set x and y font sizes.
% %%ax.XAxis.FontSize = 20;
% %%ax.YAxis.FontSize = 20;
% % The below would set everything: title, x axis, y axis, and tick mark label font sizes.
% ax.FontSize = 24;
% % Bold all labels.
% ax.FontWeight = 'bold';
% xlim([-0.01 1])
% print('EMTscore_bifur_parase4562','-dpng','-r300')%300 is the resolution in dpi
% print('EMTscore_bifur_parase4562ifur','-depsc','-tiff','-r300')%300 is the resolution in dpi
% savefig('EMTscore_bifur_parase4562.fig')
% %===================================================
% figure%to plot mirgene bifurcation
% %plot
% plot(stable_upper_k_ZeB1,stable_upper_miR200,'k','Linewidth',3);
% hold on
% plot(unstable_k_ZeB1,unstable_miR200,'color','#FF0A37','Linestyle','-.','Linewidth',3);
% hold on
% plot(stable_lower_k_ZeB1,stable_lower_miR200,'k','Linewidth',3);
% xlabel('k ZeB1');
% ylabel('miR200');
% ax=gca;
% %set axis color to dark black
% ax.YColor = 'k';
% ax.XColor = 'k';
% % Set x and y font sizes.
% %%ax.XAxis.FontSize = 20;
% %%ax.YAxis.FontSize = 20;
% % The below would set everything: title, x axis, y axis, and tick mark label font sizes.
% ax.FontSize = 24;
% % Bold all labels.
% ax.FontWeight = 'bold';
% xlim([-0.01 1])
% %ylim([0 2])
% print('mir_bifur_parase4562','-dpng','-r300')%300 is the resolution in dpi
% print('mir_bifur_parase4562ifur','-depsc','-tiff','-r300')%300 is the resolution in dpi
% savefig('mir_bifur_parase4562.fig')
% %===============================
% 
% %------------------------------------------------------
% %plot zeb gene bifurcation
% figure;
% plot(stable_upper_k_ZeB1,stable_upper_ZeB1,'k','Linewidth',3);
% hold on
% plot(unstable_k_ZeB1,unstable_ZeB1,'color','#FF0A37','Linestyle','-.','Linewidth',3);
% hold on
% plot(stable_lower_k_ZeB1,stable_lower_ZeB1,'k','Linewidth',3);
% xlabel('k ZeB1');
% ylabel('ZeB1');
% ax=gca;
% %set axis color to dark black
% ax.YColor = 'k';
% ax.XColor = 'k';
% % Set x and y font sizes.
% %%ax.XAxis.FontSize = 20;
% %%ax.YAxis.FontSize = 20;
% % The below would set everything: title, x axis, y axis, and tick mark label font sizes.
% ax.FontSize = 24;
% % Bold all labels.
% ax.FontWeight = 'bold';
% xlim([-0.01 1])
% %ylim([-0.1 3])
% %savefig('bifurcation_paraset4562.fig')
% print('zeb_bifur_parase4562','-dpng','-r300')%300 is the resolution in dpi
% print('zeb_bifur_parase4562ifur','-depsc','-tiff','-r300')%300 is the resolution in dpi
% savefig('zeb_bifur_parase4562.fig')
% %-------------------------------------------
% %plot slug gene bifurcation
% figure;
% plot(stable_upper_k_ZeB1,stable_upper_SLUG,'k','Linewidth',3);
% hold on
% plot(unstable_k_ZeB1,unstable_SLUG,'color','#FF0A37','Linestyle','-.','Linewidth',3);
% hold on
% plot(stable_lower_k_ZeB1,stable_lower_SLUG,'k','Linewidth',3);
% xlabel('k ZeB1');
% ylabel('SLUG');
% ax=gca;
% %set axis color to dark black
% ax.YColor = 'k';
% ax.XColor = 'k';
% % Set x and y font sizes.
% %%ax.XAxis.FontSize = 20;
% %%ax.YAxis.FontSize = 20;
% % The below would set everything: title, x axis, y axis, and tick mark label font sizes.
% ax.FontSize = 24;
% % Bold all labels.
% ax.FontWeight = 'bold';
% xlim([-0.01 1])
% print('slug_bifur_parase4562','-dpng','-r300')%300 is the resolution in dpi
% print('slug_bifur_parase4562ifur','-depsc','-tiff','-r300')%300 is the resolution in dpi
% savefig('slug_bifur_parase4562.fig')
% %---------------------------------------------
%plot CDH1 gene bifurcation
% figure;
% plot(stable_upper_k_ZeB1,stable_upper_CDH1,'k','Linewidth',3);
% hold on
% plot(unstable_k_ZeB1,unstable_CDH1,'color','#FF0A37','Linestyle','-.','Linewidth',3);
% hold on
% plot(stable_lower_k_ZeB1,stable_lower_CDH1,'k','Linewidth',3);
% xlabel('k ZeB1');
% ylabel('CDH1');
% ax=gca;
% %set axis color to dark black
% ax.YColor = 'k';
% ax.XColor = 'k';
% % Set x and y font sizes.
% %%ax.XAxis.FontSize = 20;
% %%ax.YAxis.FontSize = 20;
% % The below would set everything: title, x axis, y axis, and tick mark label font sizes.
% ax.FontSize = 24;
% % Bold all labels.
% ax.FontWeight = 'bold';
% xlim([-0.01 1])
% print('CDH1_bifur_parase4562','-dpng','-r300')%300 is the resolution in dpi
% print('CDH1_bifur_parase4562ifur','-depsc','-tiff','-r300')%300 is the resolution in dpi
% savefig('CDH1_bifur_parase4562.fig')
%-----------------------------------------------------
%plot PDL1 gene bifurcation
figure;
plot(stable_upper_k_ZeB1,stable_upper_PDL1,'k','Linewidth',3);
hold on
plot(unstable_k_ZeB1,unstable_PDL1,'color','#FF0A37','Linestyle','-.','Linewidth',3);
hold on
plot(stable_lower_k_ZeB1,stable_lower_PDL1,'k','Linewidth',3);
xlabel('k ZeB1');
ylabel('PDL1');
ax=gca;
%set axis color to dark black
ax.YColor = 'k';
ax.XColor = 'k';
% Set x and y font sizes.
%%ax.XAxis.FontSize = 20;
%%ax.YAxis.FontSize = 20;
% The below would set everything: title, x axis, y axis, and tick mark label font sizes.
ax.FontSize = 24;
% Bold all labels.
ax.FontWeight = 'bold';
xlim([-0.01 0.5])
print('PDL1_bifur_parase4562','-dpng','-r300')%300 is the resolution in dpi
print('PDL1_bifur_parase4562ifur','-depsc','-tiff','-r300')%300 is the resolution in dpi
savefig('PDL1_bifur_parase4562.fig')