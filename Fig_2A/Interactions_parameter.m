function dy = PDL1_eqn(t,y)
dy = zeros(5,1);    % a column vector
%parameter value
g_ZEB1 = 51.205355;
g_miR200 = 73.744958;
g_PDL1 = 43.968037;
g_CDH1 = 89.201165;
g_SLUG = 29.233455;


k_ZEB1 = 0.62873;
k_miR200 = 0.861129;
k_PDL1 = 0.837638;
k_CDH1 = 0.192489;
k_SLUG = 0.549676;


x_ZEB1_miR200 = 1.573366;
x_miR200_ZEB1 = 6.155112;
x_ZEB1_ZEB1 = 0.662827;
x_miR200_PDL1 = 9.788561;
x_PDL1_CDH1 = 34.068816;
x_ZEB1_CDH1 = 0.80368;
x_CDH1_ZEB1 = 0.547427;
x_SLUG_miR200 = 5.575322;
x_miR200_SLUG = 8.200188;
x_SLUG_ZEB1 = 10.863898;
x_SLUG_SLUG = 13.8752;
x_SLUG_CDH1 = 11.379845;


l_ZEB1_miR200 = 0.018446;
l_miR200_ZEB1 = 0.014161;
l_ZEB1_ZEB1 = 54.054369;
l_miR200_PDL1 = 0.01356;
l_PDL1_CDH1 = 0.010248;
l_ZEB1_CDH1 = 0.041006;
l_CDH1_ZEB1 = 0.01252;
l_SLUG_miR200 = 0.01051;
l_miR200_SLUG = 0.16813;
l_SLUG_ZEB1 = 72.009524;
l_SLUG_SLUG = 25.647124;
l_SLUG_CDH1 = 0.06103;


n_ZEB1_miR200 = 1;
n_miR200_ZEB1 = 2;
n_ZEB1_ZEB1 = 2;
n_miR200_PDL1 = 6;
n_PDL1_CDH1 = 2;
n_ZEB1_CDH1 = 3;
n_CDH1_ZEB1 = 2;
n_SLUG_miR200 = 6;
n_miR200_SLUG = 3;
n_SLUG_ZEB1 = 5;
n_SLUG_SLUG = 3;
n_SLUG_CDH1 = 6;


%ode eqn                  
dy(1) = g_ZEB1*hill(y(2),x_miR200_ZEB1,l_miR200_ZEB1,n_miR200_ZEB1)*hilla(y(1),x_ZEB1_ZEB1,l_ZEB1_ZEB1,n_ZEB1_ZEB1)*hill(y(4),x_CDH1_ZEB1,l_CDH1_ZEB1,n_CDH1_ZEB1)*hilla(y(5),x_SLUG_ZEB1,l_SLUG_ZEB1,n_SLUG_ZEB1) - k_ZEB1*y(1);
dy(2) = g_miR200*hill(y(1),x_ZEB1_miR200,l_ZEB1_miR200,n_ZEB1_miR200)*hill(y(5),x_SLUG_miR200,l_SLUG_miR200,n_SLUG_miR200) - k_miR200*y(2);
dy(3) = g_PDL1*hill(y(2),x_miR200_PDL1,l_miR200_PDL1,n_miR200_PDL1) - k_PDL1*y(3);
dy(4) = g_CDH1*hill(y(3),x_PDL1_CDH1,l_PDL1_CDH1,n_PDL1_CDH1)*hill(y(1),x_ZEB1_CDH1,l_ZEB1_CDH1,n_ZEB1_CDH1)*hill(y(5),x_SLUG_CDH1,l_SLUG_CDH1,n_SLUG_CDH1) - k_CDH1*y(4);
dy(5) = g_SLUG*hill(y(2),x_miR200_SLUG,l_miR200_SLUG,n_miR200_SLUG)*hilla(y(5),x_SLUG_SLUG,l_SLUG_SLUG,n_SLUG_SLUG) - k_SLUG*y(5);

