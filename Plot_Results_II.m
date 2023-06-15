lw = 1.5; %line width
pw = 4;   %point size

% Time vector
t = ti:ts:tf;

%% ///////// Interaction force /////////
% Compute Forward Kinematics from Simulation results
[BM_FK.x, BM_FK.y, BM_FK.z, BM_FK.phi, BM_FK.theta, BM_FK.psi, BM_FK.RL] = ForwardKinematics(model,y_BM(1:2,:)');

figure(9); 
subplot(3,1,1);
% title(['External torques'],'Interpreter','latex','FontSize',12);
xlabel('time (s)','Interpreter','latex','FontSize',12);
ylabel(['External' newline 'Torque (Nm)']'','Interpreter','latex','FontSize',12); 
% set(get(gca,'ylabel'),'rotation',0)
hold on; 
plot(t,y_BM(77,:),'LineWidth',lw,'Color',[0, 0.4470, 0.7410]);
plot(t,y_BM(78,:),'LineWidth',lw,'Color',[0.8500, 0.3250, 0.0980]);
plot(t,y_BM(79,:),":",'LineWidth',lw,'Color',[0.9290, 0.6940, 0.1250]);
hold off;
grid on;
leg1 = legend('$n_{x}$','$n_{y}$','$n_{z}$'); 
set(leg1,'Interpreter','latex','FontSize',11,'Orientation','horizontal'); 
xlim([0 5]);
ylim([-0.2 5.5]);

subplot(3,1,2);
% title(['External forces'],'Interpreter','latex','FontSize',12);
xlabel('time (s)','Interpreter','latex','FontSize',12);
ylabel(['External' newline 'Force (N)'],'Interpreter','latex','FontSize',12); 
% set(get(gca,'ylabel'),'rotation',0)
hold on; 
plot(t,y_BM(80,:),'LineWidth',lw,'Color',[0.4940, 0.1840, 0.5560]);
plot(t,y_BM(81,:),'LineWidth',lw,'Color',[0.3010, 0.7450, 0.9330]);
plot(t,y_BM(82,:),'LineWidth',lw,'Color',[0.6350, 0.0780, 0.1840]);
hold off;
grid on;
leg1 = legend('$f_{x}$','$f_{y}$','$f_{z}$'); 
set(leg1,'Interpreter','latex','FontSize',11,'Orientation','horizontal'); 
xlim([0 5]);
ylim([-0.2 4.5]);

subplot(3,1,3);
% title('Wrist position','Interpreter','latex','FontSize',14);
xlabel('time (s)','Interpreter','latex','FontSize',12);
ylabel(['Wrist' newline 'position (m)'],'Interpreter','latex','FontSize',12); 
% set(get(gca,'ylabel'),'rotation',0)
xlim([0 5]);
ylim([-0.25 1.1]);
hold on; 
plot(t,BM_FK.x(3,:),'LineWidth',lw,'Color',[0, 0.4470, 0.7410]);
plot(t,BM_FK.y(3,:),'LineWidth',lw,'Color',[0.8500, 0.3250, 0.0980]); 
plot(t,BM_FK.z(3,:),'LineWidth',lw,'Color',[0.4660, 0.6740, 0.1880]);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('$x_{w}$','$y_{w}$','$z_{w}$'); 
set(leg1,'Interpreter','latex','FontSize',11,'Orientation','horizontal'); 