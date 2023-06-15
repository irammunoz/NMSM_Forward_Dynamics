lw = 1.5; %line width
pw = 4;   %point size

% Time vector
t = ti:ts:tf;

%% ///////// Exitation and Activation signals ///////// 
figure(1); 
subplot(2,1,1); title('\textbf{Muscle excitation}','Interpreter','latex','FontSize',14);
xlabel('time (s)','Interpreter','latex','FontSize',14);
xlim([ti tf]);
ylim([0 1.2]);
hold on; 
plot(t,y_BM(86,:),'LineWidth',1.5);
plot(t,y_BM(87,:),'LineWidth',1.5);
hold off;
leg1 = legend('brachialis','triceps medial head','Orientation','horizontal'); 
set(leg1,'Interpreter','latex','FontSize',10);grid on;

subplot(2,1,2); title('\textbf{Muscle activation}','Interpreter','latex','FontSize',14);
xlabel('time (s)','Interpreter','latex','FontSize',14);
xlim([ti tf]);
ylim([0 1.2]);
hold on; 
plot(t,y_BM(5,:),'LineWidth',1.5);
plot(t,y_BM(6,:),'LineWidth',1.5);
hold off;
leg1 = legend('brachialis','triceps medial head','Orientation','horizontal'); 
set(leg1,'Interpreter','latex','FontSize',10);grid on;

%% ///////// Muscle length and velocity ///////// 
figure(2);
subplot(2,1,1);
xlabel('time (s)','Interpreter','latex','FontSize',14);
ylabel(['Muscle-tendon' newline 'length (m)'],'Interpreter','latex','FontSize',14); 
xlim([ti tf]);
ylim([0.1 0.24]);
hold on;
plot(t,y_BM(61,:),'LineWidth',lw,'Color',[0.00,0.45,0.74]);
plot(t,y_BM(62,:),'LineWidth',lw,'Color',[0.85,0.33,0.10]);
plot(states.data(:,1),Length.data(:,3),'-.','LineWidth',lw,'Color',[0.00,0.45,0.74]); 
plot(states.data(:,1),Length.data(:,2),'-.','LineWidth',lw,'Color',[0.85,0.33,0.10]);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('brachialis','triceps medial head','Orientation','horizontal'); 
set(leg1,'Interpreter','latex','FontSize',10); grid on;

subplot(2,1,2);
xlabel('time (s)','Interpreter','latex','FontSize',14);
ylabel(['Muscle-tendon' newline 'velocity (m/s)'],'Interpreter','latex','FontSize',14); 
xlim([ti tf]);
ylim([-0.15 0.15]);
hold on; 
plot(t,y_BM(63,:),'LineWidth',lw,'Color',[0.00,0.45,0.74]); 
plot(t,y_BM(64,:),'LineWidth',lw,'Color',[0.85,0.33,0.10]); 
plot(states.data(:,1),FiberVelocity.data(:,3),'-.','LineWidth',lw,'Color',[0.00,0.45,0.74]); 
plot(states.data(:,1),FiberVelocity.data(:,2),'-.','LineWidth',lw,'Color',[0.85,0.33,0.10]);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('brachialis','triceps medial head','Orientation','horizontal'); 
set(leg1,'Interpreter','latex','FontSize',10); grid on;

%% ///////// Joint torques /////////
figure(3);
subplot(2,2,1);
title(['Brachialis' newline 'contribution'],'Interpreter','latex','FontSize',12);
xlabel('time (s)','Interpreter','latex','FontSize',12);
ylabel('Joint torques (Nm)','Interpreter','latex','FontSize',12); 
xlim([ti tf]);
ylim([0 14]);
hold on; 
plot(t,y_BM(71,:),'LineWidth',lw,'Color',[0.9290, 0.6940, 0.1250]);
plot(t,y_BM(72,:),'LineWidth',lw,'Color',[0.4660, 0.6740, 0.1880]); 
plot(states.data(:,1),Moment_r_shoulder_elev.data(:,3),'-.','LineWidth',lw,'Color',[0.9290, 0.6940, 0.1250]);
plot(states.data(:,1),Moment_r_elbow_flex.data(:,3),'-.','LineWidth',lw,'Color',[0.4660, 0.6740, 0.1880]);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('$\tau_{s}^{bra}$','$\tau_{e}^{bra}$');
set(leg1,'Interpreter','latex','FontSize',11);

subplot(2,2,2);
title(['Triceps medial head' newline 'contribution'],'Interpreter','latex','FontSize',12);
xlabel('time (s)','Interpreter','latex','FontSize',12);
ylabel('Joint torques (Nm)','Interpreter','latex','FontSize',12); 
xlim([ti tf]);
ylim([-2.0 1.4]);
hold on; 
plot(t,y_BM(73,:),'LineWidth',lw,'Color',[0.9290, 0.6940, 0.1250]);
plot(t,y_BM(74,:),'LineWidth',lw,'Color',[0.4660, 0.6740, 0.1880]); 
plot(states.data(:,1),Moment_r_shoulder_elev.data(:,2),'-.','LineWidth',lw,'Color',[0.9290, 0.6940, 0.1250]); 
plot(states.data(:,1),Moment_r_elbow_flex.data(:,2),'-.','LineWidth',lw,'Color',[0.4660, 0.6740, 0.1880]);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('$\tau_{s}^{tmh}$','$\tau_{e}^{tmh}$');
set(leg1,'Interpreter','latex','FontSize',11);

subplot(2,2,[3,4]);
xlabel('time (s)','Interpreter','latex','FontSize',12);
ylabel('Joint torques (Nm)','Interpreter','latex','FontSize',12); 
xlim([ti tf]);
ylim([-0.2 11]);
hold on; 
plot(t,y_BM(71,:)+y_BM(73,:),'LineWidth',lw,'Color',[0.9290, 0.6940, 0.1250]);
plot(t,y_BM(72,:)+y_BM(74,:),'LineWidth',lw,'Color',[0.4660, 0.6740, 0.1880]); 
plot(states.data(:,1),Moment_r_shoulder_elev.data(:,2)+Moment_r_shoulder_elev.data(:,3),'-.','LineWidth',lw,'Color',[0.9290, 0.6940, 0.1250]); 
plot(states.data(:,1),Moment_r_elbow_flex.data(:,2)+Moment_r_elbow_flex.data(:,3),'-.','LineWidth',lw,'Color',[0.4660, 0.6740, 0.1880]); 
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('$\tau_{s}$','$\tau_{e}$');
set(leg1,'Interpreter','latex','FontSize',11);

%% ///////// Generalized positions and velocities /////////
figure(4);
subplot(2,1,1);
xlabel('time (s)','Interpreter','latex','FontSize',14);
ylabel(['Generalized' newline 'positions (rad)'],'Interpreter','latex','FontSize',14); 
xlim([ti tf]);
ylim([-0.5 4.0]);
hold on; 
plot(t,y_BM(1,:),'LineWidth',lw,'Color',[1, 0, 0]);
plot(t,y_BM(2,:),'LineWidth',lw,'Color',[0, 0, 0]); 
plot(states.data(:,1),states.data(:,2),'-.','LineWidth',lw,'Color',[1, 0, 0]); 
plot(states.data(:,1),states.data(:,4),'-.','LineWidth',lw,'Color',[0, 0, 0]);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('$q_{s}$','$q_{e}$');
set(leg1,'Interpreter','latex','FontSize',11);

subplot(2,1,2);
xlabel('time (s)','Interpreter','latex','FontSize',14);
ylabel(['Generalized' newline 'velocities (rad/s)'],'Interpreter','latex','FontSize',14); 
xlim([ti tf]);
% ylim([0 1.0]);
hold on; 
plot(t,y_BM(3,:),'LineWidth',lw,'Color',[1, 0, 0]);
plot(t,y_BM(4,:),'LineWidth',lw,'Color',[0, 0, 0]);
plot(states.data(:,1),states.data(:,3),'-.','LineWidth',lw,'Color',[1, 0, 0]); 
plot(states.data(:,1),states.data(:,5),'-.','LineWidth',lw,'Color',[0, 0, 0]);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('$\dot{q}_{s}$','$\dot{q}_{e}$');
set(leg1,'Interpreter','latex','FontSize',11);


%% ///////// Wrapping error /////////
axgrid = [2,2];  % [#rows, #cols]
titles = {'\textbf{Difference between SGW method and Obstacle-Set method}', '\textbf{Difference between SGW method and Muscle Analysis Tool (OpenSim)}'}; 
figure(5);
tclMain = tiledlayout(axgrid(1),1); 
tcl = gobjects(1,axgrid(1));
ax = gobjects(axgrid);

tcl(1) = tiledlayout(tclMain,1,axgrid(2));
tcl(1).Layout.Tile = 1; 
ax(1,1) = nexttile(tcl(1));
xlabel('time (s)','Interpreter','latex','FontSize',11);
ylabel(['Difference in' newline 'muscle-tendon' newline 'length (m)'],'Interpreter','latex','FontSize',11); 
xlim([ti tf]);
ylim([-1e-13 14e-13]);
hold on; 
plot(t,y_BM(61,:)-y_OS(61,:),'LineWidth',1.5);
plot(t,y_BM(62,:)-y_OS(62,:),'LineWidth',1.5);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('brachialis','triceps medial head','Orientation','vertical'); 
set(leg1,'Interpreter','latex','FontSize',9);

ax(1,2) = nexttile(tcl(1));
xlabel('time (s)','Interpreter','latex','FontSize',11);
ylabel(['Difference in' newline 'muscle-tendon' newline 'velocity (m/s)'],'Interpreter','latex','FontSize',11); 
xlim([ti tf]);
ylim([-9e-13 12e-13]);
hold on; 
plot(t,y_BM(63,:)-y_OS(63,:),'LineWidth',1.5);
plot(t,y_BM(64,:)-y_OS(64,:),'LineWidth',1.5);
shaded_patch_significant_timepoints(t, y_OS(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_OS(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('brachialis','triceps medial head','Orientation','vertical'); 
set(leg1,'Interpreter','latex','FontSize',9);

title(tcl(1),titles{1},'Interpreter','latex','FontSize',10)


tcl(2) = tiledlayout(tclMain,1,axgrid(2));
tcl(2).Layout.Tile = 2; 
ax(2,1) = nexttile(tcl(2));
xlabel('time (s)','Interpreter','latex','FontSize',11);
ylabel(['Difference in' newline 'muscle-tendon' newline 'length (m)'],'Interpreter','latex','FontSize',11);
xlim([ti tf]);
ylim([-2e-3 4e-3]);
hold on; 
plot(t,y_BM(61,:)-Length.data(:,3)','LineWidth',1.5);
plot(t,y_BM(62,:)-Length.data(:,2)','LineWidth',1.5);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('brachialis','triceps medial head','Orientation','vertical'); 
set(leg1,'Interpreter','latex','FontSize',9);

ax(2,2) = nexttile(tcl(2));
xlabel('time (s)','Interpreter','latex','FontSize',11);
ylabel(['Difference in' newline 'muscle-tendon' newline 'velocity (m/s)'],'Interpreter','latex','FontSize',11);
xlim([ti tf]);
ylim([-0.07 0.12]);
hold on; 
plot(t,y_BM(63,:)-FiberVelocity.data(:,3)','LineWidth',1.5);
plot(t,y_BM(64,:)-FiberVelocity.data(:,2)','LineWidth',1.5);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('brachialis','triceps medial head','Orientation','vertical'); 
set(leg1,'Interpreter','latex','FontSize',9);

title(tcl(2),titles{2},'Interpreter','latex','FontSize',10)

%% ///////// Forward Kinematics, Work and Power /////////
% Compute Forward Kinematics from Simulation results
[BM_FK.x, BM_FK.y, BM_FK.z, BM_FK.phi, BM_FK.theta, BM_FK.psi, BM_FK.RL] = ForwardKinematics(model,y_BM(1:2,:)');

figure(6);
subplot(2,2,[1,2]);
xlabel('time (s)','Interpreter','latex','FontSize',14);
ylabel('Wrist position (m)','Interpreter','latex','FontSize',14); 
xlim([ti tf]);
ylim([-0.25 1.0]);
hold on; 
plot(t,BM_FK.x(3,:),'LineWidth',lw,'Color',[0, 0.4470, 0.7410]);
plot(t,BM_FK.y(3,:),'LineWidth',lw,'Color',[0.8500, 0.3250, 0.0980]); 
plot(t,BM_FK.z(3,:),'LineWidth',lw,'Color',[0.4660, 0.6740, 0.1880]); 
plot(states.data(:,1),OS_FK.x(3,:),'--','LineWidth',lw,'Color',[0, 0.4470, 0.7410]);
plot(states.data(:,1),OS_FK.y(3,:),'--','LineWidth',lw,'Color',[0.8500, 0.3250, 0.0980]);
plot(states.data(:,1),OS_FK.z(3,:),'--','LineWidth',lw,'Color',[0.4660, 0.6740, 0.1880]);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;
leg1 = legend('$x_{w}$','$y_{w}$','$z_{w}$'); 
set(leg1,'Interpreter','latex','FontSize',11,'Orientation','horizontal'); 

subplot(2,2,3);
xlabel('time (s)','Interpreter','latex','FontSize',14);
ylabel('Power (W)','Interpreter','latex','FontSize',14); 
xlim([ti tf]);
ylim([-12 65]);
hold on; 
plot(t,y_BM(85,:),'LineWidth',lw,'Color',[1, 0, 0]);
plot(states.data(:,1),OS_pwr','--','LineWidth',lw,'Color',[1, 0, 0]);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;

subplot(2,2,4);
xlabel('time (s)','Interpreter','latex','FontSize',14);
ylabel('Work (J)','Interpreter','latex','FontSize',14); 
xlim([ti tf]);
ylim([0 21]);
hold on; 
plot(t,y_BM(9,:),'LineWidth',lw,'Color',[0, 0, 0]); 
plot(states.data(:,1),OS_wrk','--','LineWidth',lw,'Color',[0, 0, 0]);
shaded_patch_significant_timepoints(t, y_BM(51,:), [0.00,0.45,0.74]);
shaded_patch_significant_timepoints(t, y_BM(52,:), [0.85,0.33,0.10]);
hold off;
grid on;