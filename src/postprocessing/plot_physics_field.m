% This plotting code was created to visualize physics evolution of a
% shear layer where the across-layer domain was resolved.
% Author: Yu-Han W. 
% Date: 25 01 2025

%% Determine the folder to save outputs
% outputFolder is defined in the modelling script
outputFolderPath = fullfile(parentFolder, 'post', outputFolder, ['Rupture',num2str(i_t/freq_post),'.png']);

%% Plotting panel setup
set(fig, 'Position', [0, 0, 1000, 1200], 'Units', 'pixels', 'Color', 'w'); % Set figure background to white
clf
ax = [];
% Slip rate versus fault length
ax(1) = subplot(4,3,1);
plot(xx/L, V, 'LineWidth', 1.5); 

xlim([-1/2, 1/2]);
if inject_switch
    set(ax(1), 'YScale', 'log')  % Apply log scale to just this subplot
    ylim([1e-9, 20]);
else
    if max(V) < 5
        ylim([min(V), 5]);
    else
        ylim([min(V), max(V)*1.025]);
    end
end

hold on;
xleft = xx(nx_left) / L;
xright = xx(nx_right) / L;
yl = ylim;
plot([xleft, xleft], yl, 'k--', 'LineWidth', 1);
plot([xright, xright], yl, 'k--', 'LineWidth', 1);
hold off;
% axis tight;
title(['$V$ along the fault'], 'Interpreter', 'latex', 'FontSize',15);
xlabel('$x/L$','Interpreter', 'latex', 'FontSize',15);
ylabel('Slip rate $V$ [m/s]','Interpreter', 'latex', 'FontSize',15);
pbaspect([45 20 1]); % Maintain the aspect ratio as 16:9
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15)

% Boundary pressure
ax(2) = subplot(4,3,2);
plot(xx/L, p(:,1)/1e6, 'LineWidth', 1.5, 'Color', 'red'); 
hold on;
plot(xx/L, p(:,end)/1e6, 'LineWidth', 1.5, 'Color', 'blue'); 
hold off;
xlim([-1/2, 1/2]);

% axis tight;
title(['$p$ at the boundary'], 'Interpreter', 'latex', 'FontSize',15);
xlabel('$x/L$','Interpreter', 'latex', 'FontSize',15);
ylabel('Pore pressure [MPa]','Interpreter', 'latex', 'FontSize',15);
pbaspect([45 20 1]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15)

% Boundary temperature
ax(3) = subplot(4,3,3);
plot(xx/L, T(:,1), 'LineWidth', 1.5, 'Color', 'red', 'LineStyle', '--'); 
hold on;
plot(xx/L, T(:,end), 'LineWidth', 1.5, 'Color', 'blue', 'LineStyle', '-.'); 
hold off;
xlim([-1/2, 1/2]);

% axis tight;
title(['$T$ at the boundary'], 'Interpreter', 'latex', 'FontSize',15);
xlabel('$x/L$','Interpreter', 'latex', 'FontSize',15);
ylabel('Temperature [K]','Interpreter', 'latex', 'FontSize',15);
pbaspect([45 20 1]); % Maintain the aspect ratio as 16:9
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15)

% Profile of shear strain rate, gamma
ax(4) = subplot(4,3,4);
evolution_plot(real(log(gamma'/gamma_c + exp(-5))), 15, ...
    yy/h, xx/L, '$x / l_f$', '$y / h$', '$\mathrm{log} \dot{\gamma}$', ...
    [-5,5], 'viridis_white', 'Shear strain rate $\mathrm{log} \dot{\gamma} / \dot{\gamma}_c$',[45 20 1]);

% Profile of pore pressure p
ax(5) = subplot(4,3,5);
evolution_plot(p'/1e6, 15, ...
    yy/h, xx/L, '$x / l_f$', '$y / h$', '$p$ [MPa]', ...
    [min(p(:))/1e6-0.5, max(p(:)/1e6)], 'viridis_white', 'Pore pressure $p$',[45 20 1]);

% Profile of temperature T
ax(6) = subplot(4,3,6);
evolution_plot(T', 15, ...
    yy/h, xx/L, '$x / l_f$', '$y / h$', '$T$ [K]', ...
    [0, max(T(:))], 'viridis_white', 'Temperature $T$',[45 20 1]);

% Sections of shear strain rate, gamma
ax(7) = subplot(4,3,7);
plot(real(log(gamma(nx/2, :)/gamma_c+exp(-5))), flip(yy_gammadot/h), 'LineWidth', 1.5); 
hold on;
plot(real(log(gamma(nx/2+round(nx/10), :)/gamma_c+ exp(-5))), flip(yy_gammadot/h), 'LineWidth', 1.5); 
plot(real(log(gamma(nx/2-round(nx/10), :)/gamma_c+ exp(-5))), flip(yy_gammadot/h), 'LineWidth', 1.5); 
xlim([-5, 5]);
ylim([-1/2, 1/2]);
% axis tight;
title(['$\dot{\gamma}/ \dot{\gamma}_c$ slice across gouge'], 'Interpreter', 'latex', 'FontSize',15);
xlabel('$\dot{\gamma}/ \dot{\gamma}_c$','Interpreter', 'latex', 'FontSize',15);
ylabel('$y/h$','Interpreter', 'latex', 'FontSize',15);
pbaspect([45 20 1]); % Maintain the aspect ratio as 16:9
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15)
hold off;

% Sections of p
ax(8) = subplot(4,3,8);
plot(p(nx/2, :)/1e6, flip(yy/h), 'LineWidth', 1.5); 
hold on;
plot(p(nx/2+round(nx/10), :)/1e6, flip(yy/h), 'LineWidth', 1.5); 
plot(p(nx/2-round(nx/10), :)/1e6, flip(yy/h), 'LineWidth', 1.5); 
xlim([min(p(:))/1e6-0.1, max(p(:)/1e6+0.1)]);
ylim([-1/2, 1/2]);
% axis tight;
title(['$p$ slice across gouge'], 'Interpreter', 'latex', 'FontSize',15);
xlabel('$p$ [MPa]','Interpreter', 'latex', 'FontSize',15)
ylabel('$y/h$','Interpreter', 'latex', 'FontSize',15)
pbaspect([45 20 1]); % Maintain the aspect ratio as 16:9
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15)
hold off;

% Sections of T
ax(9) = subplot(4,3,9);
plot(flip(T(nx/2, :)), yy/h, 'LineWidth', 1.5); 
hold on;
plot(flip(T(nx/2+round(nx/10), :)), yy/h, 'LineWidth', 1.5); 
plot(flip(T(nx/2-round(nx/10), :)), yy/h, 'LineWidth', 1.5); 
xlim([0, max(T(:))+0.5]);
ylim([-1/2, 1/2]);
% axis tight;
title(['$T$ slice across gouge'], 'Interpreter', 'latex', 'FontSize',15);
xlabel('$T$','Interpreter', 'latex', 'FontSize',15)
ylabel('$y/h$','Interpreter', 'latex', 'FontSize',15)
pbaspect([45 20 1]); % Maintain the aspect ratio as 16:9
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15)
hold off;

% Shear stress
ax(10) = subplot(4,3,10);
plot(xx/L, tauf/1e6, 'LineWidth', 1.5, 'Color', 'red'); 
hold on;
plot(xx/L, tauf_/1e6, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '--');
if max((tauf - tauf_)./tauf_) > 1e-3
    warning("resultant shear stress is not balanced by friction")
end
xlim([-1/2, 1/2]);
ylim([0, max(max(tauf))/1e6+0.1]);
% axis tight;
title(['Shear stress $\tau$'], 'Interpreter', 'latex', 'FontSize',15);
xlabel('$x/L$','Interpreter', 'latex', 'FontSize',15);
ylabel('$\tau$ [MPa]','Interpreter', 'latex', 'FontSize',15);
pbaspect([45 20 1]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15)

% Effective stress
ax(11) = subplot(4,3,11);
evolution_plot((sigma_normal)'/1e6, 15, ...
    yy/h, xx/L, '$x / l_f$', '$y / h$', '$f$', ...
    [min((sigma_normal),[],'all'), max((sigma_normal),[],'all')]./1e6, ...
    'viridis_white', 'Effective stress $\sigma_{eff}$',[45 20 1]);

% friction
ax(12) = subplot(4,3,12);
friction = (tau_shear./sigma_normal)';
evolution_plot(friction, 15, ...
    yy/h, xx/L, '$x / l_f$', '$y / h$', '$f$', ...
    [min(friction,[],'all'), max(friction,[],'all')], ...
    'viridis_white', 'Friction $f$',[45 20 1]);

% Adjust margin between figures
pos=[];
for i=1:1:12
    pos(i,:) = get(ax(i), 'Position');
end

% Reduce space between the two subplots
for i=1:3
    set(ax(i), 'Position', [pos(i,1), pos(i,2)-0.005, pos(i,3), pos(i,4)]);
  
    set(ax(i+3), 'Position', [pos(i+3,1), pos(i+3,2)+0.005, pos(i+6,3), pos(i+6,4)]);
  
    set(ax(i+6), 'Position', [pos(i+6,1), pos(i+6,2)+0.1, pos(i+6,3), pos(i+6,4)]);

    if i+9 ~= 10
      set(ax(i+9), 'Position', [pos(i+9,1), pos(i+9,2)+0.14, pos(i+9,3), pos(i+9,4)]);
    else
      set(ax(i+9), 'Position', [pos(i+9,1), pos(i+9,2)+0.175, pos(i+9,3), pos(i+9,4)]);
    end
end 

sgtitle(['Scaled time $ ', 't_{\dot{\gamma}_c}=', sprintf('%.4f',timev(end)*vs*fr0*sigma_zz0/h/G), ...
    ', \delta_{x} / \delta_{c}=', sprintf('%.4f',max(delta_x./delta_xc)), '$'], ...
     'Interpreter', 'latex', 'FontSize',15);

print(gcf,outputFolderPath,'-dpng','-r600');
