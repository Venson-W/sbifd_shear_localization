function evolution_plot(physics_field, fontsize, x_coord, y_coord, x_label, y_label, z_label, z_caxis, colormap_name, figure_title, panel_size)
    % Function to visualize a 2D physics field
    % Arguments:
    % physics_field: A 2D matrix representing the physics field
    % fontsize: Font size for the axis labels and colorbar label
    % x_coord: Values for the x-axis (time step)
    % y_coord: Values for the y-axis (fault length)

    % Create the image plot with specified x and y coordinates
    imagesc(x_coord, y_coord, physics_field);
    
    % Add xlabel with LaTeX interpreter and specified font size
    xlabel(x_label, 'Interpreter', 'latex', 'FontSize', fontsize);
    
    % Add ylabel with LaTeX interpreter and specified font size
    ylabel(y_label, 'Interpreter', 'latex', 'FontSize', fontsize);
    
    % Add colorbar and label it with LaTeX interpreter and specified font size
    h = colorbar;
        % Set LaTeX interpreter and font size for colorbar ticks
    set(h, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize, 'Location', 'southoutside');

    ylabel(h, z_label, 'Interpreter', 'latex', 'FontSize', fontsize);
    % The scale range
    clim(z_caxis);

    % Set the colormap (optional, for better visualization)
    colormap(feval(colormap_name)); % You can choose other colormaps like 'parula', 'hot', etc.
    
    % Apply LaTeX formatting to the tick labels for both axes
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize, 'Color', [0 0 0] );
    %axis equal;
    axis tight;

    pbaspect(panel_size); % Maintain the aspect ratio as 16:9

    % Add title if needed (optional)
    title(figure_title, 'Interpreter', 'latex', 'FontSize', fontsize);
end
