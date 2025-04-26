function vivid = vivid_colormap()
    % Define a custom vivid colormap
    vivid = [
        1, 0, 0;      % Red
        1, 0.5, 0;    % Orange
        1, 1, 0;      % Yellow
        0, 1, 0;      % Green
        0, 1, 1;      % Cyan
        0, 0, 1;      % Blue
        0.5, 0, 1;    % Purple
        1, 0, 1       % Magenta
    ];
    
    % Interpolate to make it smoother
    n = 256; % Number of colors in the colormap
    x = linspace(1, size(vivid, 1), n);
    vivid = interp1(1:size(vivid, 1), vivid, x);
end