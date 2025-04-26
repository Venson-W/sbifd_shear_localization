function result = fast_trapz(y, dy)
    % FAST_TRAPZ - High-performance trapezoidal integration for uniform grids
    %
    % Usage:
    %   result = fast_trapz(y, dy)       % integrate along 2nd dimension with spacing dy
    %   result = fast_trapz(y, dy, dim)  % integrate along dimension dim with spacing dy
    
    % Default dimension (along columns/2nd dimension)
    total = sum(y, 2);
    
    % Subtract half the endpoints (first and last columns)
    total = total - 0.5 * (y(:,1) + y(:,end));
    
    % Multiply by spacing
    result = total * dy;
end