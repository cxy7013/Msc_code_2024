function plot_beam_D_C(x_vec, y_vec, I, figure_title)
    % Calculate radius and centroid
    [radius, centroid] = calculate_d_c(I, x_vec, y_vec);
    
    % Plotting the beam image
    imagesc(x_vec, y_vec, I);
    hold on; axis xy;
    colormap('hot'); colorbar;
    xlabel('x (m)'); ylabel('y (m)');
    title(figure_title, 'FontWeight', 'bold');
    axis square;


    % Plotting the D86 circle
    viscircles(centroid, radius, 'Color', 'b', 'LineWidth', 0.5);
    
    % mark centroid
    plot(centroid(1), centroid(2), 'bx', 'MarkerSize', 10, 'LineWidth', 0.5);

    % display value
    text(min(x_vec), max(y_vec), sprintf('R_{D86} = %.2f mm', radius * 1e3), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Color', 'white', 'FontSize', 10);
    text(min(x_vec), min(y_vec), sprintf('Centroid = (%.2f, %.2f) mm', centroid(1) * 1e3, centroid(2) * 1e3), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'white', 'FontSize', 10);
    hold off;
end
