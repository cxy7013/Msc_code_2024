function plot_beam_D_C_big(x_vec, y_vec, I, figure_title, enlarge_factor)
    % enlarge_factor 
    if nargin < 5
        enlarge_factor = 2; 
    end
    
    [radius, centroid] = calculate_d_c(I, x_vec, y_vec);
    
    imagesc(x_vec, y_vec, I);
    hold on; axis xy;
    colormap('hot'); colorbar;
    xlabel('x (m)'); ylabel('y (m)');
    title(figure_title, 'FontWeight', 'bold');

    viscircles(centroid, radius, 'Color', 'b', 'LineWidth', 0.5);
    plot(centroid(1), centroid(2), 'bx', 'MarkerSize', 10, 'LineWidth', 0.5);

    x_range = (max(x_vec) - min(x_vec)) / enlarge_factor;
    y_range = (max(y_vec) - min(y_vec)) / enlarge_factor;
    
    xlim([centroid(1) - x_range/2, centroid(1) + x_range/2]);
    ylim([centroid(2) - y_range/2, centroid(2) + y_range/2]);

    text(centroid(1) - x_range/2, centroid(2) + y_range/2, sprintf('R_{D86} = %.2f mm', radius * 1e3), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Color', 'white', 'FontSize', 10);
    text(centroid(1) - x_range/2, centroid(2) - y_range/2, sprintf('Centroid = (%.2f, %.2f) mm', centroid(1) * 1e3, centroid(2) * 1e3), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', 'white', 'FontSize', 10);
    
    hold off;
end
