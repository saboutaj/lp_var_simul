function plot_loss(horzs, results, add_line, plot_name, plot_legend)

    % Function for plotting results across estimation methods

    line_colors = repmat(lines(7),2,1);
    line_styles = repmat({'-'},14,1);
    line_width  = repmat(3.5,14,1);
    
    figure;
    set(gca,'FontSize',14)
    set(gca,'TickLabelInterpreter','latex')
    grid on
    
    hold on;
    for i=1:size(results,2)
        plot(horzs, results(:,i)', 'Color', line_colors(i,:), 'LineStyle', line_styles{i}, 'LineWidth', line_width(i));
        xlim([min(horzs) max(horzs)])
    end
    hold off;
    
    if ~isempty(add_line)
        hold on;
        the_xlim = xlim;
        line(the_xlim, add_line*ones(1,2), 'Color', 'k', 'LineStyle', ':');
        hold off;
        xlim(the_xlim);
    end
    title(plot_name,'interpreter','latex','FontSize',18);
    xlabel('Horizon','interpreter','latex','FontSize',16);
    legend(plot_legend, 'Location', 'southoutside', 'NumColumns', 3, 'interpreter', 'latex','FontSize',16);

end