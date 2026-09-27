function plotdata_SKT(u1, u2, phi1, phi2, para, n)

% This produces plots of the steady state (u1,u2), and of the eigenfunction
% (phi1,phi2).

if nargin < 6
    n = [200, 200]; % spatial discretization for the plots
end

is1D = not(isfield(para,'Ly')); % check whether the solution is 1D or 2D
if is1D
    if length(n) == 2
        n(2) = 1;
    else
        n = [n, 1];
    end
end

xresc = linspace(0,pi,n(1));
yresc = linspace(0,pi,n(2));
[X, Y] = meshgrid(xresc, yresc);

u1eval = eval_cos(u1, X, Y);
u2eval = eval_cos(u2, X, Y);
phi1eval = eval_cos(phi1, X, Y);
phi2eval = eval_cos(phi2, X, Y);

x = xresc*para.Lx/pi;

if is1D
    figure
    plot(x, u1eval, 'b', 'Linewidth', 2)
    hold on
    plot(x, u2eval, 'r', 'Linewidth', 2)
    plot(x, phi1eval, '--c', 'Linewidth', 2)
    plot(x, phi2eval, '--m', 'Linewidth', 2)
    xlabel('$x$','Interpreter','latex')
    set(gca,'FontSize',15) 
    axis tight
    legend('$u_1$','$u_2$','$\varphi_1$','$\varphi_2$','Interpreter','latex','Location','best')

else
    y = yresc*para.Ly/pi;

    figure
    tiledlayout(1,4,'Padding','compact','TileSpacing','compact');
    
    % Tile 1
    ax1 = nexttile;
    imagesc(x, y, u1eval);
    axis(ax1, 'xy', 'equal', 'tight');
    colormap(ax1, parula);
    colorbar(ax1);
    xlabel(ax1, '$x$', 'Interpreter', 'latex'); 
    ylabel(ax1, '$y$', 'Interpreter', 'latex');
    title(ax1, '$u_1$', 'Interpreter', 'latex');
    set(gca,'FontSize', 20) 
    axis tight
    
    % Tile 2
    ax2 = nexttile;
    imagesc(x, y, u2eval);
    axis(ax2, 'xy', 'equal', 'tight');
    colormap(ax2, parula);
    colorbar(ax2);
    xlabel(ax2, '$x$', 'Interpreter', 'latex'); 
    ylabel(ax2, '$y$', 'Interpreter', 'latex');
    title(ax2, '$u_2$', 'Interpreter', 'latex');
    set(gca,'FontSize', 20) 
    axis tight
    
    % Tile 3
    ax3 = nexttile;
    imagesc(x, y, phi1eval);
    axis(ax3, 'xy', 'equal', 'tight');
    colormap(ax3, parula);
    colorbar(ax3);
    xlabel(ax3, '$x$', 'Interpreter', 'latex'); 
    ylabel(ax3, '$y$', 'Interpreter', 'latex');
    title(ax3, '$\varphi_1$', 'Interpreter', 'latex');
    set(gca,'FontSize', 20) 
    axis tight
    
    % Tile 4
    ax4 = nexttile;
    imagesc(x, y, phi2eval);
    axis(ax4, 'xy', 'equal', 'tight');
    colormap(ax4, parula);
    colorbar(ax4);
    xlabel(ax4, '$x$', 'Interpreter', 'latex'); 
    ylabel(ax4, '$y$', 'Interpreter', 'latex');
    title(ax4, '$\varphi_2$', 'Interpreter', 'latex');
    set(gca,'FontSize', 20) 
    axis tight
end
