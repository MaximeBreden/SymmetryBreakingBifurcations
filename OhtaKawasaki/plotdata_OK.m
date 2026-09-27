function plotdata_OK(u, phi, para, n)

% This produces plots of the steady state u, and of the eigenfunction phi.

if nargin < 4
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

ueval = eval_cos(u, X, Y);
phieval = eval_cos(phi, X, Y);

x = xresc*para.Lx/pi;

if is1D
    figure
    plot(x, ueval, 'b', 'Linewidth', 2)
    hold on
    plot(x, phieval, '--c', 'Linewidth', 2)
    xlabel('$x$','Interpreter','latex')
    set(gca,'FontSize',15) 
    axis tight
    legend('$u$','$\varphi$','Interpreter','latex','Location','best')
else
    y = yresc*para.Ly/pi;

    figure;
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
    
    % Tile 1
    ax1 = nexttile;
    imagesc(x, y, ueval);
    axis(ax1, 'xy', 'equal', 'tight');
    colormap(ax1, parula);
    colorbar(ax1);
    xlabel(ax1, '$x$', 'Interpreter', 'latex'); 
    ylabel(ax1, '$y$', 'Interpreter', 'latex');
    title(ax1, '$u$', 'Interpreter', 'latex');
    set(gca,'FontSize', 20) 
    axis tight
    
    % Tile 2
    ax2 = nexttile;
    imagesc(x, y, phieval);
    axis(ax2, 'xy', 'equal', 'tight');
    colormap(ax2, parula);
    colorbar(ax2);
    xlabel(ax2, '$x$', 'Interpreter', 'latex'); 
    ylabel(ax2, '$y$', 'Interpreter', 'latex');
    title(ax2, '$\varphi$', 'Interpreter', 'latex');
    set(gca,'FontSize', 20) 
    axis tight
end
