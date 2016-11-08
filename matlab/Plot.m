%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

function Plot (u, t)
	global l M Nv V X

	for i=1:M
		surf([X(:,1,i), X(:,1,i)], [X(:,2,i), X(:,2,i)], [u(:,i), u(:,i)],...
			[u(:,i), u(:,i)], 'EdgeColor', 'flat',... 
			'FaceColor', 'none', 'LineWidth', 3.0), hold on
	end % for i
	for i=1:Nv
		plot3(V(1,i), V(2,i), 10.0, 'bo', 'LineWidth', 5.0), hold on
	end % for i
	axis square, grid off, hold off
	colormap(cbrewer('div', 'RdBu', 11));
	colorbar;
	xlabel(colorbar, 'u(x,t)')
	xlabel('$x$', 'interpreter', 'latex', 'fontsize', 14);
	ylabel('$y$', 'interpreter', 'latex', 'fontsize', 14, 'Rotation', 1);
	zlabel('$u(x,y,t)$', 'interpreter', 'latex', 'fontsize', 14);
	title(['sine-Gordon solution on a graph at t = ', num2str(t,'%4.2f')],...
		'interpreter', 'latex', 'fontsize', 12);
	axis([-1.1*l 1.1*l -1.1*l 1.1*l 0 2.0*pi]);
	view([0 90]);
	set(gcf, 'Color', 'w');
	drawnow
	
end % Plot ()