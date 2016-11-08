%%% -------------------------------------------------- %%%
%%% sine-Gordon equation solver on graphs              %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

close
clear all
format longE

addpath 'tools/'

% declaration of global variables:
global l M Nv V X

%%% Numerical parameters:
N  = 500;			% number of nodes in each lattice
l  = 10.0;			% physical length of the lattice
x  = linspace(0,l,N);% 1D spatial coordinate
dx = x(2) - x(1);	% discretization step in space
dt = 0.5*l/N;		% discrete time step
% Nt = 4000;			% number of time steps we will perform
Nt = 3300;
t  = 0.0;			% initial moment of time
m  = 5;				% we plot every m time steps

dt2   = dt*dt;		% dt^2
dtdx2 = (dt/dx)^2;	% coefficient in the leap-frog scheme

% Discrete Laplacian:
e = ones(N,1);		% vector of ones
L = dtdx2*spdiags([e -2*e e], -1:1, N, N);

%%% Graph definition and declaration of variables on it:
% (we consider the Mercedes-Benz graph from the draft)
M = 6;				% number of edges in the graph
Nv = 4;				% number of vertices

% we are going to do some leapfrogging, so three levels are needed to store:
u0 = zeros(N,M);	% M lattices consisting of N points
u1 = zeros(N,M);	% M lattices consisting of N points
u2 = zeros(N,M);	% M lattices consisting of N points

E = zeros(M, 2);	% condensed incidence matrix

% Edge 1:
E(1,:) = [1; 2]; % v1 -> v2

% Edge 2:
E(2,:) = [2; 3]; % v2 -> v3

% Edge 3:
E(3,:) = [2; 4]; % v2 -> v4

% Edge 4:
E(4,:) = [3; 4]; % v3 -> v4

% Edge 5:
E(5,:) = [4; 1]; % v4 -> v1

% Edge 6:
E(6,:) = [1; 3]; % v1 -> v3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initial condition specification:
c0 = 0.95;			% kink's speed

% Edge 1: kink propagating forward:
u0(:,1) = 4.0*atan(exp((x - x(N/2))/sqrt(1 - c0*c0)));
u1(:,1) = 4.0*atan(exp((x - x(N/2) - c0*dt)/sqrt(1 - c0*c0)));

% Edge 2: constant state (by continuity):
u0(:,2) = 2*pi;
u1(:,2) = 2*pi;

% Edge 3: constant state (by continuity):
u0(:,3) = 2*pi;
u1(:,3) = 2*pi;

% Edge 4: constant state (by continuity):
u0(:,4) = 2*pi;
u1(:,4) = 2*pi;

% Edge 5: forward-propagating kink:
% (kink's front is initially placed at the edge center)
u0(end:-1:1,5) = 4.0*atan(exp((x - x(N/2))/sqrt(1 - c0*c0)));
u1(end:-1:1,5) = 4.0*atan(exp((x - x(N/2) - c0*dt)/sqrt(1 - c0*c0)));

% Edge 6: backward-propagating kink:
% (kink's front is initially placed at the edge center)
u0(:,6) = 4.0*atan(exp((x - x(N/2))/sqrt(1 - c0*c0)));
u1(:,6) = 4.0*atan(exp((x - x(N/2) - c0*dt)/sqrt(1 - c0*c0)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Graph definition:
X = zeros(N,2,M);	% vector of (x,y) coordinates
V = zeros(2,Nv);	% vertices in the graph
th = linspace(0.0, 1.0, N);

% Vertice 1:
V(1,1) = 0.0;
V(2,1) = -l;

% Vertice 2:
V(1,2) = 0.0;
V(2,2) = 0.0;

% Vertice 3:
V(1,3) = l*cos(pi/6);
V(2,3) = l*sin(pi/6);

% Vertice 4:
V(1,4) = l*cos(5*pi/6);
V(2,4) = l*sin(5*pi/6);

% Vertice connectivity to lattice points:
Vn = zeros(Nv,3,2); % Nv vertices, 3 neighbors, each neighbor is coded with 2 numbers

%%% Vertice 1 and its neighbours:
Vn(1,1,:) = [2; 1]; % a1
Vn(1,2,:) = [2; 6]; % a6
Vn(1,3,:) = [N-1; 5]; % b5

%%% Vertice 2 and its neighbours:
Vn(2,1,:) = [2; 2]; % a2
Vn(2,2,:) = [2; 3]; % a3
Vn(2,3,:) = [N-1; 1]; % b1

%%% Vertice 3 and its neighbours:
Vn(3,1,:) = [2; 4]; % a4
Vn(3,2,:) = [N-1; 2]; % b2
Vn(3,3,:) = [N-1; 6]; % b6

%%% Vertice 4 and its neighbours:
Vn(4,1,:) = [2; 5]; % a5
Vn(4,2,:) = [N-1; 3]; % b3
Vn(4,3,:) = [N-1; 4]; % b4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Embedding of the graph in R^2:

% Edge 1:
X(:,1,1) = 0.0;
X(:,2,1) = -l*(1 - th);

% Edge 2:
X(:,1,2) = l*cos(pi/6)*th;
X(:,2,2) = l*sin(pi/6)*th;

% Edge 3:
X(:,1,3) = l*cos(5*pi/6)*th;
X(:,2,3) = l*sin(5*pi/6)*th;

% Edge 4:
X(:,1,4) = l*cos(pi/6*(1-th) + 5*pi/6*th);
X(:,2,4) = l*sin(pi/6*(1-th) + 5*pi/6*th);

% Edge 5:
X(:,1,5) = l*cos(5*pi/6*(1-th) + 3*pi/2*th);
X(:,2,5) = l*sin(5*pi/6*(1-th) + 3*pi/2*th);

% Edge 6:
X(:,1,6) = l*cos(3*pi/2*(1-th) + 13*pi/6*th);
X(:,2,6) = l*sin(3*pi/2*(1-th) + 13*pi/6*th);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Let's plot the initial condition:
FigHandle = figure(1);
set(FigHandle, 'Renderer', 'zbuffer');
set(FigHandle, 'Position', [100, 100, 600, 500]);

Plot(u0, t);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Nt % Main loop in time:
	t = t + dt;	% update the time variable
	for j=1:M % Loop over the graph edges:
		% reconstruct the value at the initial point:
		List = [u0(Vn(E(j,1),1,1),Vn(E(j,1),1,2));...
		 u0(Vn(E(j,1),2,1),Vn(E(j,1),2,2));...
		 u0(Vn(E(j,1),3,1),Vn(E(j,1),3,2))];
		u0(1,j) = sum(List)/length(List);

		List = [u1(Vn(E(j,1),1,1),Vn(E(j,1),1,2));...
		 u1(Vn(E(j,1),2,1),Vn(E(j,1),2,2));...
		 u1(Vn(E(j,1),3,1),Vn(E(j,1),3,2))];
		u1(1,j) = sum(List)/length(List);

		% reconstruct the value at the terminal point:
		List = [u0(Vn(E(j,2),1,1),Vn(E(j,2),1,2));...
		 u0(Vn(E(j,2),2,1),Vn(E(j,2),2,2));...
		 u0(Vn(E(j,2),3,1),Vn(E(j,2),3,2))];
		u0(N,j) = sum(List)/length(List);

		List = [u1(Vn(E(j,2),1,1),Vn(E(j,2),1,2));...
		 u1(Vn(E(j,2),2,1),Vn(E(j,2),2,2));...
		 u1(Vn(E(j,2),3,1),Vn(E(j,2),3,2))];
		u1(N,j) = sum(List)/length(List);

		u2(:,j) = 2*u1(:,j) - u0(:,j) + L*u1(:,j) - dt2*sin(u1(:,j));
	end % for j

	u0 = u1; u1 = u2;

	if (mod(i,m) == 0)
		Plot(u0, t);
	end % if ()
end % for i