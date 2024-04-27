% A generic template for computing Cart-Pole Dynamics using QP Control

syms l M m g real  % Parameters/Parametric variables of the dynamic system
syms x theta dx dtheta u real % State and input variables of the system

[A,B,ddq] = computeRobotDynamicsSymbolic();

% Substituting the system parameters with physical measurement values
% Mass of cart(M) = 1.0 [kg]
% Mass of pole(m) = 0.2 [kg]
% length of pole(l) = 0.3 [m]
% acceleration due to gravity(g) = 9.81 [m/s^2]
global params

A = subs(A,[M,m,l,g],[1.0,0.2,0.3,9.81]);

B = subs(B,[M,m,l,g],[1.0,0.2,0.3,9.81]);

ddq = subs(ddq,[M,m,l,g],[1.0,0.2,0.3,9.81]);

A_fun = matlabFunction(A);

params.A_Mat = A_fun(0,0,0); % Desired States dtheta, dx, x, theta

%disp(A_afterComputation);

B_fun = matlabFunction(B);

params.B_vec = B_fun(0);

ddq_fun = matlabFunction(ddq);


% State cost : Minimizing deviation from x and x_dot is more important
% given 100 i.e., higher value compared to 0.01 for theta and theta_dot.
Qmpc = diag([100,0.01,100,0.01]); 

% Control cost : 
Rmpc = 1; 


%params.Ka = lqr(A_Mat,B_vec,Qmpc,Rmpc);

x0=[0;pi/6;0;0]; % x0 is the intial state of the system

tspan= 0:0.01:2; % simulation time

[t,x]=ode45(@(t,x) sys_dynamics(t,x,ddq_fun),tspan,x0);

% recreate control inputs
 for i=1:length(t)
     u(:,i)=controller(t(i),x(i,:)');
 end

% Plot x(t) and theta(t)
figure;
subplot(3, 1, 1);
plot(t, x(:,1));
xlabel('Time (s)');
ylabel('x (m)');
title('Cart Position over Time for case 1c with MPC controller.');
grid on;

subplot(3, 1, 2);
plot(t, x(:,2));
xlabel('Time (s)');
ylabel('\theta (rad)');
title('Pendulum Angle over Time for case 1c with MPC controller.');
grid on;

subplot(3,1,3);
plot(t,u);
title('control input');
title('Control input for case 1c with MPC controller.');
xlabel('Time (s)');
ylabel('u(N)')
grid on;

function dx=sys_dynamics(t,x,ddq_fun)

u=controller(t,x);

dq = x(3:4);

ddq = ddq_fun(x(4),x(2),u);

dx= [dq;ddq];
end


function u=controller(t,x)

global params

% Parameters for MPC
N = 20; % Number of horizons
dt = 0.05; % sampling time / time step
A_mat = params.A_Mat;
params.Qmpc = diag([100,0.01,100,0.01]);
params.Rmpc = 1;

% Discretize
Ak = A_mat * dt + eye(size(A_mat));
Bk = params.B_vec * dt ;

% Setup cost function
H = blkdiag(kron(eye(N), params.Qmpc), kron(eye(N), params.Rmpc)); % [100 * 100] 
f = zeros(1, max(size(H)));

% Equality Constraints.
Aeq2 = kron(eye(N),-Bk);
Aeq1 = kron(diag(ones(N-1,1),-1),-Ak)+eye(4*N);
Aeq = [Aeq1,Aeq2];

beq = [Ak * x;zeros(4*(N-1),1)];

% get inequality constraints
x_constr = 0.8;
th_constr = pi/4;
u_constr = 10;

Aqp1 = kron(eye(N),[diag([1 1 0 0]);diag([-1 -1 0 0])]);
Aqp2 = kron(eye(N),[1;-1]);
Aqp = blkdiag(Aqp1,Aqp2);

bqp = [repmat([x_constr;th_constr;0;0],2*N,1);u_constr*ones(2*N,1)];
options = optimoptions('quadprog','Display','off');
X_star = quadprog(H,f,Aqp,bqp,Aeq,beq,[],[],[],options);
u = X_star(4*N+1);
end



function [A,B,ddq] = computeRobotDynamicsSymbolic()

% Symbolically note down the variables

syms l M m g real % Parameters/Parametric variables of the dynamic system
syms x theta dx dtheta u real % State and input variables of the system

% State vectors; State of the dynamic system: X = [q;dq]
% =[x;theta;dx;dtheta]

q = [x; theta];
dq = [dx; dtheta];
u_vec = [u;0];

% x_cm_pole = x - l*sin(theta); y_cm_pole = l*cos(theta); 

vxpole = dx - l * cos(theta) * dtheta;

vypole = -l*sin(theta)*dtheta;

% Kinetic Energy
T = (0.5 * M * dx^2) + (0.5 * m * (vxpole^2+vypole^2)); % K.E of cart + K.E. of pole

% Potential Energy(P.E)
% P.E of cart = 0, so U = 0 + P.E of pole 
U = m*g*l*cos(theta);

% Lagrangian
L_q_dq = simplify(T - U);

f_q_dq = simplify(jacobian(L_q_dq,dq));

D_q = simplify(jacobian(f_q_dq,dq)); % Inertia Matrix % Say Internal components affecting states.

N_q_dq = simplify(jacobian(f_q_dq,q)*dq-jacobian(L_q_dq,q)');% Say External components affecting states.

ddq = simplify(D_q\(u_vec-N_q_dq));

f_x_u = [dq;ddq];

X = [q;dq];

A = simplify(jacobian(f_x_u,X));

B = simplify(jacobian(f_x_u,u));

end
