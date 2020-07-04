dbs = 10;
dse = 80;
dew = 70;
dwh = 40;

a = [0 0 0 0 0 0 0];
alpha = [-pi/2 pi/2 -pi/2 pi/2 -pi/2 pi/2 0];
d = [dbs 0 dse 0 dew 0 dwh];

for i = 1:7
	L(i) = Link([0 d(i) a(i) alpha(i)], 'standard');
end

robot = SerialLink(L, 'name', '7DOF-Arm');

n_joints = size(robot.qlim, 1);

q_start = 2*[-10; -10; 10; 20; 10; 10; 10];

T_start = robot.fkine(q_start);
T_end = transl([120 -90 20]) * troty(-pi/8)*trotz(pi/8)*trotx(pi/12);

n_steps = 50;
time = 5;
dt = time / n_steps;
Ttraj = ctraj(T_start, T_end, n_steps);

q_traj(1,1:7) = q_start;
q_traj_noopt(1, 1:7) = q_start;

q_current = q_start;


for i = 1:(n_steps-1)
	
	Tprev = Ttraj(:,:,i);
	Tcurrent = Ttraj(:,:,i+1);

	prev_pos(1:3) = transl(Tprev);
	prev_pos(4:6) = tr2rpy(Tprev);
	cur_pos(1:3) = transl(Tcurrent);
	cur_pos(4:6) = tr2rpy(Tcurrent);

	dVe = (cur_pos - prev_pos) / dt;

	I = eye(size(robot.qlim, 1));
	J = robot.jacob0(q_current);
	J_pinv = J' * ((J*J')^-1);
	P = I - J_pinv*J;

	dq0_dot_current = 0;
	 
	dq_dot_curr = J_pinv * dVe';
	q_current = q_current + (dq_dot_curr * dt);
	q_traj_noopt(i+1, 1:7) = q_current;
	
	%A = robot.fkine(q_current);

	if false
		lb = robot.qlim(:,1) - q_current;
		ub = robot.qlim(:,2) - q_current;
		
		options = optimoptions('fmincon', 'Display', 'off');
		obj_fun = @(x) -manipulability(x, q_current, robot, dt, P);
		dq0_dot = fmincon(obj_fun, zeros(n_joints, 1), [], [], [], [], lb, ub, [], options);
		q_current = q_current + (P * dq0_dot) * dt;
	elseif true
		lb = robot.qlim(:,1) - q_current;
		ub = robot.qlim(:,2) - q_current;
		
		options = optimoptions('fmincon', 'Display', 'off');
		obj_fun = @(x) -dist_from_limit(x, q_current, robot, dt, P);
		dq0_dot = fmincon(obj_fun, zeros(n_joints, 1), [], [], [], [], lb, ub, [], options);
		q_current = q_current + (P * dq0_dot) * dt;
	end
	
	q_traj(i+1, 1:7) = q_current;
	
	%robot.fkine(q_current) - A
end

%T_traj = robot.fkine(q_traj);
%T_traj_noopt = robot.fkine(q_traj_noopt);


robot.plot(q_traj);


%% OBJECTIVE FUNCTIONS


function m = manipulability(dq0_dot, q, robot, dt, P)

	q_updated = q + (P*dq0_dot) * dt;

	J = robot.jacob0(q_updated);
	m = sqrt(det(J*J'));
end


function lim_dist = dist_from_limit(dq0_dot, q, robot, dt, P)

	q_updated = q + (P*dq0_dot * dt);

	avg_sum = 0;
	for i = 1 : size(q_updated, 1)
		lb = robot.qlim(i,1);
		ub = robot.qlim(i,2);
		middle_val = (ub - lb) / 2;
		
		avg_sum = avg_sum + ( (q_updated(i) - middle_val) / (ub - lb) )^2; 
	end
	
	lim_dist = avg_sum / (2 * size(q_updated, 1));
end

