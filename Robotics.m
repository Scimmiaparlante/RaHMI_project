%% ROBOT CREATION

a = [0 50 0 70 0 50 0];
alpha = [pi/2 0 pi/2 pi/2 pi/2 pi/2 0];
d = [30 0 0 0 0 0 0];
limits = [-pi/2 pi/2; -pi/2 pi/2; -pi/2 pi/2; -pi/2 pi/2; -pi/2 pi/2; -pi/2 pi/2; -pi/2 pi/2];

n_links = size(limits, 1);

for i = 1:n_links
	L(i) = Link([0 d(i) a(i) alpha(i)], 'standard');
	L(i).qlim = limits(i,:);
end

robot = SerialLink(L, 'name', 'Robot');

%robot.plot([-pi/2 pi/4 -pi/8 0 pi/8 0 0]);
%robot.qlim

%% MOTION PLANNING

%"manip" for manipulability optimization
%"limits" for workspace limits optimization
%"noopt" for no optimization
optimization = "limits";


n_joints = size(robot.qlim, 1);

q_start = [-pi/2 pi/4 -pi/8 0 pi/8 pi/4 0]';

T_start = robot.fkine(q_start);
T_end = transl([-100 -90 20]) *trotx(-pi);

n_steps = 50;
time = 5;
dt = time / n_steps;
T_traj = ctraj(T_start, T_end, n_steps);


q_traj(1,1:7) = q_start;
q_traj_noopt(1, 1:7) = q_start;

manip_traj(1) = manipulability(0, q_start, robot, 0, 0);
manip_traj_noopt(1) = manipulability(0, q_start, robot, 0, 0);

limit_traj(1) = dist_from_limit(0, q_start, robot, 0, 0);
limit_traj_noopt(1) = dist_from_limit(0, q_start, robot, 0, 0);

manip_matrix_traj(1,:,:) = manipulab_ellipsoid_matrix(robot, q_start);


q_current = q_start;


for i = 1:(n_steps-1)
	
	T_prev = T_traj(:,:,i);
	T_current = T_traj(:,:,i+1);

	prev_pos(1:3) = transl(T_prev);
	prev_pos(4:6) = tr2rpy(T_prev);
	cur_pos(1:3) = transl(T_current);
	cur_pos(4:6) = tr2rpy(T_current);

	dS = (cur_pos - prev_pos);
	
	%remove the pi multiples from the angular part of dVe
	for j = 4:6
		while abs(dS(j)) > (2*pi - 0.1)  %0.1 is a tolerance threshold
			dS(j) = sign(dS(j))*2*pi - dS(j);
		end
	end
	
	dv = dS / dt;
	
	%transform the rpy angle ratio to angular velocity, since we use the geometric jacobian 
	dv(4:6) = rpy2jac(prev_pos(4:6)) * dv(4:6)';
	
	I = eye(size(robot.qlim, 1));
	J = robot.jacob0(q_current);
	J_pinv = J' * ((J*J')^-1);
	P = I - J_pinv*J;

	dq0_dot_current = 0;
	 
	dq_dot_curr = J_pinv * dv';
	q_current = q_current + (dq_dot_curr * dt);

	if optimization ~= "noopt"
		
		lb = robot.qlim(:,1) - q_current;
		ub = robot.qlim(:,2) - q_current;
		
		options = optimoptions('fmincon', 'Display', 'off');

		if optimization == "manip"
			obj_fun = @(x) -manipulability(x, q_current, robot, dt, P);
		elseif optimization == "limits"
			obj_fun = @(x) -dist_from_limit(x, q_current, robot, dt, P);
		end
		
		dq0_dot = fmincon(obj_fun, zeros(n_joints, 1), [], [], [], [], lb, ub, [], options);
		q_current = q_current + (P * dq0_dot) * dt;
	
	else
		q_traj_noopt(i+1, 1:7) = q_current;
		manip_traj_noopt(i+1) = manipulability(0, q_current, robot, 0, 0);
		limit_traj_noopt(i+1) = dist_from_limit(0, q_current, robot, 0, 0);
	end
	
	q_traj(i+1, 1:7) = q_current;
	manip_traj(i+1) = manipulability(0, q_current, robot, 0, 0);
	limit_traj(i+1) = dist_from_limit(0, q_current, robot, 0, 0);
	manip_matrix_traj(i+1,:,:) = manipulab_ellipsoid_matrix(robot, q_current);

end

%T_traj = robot.fkine(q_traj)
%T_traj_noopt = robot.fkine(q_traj_noopt)


%robot.plot(q_traj);

%% SIMULATION WITH ELLIPSOID

for i = 1:(n_steps-1)

	robot.plot(q_traj(i,:));
	delete(h);
	hold on
	ellipsoid_matrix(:,:) = manip_matrix_traj(i,:,:);
	h = plot_ellipse(10000000*ellipsoid_matrix, transl(robot.fkine(q_traj(i,:))), 'edgecolor', 'y');
	pause(0.5);
	
end

%% Q PLOTS (OPT vs NO_OPT)

for i = 1:n_links
	figure(i)
	plot(1:50, q_traj(:,i), 'b');
	hold on
	plot(1:50, q_traj_noopt(:,i), 'r');
end

%% MANIPULABILITY PLOT

plot(1:50, manip_traj(:), 'b');
hold on
plot(1:50, manip_traj_noopt(:), 'r');

%% LIMITS DISTANCE PLOT

plot(1:50, limit_traj(:), 'b');
hold on
plot(1:50, limit_traj_noopt(:), 'r');


%% PATH PLOT

points = transl(T_traj);

robot.plot(q_start');
hold on
plot3(points(:,1), points(:,2), points(:,3), 'r');

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


%% OTHER UTILITY FUNCTIONS

function M = manipulab_ellipsoid_matrix(robot, q)

	J = robot.jacob0(q);
	
	M = inv(J(1:3,:)*J(1:3,:)');

end