%
% Complement to get geometric jacobian
% Inputs: dq_kinematics and joint configuration
% Output: 6xn geometric jacobian 
%
function [J] = human_geom(kine, dq_jaco, q, dof)
J = zeros(6,dof);
v = [-1, 1, 1, 1, -1, 1, 1, 1];
C8 = diag(v);
C4m = -C8(1:4, 1:4);
CJ4_2_J3  = [0, 1, 0, 0;0, 0, 1, 0;0, 0, 0, 1];
xm = kine.fkm(q,dof-1);
J(1:3, 1:dof) = CJ4_2_J3 * 2 * xm.P().conj().haminus4() * poseJacobian.topRows(4);
end