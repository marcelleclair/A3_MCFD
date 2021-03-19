function [Fx, Fy] = getForces(x,y,Ex,Ey,ds)

q_0 = 1.60217653e-19; % elementary charge

x_index = uint16(x./ds) + 1;
y_index = uint16(y./ds) + 1;
lin_index = y_index + length(y_index)*(1-x_index);
Fx = Ex(lin_index) .* q_0;
Fy = Ey(lin_index) .* q_0;
end