function detect_array=sphere_interface_scatter(center_x,center_y,center_z,r,n)
theta = linspace(0, 2*pi, n);    % 方位角（0到2π）
phi = linspace(0, pi, n);        % 极角（0到π）
[Theta, Phi] = meshgrid(theta, phi);  % 生成网格
x = r * sin(Phi) .* cos(Theta);
y = r * sin(Phi) .* sin(Theta);
z = r * cos(Phi);
detect_array=[x(:)+center_x y(:)+center_y z(:)+center_z];