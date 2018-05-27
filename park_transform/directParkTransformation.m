%% Realiza transformação de Park partindo do referencial 'abc' para o 'odq'

function u_odq = directParkTransformation(u_abc, theta_r)

% Matriz de transformação de Park
P = sqrt(2/3)*[1/sqrt(2) 1/sqrt(2) 1/sqrt(2); ...
    cos(theta_r) cos(theta_r - 2*pi/3) cos(theta_r - 4*pi/3); ...
    -sin(theta_r) -sin(theta_r - 2*pi/3) -sin(theta_r - 4*pi/3)];

% Realiza transformação
u_odq = P*u_abc;

end

