%% Realiza transformação de Park partindo do referencial 'odq' para o 'abc'

function u_abc = inverseParkTransformation(u_odq, theta_r)

% Matriz de transformação inversa de Park
Q = sqrt(2/3)*[1/sqrt(2) cos(theta_r) -sin(theta_r); ...
    1/sqrt(2) cos(theta_r - 2*pi/3) -sin(theta_r - 2*pi/3); ...
    1/sqrt(2) cos(theta_r - 4*pi/3) -sin(theta_r - 4*pi/3)];

% Realiza transformação
u_abc = Q*u_odq;

end

