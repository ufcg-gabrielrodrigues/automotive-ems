%% Modelo matemático em referencial DQ para alternador Lundell

function [i_d, i_q] = alternatorDQModel(e_d, e_q, v_d, v_q, r_s, l_d, l_q, ...
    omega_r, dt)

%
persistent past_i_d past_i_q past_l_d past_l_q;

%
if (isempty(past_i_d))
    past_i_d = 0.0;
    past_i_q = 0.0;
    past_l_d = l_d;
    past_l_q = l_q;
end

%
Z = [(r_s + (l_d - past_l_d)/dt + l_d/dt) (-omega_r*l_q); ...
    (omega_r*l_d) (r_s + (l_q - past_l_q)/dt + l_q/dt)];

%
V = [(e_d - v_d + (l_d*past_i_d)/dt); (e_q - v_q + (l_q*past_i_q)/dt)];

%
I = Z\V;

%
i_d = I(1);
i_q = I(2);

end

