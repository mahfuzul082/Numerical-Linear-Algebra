function x = forward_backward_sub(R, b)

n = size(R, 1); % no. of rows
R_t = R'; % transpose of R

% initialize y and x vectors
y = zeros(n, 1); 
x = zeros(n, 1);

% forward substitution to get y
for i = 1:n
    y(i) = (b(i) - R_t(i, 1:i-1)*y(1:i-1))/R_t(i, i);
end

fprintf("y = \n");
disp(y);

% backward substitution to get x
for i = n:-1:1
    x(i) = (y(i) - R(i,i+1:n)*x(i+1:n)) / R(i,i);
end
end