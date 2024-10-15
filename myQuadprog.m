function xstar = myQuadprog(H,f,A,b)
% my quadprog solver that uses manual matrix solving methods
Q = [H A.'; A zeros(length(A(:,1)))];

if ~iscolumn(f)
    f = f(:);
end
if ~iscolumn(b)
    b = b(:);
end
z = [f;b];

% solve Qx=z

% Linear solve matlab
% xstar = Q\z;
% xstar = xstar(1:length(H(:,1)));

% LU factorization
% TODO write my own LU
[L,U,P] = lu(Q);
y = forwardSubstitution(L,P*z);
xstar = backwardSubstitution(U,y);
xstar = xstar(1:length(H(:,1)));

% QR factorization
% TODO write my own QR
% [Q,R] = qr(Q);
% y = Q.'*z;
% xstar = backwardSubstitution(R,y);
% xstar = xstar(1:length(H(:,1)));
end

function x = forwardSubstitution(L, b) % GENERATED WITH CHATGPT
% forward substitution function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: This function was generated using chatgpt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Forward substitution for lower triangular system L*x = b
    
    % Check if L is a square matrix
    [m, n] = size(L);
    if m ~= n
        error('Matrix L must be square.');
    end
    
    % Check if L is lower triangular
    if ~istril(L)
        keyboard
        error('Matrix L must be lower triangular.');
    end
    
    % Check if the size of L and b are compatible
    if size(b, 1) ~= n
        error('Incompatible dimensions between L and b.');
    end
    
    % Initialize solution vector x
    x = zeros(n, 1);
    
    % Perform forward substitution
    for i = 1:n
        x(i) = (b(i) - L(i, 1:i-1)*x(1:i-1)) / L(i, i);
    end
end

function x = backwardSubstitution(U, b) % GENERATED WITH CHATGPT
    % Backward substitution for upper triangular system U*x = b

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE: This function was generated using chatgpt
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Check if U is a square matrix
    [m, n] = size(U);
    if m ~= n
        error('Matrix U must be square.');
    end
    
    % Check if U is upper triangular
    if ~istriu(U)
        error('Matrix U must be upper triangular.');
    end
    
    % Check if the size of U and b are compatible
    if size(b, 1) ~= n
        error('Incompatible dimensions between U and b.');
    end
    
    % Initialize solution vector x
    x = zeros(n, 1);
    
    % Perform backward substitution
    for i = n:-1:1
        x(i) = (b(i) - U(i, i+1:end)*x(i+1:end)) / U(i, i);
    end
end

function [Q,R] = myHouseholderQR(A)

N = length(A(:,1));
Q = eye(N);

for j = 1:N
    x = A(j:end,j);
    alpha = norm(x);
    sgn = sign(x(1));
    alpha = -sgn*alpha;
    e = zeros(length(x),1);
    e(1) = 1;
    u = x-alpha*e;
    v = u/norm(u);
    Qkp = eye(length(x))-2*(v*v');

    if length(Qkp(:,1))<N
        Qk = [eye(N-length(x)) zeros(N-length(x),length(x));
            zeros(length(x),N-length(x)) Qkp];
    else
        Qk = Qkp;
    end
    Q = Qk*Q;
    A = Qk*A;
end
R = A;
Q = Q';
end

