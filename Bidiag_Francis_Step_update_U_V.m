function [U, Bi, V] = Bidiag_Francis_Step(Bi)
% Perform a single step of the "bidiagonal Francis Step" algorithm on the
% bidiagonal matrix Bi to introduce the bulge and chase it out.

[m, n] = size(Bi);
if m <= 2 || n <= 2
    % Nothing to do
    U = eye(m);
    V = eye(n);
    return
end

% Construct Givens rotations to zero the (m-2, m-1) element of Bi
c = Bi(m-1, m-1) / norm([Bi(m-2, m-1), Bi(m-1, m-1)]);
s = Bi(m-2, m-1) / norm([Bi(m-2, m-1), Bi(m-1, m-1)]);
G = [c, s; -s, c];

% Apply G to Bi from the left
Bi(1:m-2, m-1:m) = G * Bi(1:m-2, m-1:m);

% Construct Givens rotations to zero the (1, 2) element of the updated Bi
c = Bi(1, 1) / norm([Bi(1, 1), Bi(2, 1)]);
s = Bi(2, 1) / norm([Bi(1, 1), Bi(2, 1)]);
H = [c, s; -s, c];

% Apply H to Bi from the right
Bi(1:2, :) = Bi(1:2, :) * H';

% Construct U and V matrices to apply G and H to
U = eye(m);
U(1:m-2, 1:m-2) = eye(m-2) - (eye(m-2) - G' * G);
V = eye(n);
V(1:2, 1:2) = H;

end
