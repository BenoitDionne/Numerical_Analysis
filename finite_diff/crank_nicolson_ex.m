% Example:
f = @(x,y) x.*y;
gb = @(x) x.*(2-x);
gl = @(y) y.*(4-y);
gr = @(y) y.*(4-y).^2;
c = 0.5;
a = 0;
b = 2;
t1 = 0;
t2 = 4;
N = 20;
M = 40;

[X,T,U] = crank_nicolson(f,gb,gl,gr,c,N,M,a,b,t1,t2);
[XX,TT] = meshgrid(X,T);

% We need to transpose the matrix U because meshgrid()
% transposes the coordinates.
surf(XX,TT,U');
xlabel('x')
ylabel('t')
zlabel('u')

% Umax = max(max(U));
% Umin = min(min(U));
% 
% figure;
% grid on
% axis([a,b,Umin,Umax])
% for k = 1:1:M+1
%     plot(X,U(:,k),'b');
%     title(['frame ',num2str(k),' of ',num2str(M+1)]);
%     pause(2);
% end
