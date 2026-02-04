
eps = 0.1;
n = 10;
delta = 0.0044;

h = (1-eps)/n;
ATilde = sparse(n+1,n+1);
node = @(i) eps + i*h;
xWrapper = @(x) x*(x>=eps)*(x<=1);
for i=0:n
        for j=0:n
            if abs(i-j) == 1
                tmp = @(x) -1/h^2*xWrapper(x);
                point = (node(i)+node(j))/2;
                ATilde(i+1,j+1)=h*tmp(point);
            end
            if i == j
                tmp = @(x) 1/h^2*xWrapper(x);
                point1 = (node(i)+node(i-1))/2;
                point2 = (node(i)+node(i+1))/2;     
                ATilde(i+1,j+1) = h*(tmp(point1)+tmp(point2));
            end
        end
end
full(ATilde)
%%
rs = linspace(eps, 1, n+1 );
truesol = 1/(2*log(eps))*log(rs.^2);

m = 8;
ds = zeros(1,m);
steps = zeros(1,m);
errors = zeros(1,m);

for i = 1:m
d = delta*2^(-i);
[u, step, error] = GradientDescent(n, zeros(n-1, 1), d, ATilde, 1e-4, 2500000, truesol');
u = [1; u ; 0];
ds(i) = d;
steps(i) = step;
errors(i) = error;
end
A = [ds; steps; errors]
array2table(A)

plot(rs, u, "DisplayName", "gradient descent");
hold on;
plot(rs, truesol, "DisplayName", "true sol");
legend();
hold off;
steps
error
any(isnan(u(:)))



function [u, steps, error] = GradientDescent(n, u0, delta, A, tolerance, maxsteps, truesol)
    
    u = u0;
   

    %f = @(v)  [1 v' 0]*A*[1 ; v ; 0];
    steps = 0;
    while steps < maxsteps


        steps = steps+1;
        %prevu = u;
        grad = 2*A*[1; u; 0];
        grad = grad(2 : n);
        u = u- delta*grad;
        
        error = norm([1; u ; 0]-truesol, inf);
        if error < tolerance
            break
        end



    end

end