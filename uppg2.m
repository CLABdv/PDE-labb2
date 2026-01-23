eps = 0.01;
n=10;
rs = linspace(eps,1,n);
vals = FEM(n, eps);
figure(1)
plot(rs,vals);
hold on;
%plot(rs, 1/(2*log(eps))*log(rs.^2));
hold off;

function vals = FEM(n, eps) 
    h = (1-eps)/n;
    ATilde = sparse(n,n);
    node = @(i) eps + i*h;
    xWrapper = @(x) x*(x>=eps)*(x<=1);
    for i=1:n
        for j=1:n
            if abs(i-j) == 1
                tmp = @(x) -1/h^2*xWrapper(x);
                point = (node(i)+node(j))/2;
                ATilde(i,j)=h*tmp(point);
            end
            if i == j
                tmp = @(x) 1/h^2*xWrapper(x);
                point1 = (node(i)+node(i-1))/2;
                point2 = (node(i)+node(i+1))/2;     
                ATilde(i,j) = h*(tmp(point1)+tmp(point2));
            end
        end
    end
    
    LTilde = zeros(n,1);
    tmp = @(x) -1/h^2*xWrapper(x);
    LTilde(1) = -h*tmp((node(1)+node(0))/2);
    vals = ATilde \ LTilde;
end

