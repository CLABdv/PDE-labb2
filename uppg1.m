a = @(x) exp(x);
f = @(x) exp(x);
g = 1;
n = 20;

X = linspace(0,1,n+1);
figure(1);
sol = FEM(a,f,g,n,true);
plot(X,sol);

minm=6;
maxm = 14;
nodes=2^minm;
diffs = zeros(maxm-minm,2^minm);
ns = zeros(maxm-minm,1);
for i = minm:maxm-1
    n = 2^(i);
    ns(i-minm+1)=n;
    tmp = FEM(a,f,g,n,false)'; %change to two
    diffs(i-minm+1, :) = tmp(n/nodes:n/nodes:n);
end
n=2^16;
tmp = FEM(a,f,g,n,false)'; %change to two
for i = maxm-minm:-1:1
   diffs(i,:) = (diffs(i,:)-tmp(n/nodes:n/nodes:n)).^2;
   %diffs(i,:) = (diffs(i,:)-diffs(i-1,:)).^2;
end
%diffs(1,:)=[];
figure(2)
loglog(ns, sum(diffs,2))

function vals = FEM(a,f,g,nNodes,twoPoints)
    h = 1/nNodes;
    ATilde = sparse(nNodes,nNodes);
    aWrapper = @(x) a(x)*(x<1);
    for i=1:nNodes
        for j=1:nNodes
            if (i-j) == -1
                tmp = @(x) -1/h^2*aWrapper(x);
                if twoPoints
                    point1 = (j+i)/2*h - h/(2*sqrt(3));
                    point2 = (j+i)/2*h + h/(2*sqrt(3));
                    ATilde(i,j) = h/2*(tmp(point1)+tmp(point2));
                else
                    point = (j+i)/2*h;
                    ATilde(i,j)=h*tmp(point);
                end
            end
            if i == j
                tmp = @(x) 1/h^2*aWrapper(x);
                if twoPoints
                    point1 = (j+i-1)/2*h - h/(2*sqrt(3));
                    point2 = (j+i-1)/2*h + h/(2*sqrt(3));
                    
                    point3 = (j+i+1)/2*h - h/(2*sqrt(3));
                    point4 = (j+i+1)/2*h + h/(2*sqrt(3));
                    ATilde(i,j) = h/2*(tmp(point1)+tmp(point2)) + h/2*(tmp(point3)+tmp(point4));
                else
                    point1 = (j+i-1)/2*h;
                    point2 = (j+i+1)/2*h;
                       
                    ATilde(i,j) = h*(tmp(point1)+tmp(point2));
                end
            end
            if (i-j) == 1
                tmp = @(x) -1/h^2*aWrapper(x);
                if twoPoints
                    point1 = (j+i)/2*h - h/(2*sqrt(3));
                    point2 = (j+i)/2*h + h/(2*sqrt(3));
                    ATilde(i,j) = h/2*(tmp(point1)+tmp(point2));
                else
                    point = (j+i)/2*h;
                    ATilde(i,j) = h*tmp(point);
                end
            end
        end
    end
    
    LTilde = zeros(nNodes,1);
    for i = 1:nNodes
                tmp = @(x) f(x)*(1-abs(x-i*h)/h);
        if twoPoints
                point1 = (i+i-1)/2*h - h/(2*sqrt(3));
                point2 = (i+i-1)/2*h + h/(2*sqrt(3));
                
                point3 = (i+i+1)/2*h - h/(2*sqrt(3));
                point4 = (i+i+1)/2*h + h/(2*sqrt(3));
                LTilde(i) = h/2*(tmp(point1)+tmp(point2)) + h/2*(tmp(point3)+tmp(point4));
        else
            point1 = (i+1)*h/2;
            point2 = (i-1)*h/2;

            LTilde(i) = h*(tmp(point1)+tmp(point2));
        end
    end
    LTilde(nNodes)=LTilde(nNodes)+g;
    vals = [0;ATilde \ LTilde];
end
