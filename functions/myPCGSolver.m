function [z] = myPCGSolver(A,b,z,mu,v1,v2,v3,g1,g2,g3,x1,x2,x3,g11,g22,g33)

    rhs   = 2*A'*b+(g1+g2+g3)+2*mu*(v1+v2+v3)+(g11+2*mu*x1)+(g22+2*mu*x2)+(g33+2*mu*x3);                                                 
    [z,~] = pcg(@(x) Fun(x),rhs,1e-6,400,[],[],z);
  
    function y = Fun(x)
        y = 2*A'*(A*x)+6*mu*x+2*(mu+mu+mu)*x;
    end

 end