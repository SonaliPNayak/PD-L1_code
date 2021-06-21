
% hill function
function H1 = hill(X,X0,lamb,n) 
% Hill function =  H- plus (lambda)* H+
H1= ((X0^n)/(X0^n+X^n)) + lamb*((X^n)/(X0^n+X^n));

function H2 = hilla(X,X0,lamb,n) 
% Hill function =  H- plus (lambda)* H+
H2= (1/lamb)*((X0^n)/(X0^n+X^n)) + ((X^n)/(X0^n+X^n));
