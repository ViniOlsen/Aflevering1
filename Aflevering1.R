#Steepest Decent, 
#Steepest Decent is defind by going in the negative gradient direction, 
#where a positive skalar alpha is called the step lenght
#the first I have done i defind the function, the algorithm later will use
#the gradient, nomal R can solve it with the function grand in the numDeriv-package
#unfortunate it do not work for me, so I had to solve the gradient myself which is 
#defind as g1 (x[1]) and g2 (x[2])
# I have choosen my epsilon (tolerance) to be 10e-10.


f <- function(x){
  100*(x[2]-x[1]^2)^2+(1-x[1])^2 }
g1 <- function(x){ (-400)*(x[2]-x[1]^2)*x[1]-2*(1-x[1])}
g2 <- function(x){ 200*(x[2]-x[1]^2)}
x=c(-1.2, 1)
E=10e-10

#I have called my algorithm for steep.
#k defind how many step my algorithm use to find a solution.
#The are to loop in the algorithm, the frist i to find a new point to start with. 
#the other loop is to find the best solution for alpha.

steep <- function(f, x, E) { 
  k=0
  a=1
  d=c(E+1,E+1)
  while(sqrt((a*d[1])^2+(a*d[2])^2) > E) {
    d= -c(g1(x),g2(x));
    r=0.5;
    c=10e-4;
    while(f(x+a*d)>f(x)+c*a* sum(t(d)*d)) {a=r*a};
    cat("k = ", k, "; f =", f(x), "@ x =", x, "\n")
    x=x+a*d;
    k=k+1
  }
}
steep(f, x, E)

#Newton metode is defind by going in the negative Hessian direction, 
#where a positive skalar alpha is called the step lenght
#the first I have done i defind the function, the algorithm later will use
#the gradient and Hessian, nomal R can solve the gradient with the function grand in the numDeriv-package
#unfortunate it do not work for me, so I had to solve the gradient and hessian myself which is 
#defind as g1 (x[1]) and g2 (x[2]) for the gradient and, h1 (X[1]x[1]), h2 (X[1]x[2])
# h3 (X[2]x[1]) and (X[2]x[2]) and put those in a matrix H.
# I have choosen my epsilon (tolerance) to be 10e-10.
h1 <- function(x){1/(400*x[1]^2-400*x[2]-2)}
h2 <- function(x){x[1]/(200*x[1]^2-200*x[2]-1)}
h3 <- function(x){x[1]/(200*x[1]^2-200*x[2]-1)}
h4 <- function(x){(-600*x[1]^2+200*x[2]+1)/(200*(-200*x[1]^2+200*x[2]+1))}
H <- matrix(c(h1(x), h2(x), h3(x), h4(x)), nrow = 2, ncol = 2, byrow=TRUE)

#I have called my algorithm for newton.
#k defind how many step my algorithm use to find a solution.
#The are to loop in the algorithm, the frist i to find a new point to start with. 
#the other loop is to find the best solution for alpha.

newton <- function(f, x, E) {
  k=0
  a=1
  d=c(E+1,E+1)
  while(sqrt((a*d[1])^2+(a*d[2])^2) > E) {
    d= -H*c(g1(x),g2(x));
    r=0.5;
    c=10e-4;
    while(f(x+a*d)>f(x)+c*a* sum(t(d)*d)) {a=r*a};
    cat("k = ", k, "; f =", f(x), "@ x =", x, "\n")
    x=x+a*d;
    k=k+1  
  }
}
newton(f, x, E)
