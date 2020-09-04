function [vsol,ksol,resd]=prac5newton(F,Ra0,K0,tol,kmax)
%All vector inputs should be vertical vectors
k=0; tolk=1; X=Ra0;
while tolk>tol && k<kmax
    m=length(Ra0);
    I=eye(m);
    h=sqrt(1);
    for j=1:m
        f1=feval(F,Ra0-I(:,j)*h,K0);
        f2=feval(F,Ra0+I(:,j)*h,K0);
        DF(:,j)=(f2-f1)/(2*h);
    end
    Fk=feval(F,X(:,k+1),K0);
    deltaxk=DF\(-Fk);
    xk=X(:,k+1)+deltaxk;
    X=[X xk];
    tolk=max(abs(X(:,k+1)-X(:,k+2)));
    k=k+1;
end
vsol=X(:,end);
ksol=k;
resd=max(abs(X(:,end)-X(:,end-1)));
end