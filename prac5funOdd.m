function D = prac5funOdd(Ra,K0)

K = K0;
lambda = (Ra/K^4)^(1/3);

q_posi=[K*sqrt(1-lambda) K*sqrt(lambda*(1+sqrt(3)*1i)/2+1) K*sqrt(lambda*(1-sqrt(3)*1i)/2+1)];

M=zeros(3);

M(1,1)=sinh(q_posi(1)/2);
M(2,1)=q_posi(1)*cosh(q_posi(1)/2);
M(3,1)=(q_posi(1)^2-K^2)^2*(sinh(q_posi(1)/2));
        
M(1,2)=sinh(q_posi(2)/2);
M(2,2)=q_posi(2)*cosh(q_posi(2)/2);
M(3,2)=(q_posi(2)^2-K^2)^2*(sinh(q_posi(2)/2));
        
M(1,3)=sinh(q_posi(3)/2);
M(2,3)=q_posi(3)*cosh(q_posi(3)/2);
M(3,3)=(q_posi(3)^2-K^2)^2*(sinh(q_posi(3)/2));
        
D = real(det(M)); % Impose real part equal to zero
end