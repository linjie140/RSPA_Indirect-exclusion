clear all;
Z=50;
M=zeros(Z+1,Z+1);
mu=0.01;s=2;
c=1;N=5;F=3;T=5;p=0.5;%%stable
%  N=5;F=3;c=1;T=2;p=0.01;%D
% N=5;F=3;c=1;T=6;p=0.8;%C stable

M(Z+1,Z+1)=1-mu;M(Z+1,Z)=mu;M(1,2)=mu;M(1,1)=1-mu;C=[];

for i=2:Z
    k=i-1;%%%population c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Pd=0;  
  Pc=0;
    for j_c=0:N-1
        for joc=0:i
       if j_c>k-1
            E1=0;
       else
            E1=nchoosek(k-1, j_c);
       end
       if (N-j_c-1)>(Z-k)
            E2=0;
       else
            E2=nchoosek(Z-k, N-j_c-1);
       end
       if j_c>k
            E3=0;
       else
            E3=nchoosek(k, j_c);
       end
       if (N-j_c-1)>(Z-k-1)
            E4=0;
       else
            E4=nchoosek(Z-k-1, N-j_c-1);
       end
       if joc>j_c
           E5=0;
       else
           E5=nchoosek(j_c, joc);
       end
       temp1=E1*E2*E5*p^joc*(1-p)^(j_c-joc)/nchoosek(Z-1, N-1);
       temp2=E3*E4*E5*p^joc*(1-p)^(j_c-joc)/nchoosek(Z-1, N-1);
        pi_c=F*c*(j_c+1)/N-c+(F*c*(joc+1)/(j_c+1)-c)*(T-1)*(1+p*(N-1))/N;
        pi_d=F*c*j_c/N+(T-1)*F*c*joc*(1+(1-p)*(N-1))/(N*(j_c+1));
        Pc=Pc+temp1*pi_c;
        Pd=Pd+temp2*pi_d;
        end
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T_jia=k*(Z-k)/((1+exp(s*(Pd-Pc)))*Z^2);
    T_jian=k*(Z-k)/((1+exp(s*(Pc-Pd)))*Z^2);
    M(i,i-1)=(1-mu)*T_jian+mu*k/Z;
    M(i,i+1)=(1-mu)*T_jia+mu*(Z-k)/Z;
    M(i,i)=1-M(i,i-1)-M(i,i+1);
end

  [V,D]=eig(M');
  [x1,x2]=find(D>0.999999);
  x=V(:,x2);
  x=x/sum(x);

plot(x)
hold on