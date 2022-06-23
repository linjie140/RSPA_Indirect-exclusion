clear all
N=5;F=3;c=1;
T=2;
%p=0.6;
%%%%%%%%%%%%%%%%%%%%%%%%
x=0:0.01:1;p_num=[];syms p;I=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(x)
  Pc=0;
  Pd=0;  
    for Nc=0:N-1
        for Noc=0:Nc
        Nd=N-Nc-1;
        temp=factorial(N-1)/(factorial(Noc)*factorial(Nd)*factorial(Nc-Noc))*p^Noc*(1-p)^(Nc-Noc)*(x(i))^Nc*(1-x(i))^Nd;

        pi_c=F*c*(Nc+1)/N-c+(F*c*(Noc+1)/(Nc+1)-c)*(T-1)*(1+p*(N-1))/N;
        pi_d=F*c*Nc/N+(T-1)*F*c*Noc*(1+(1-p)*(N-1))/(N*(Nc+1));

        Pc=Pc+temp*pi_c;
        Pd=Pd+temp*pi_d;
        end
    end   
dx1=Pc-Pd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   stupid = solve(dx1,p);
     if isempty(stupid) || i==1|| i==101
        continue;
    else
        I(i)=x(i);
         p_num(i,:)=stupid;
    end
end

% I=I/1;
plot(p_num(2:99),I(2:99),'*')
hold on
