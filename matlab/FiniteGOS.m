clear all
c=1;F=3;
Z=50;
% p=0.01;T=2;
p=0.6;T=5;
% p=0.8;T=6;
N=5;%M:total_population_size:N:participations_in_the_games;
mu=0.01;%mu:exporation_rate;s:strong imitation;sigmun:nonparticipators;
beta=2;
i_c=0:1:Z;

deltaic=[];
for i=1:size(i_c,2)
  Pd=0;  
  Pc=0;
    for j_c=0:N-1
        for joc=0:i
       if j_c>i_c(i)-1
            E1=0;
       else
            E1=nchoosek(i_c(i)-1, j_c);
       end
       if (N-j_c-1)>(Z-i_c(i))
            E2=0;
       else
            E2=nchoosek(Z-i_c(i), N-j_c-1);
       end
       if j_c>i_c(i)
            E3=0;
       else
            E3=nchoosek(i_c(i), j_c);
       end
       if (N-j_c-1)>(Z-i_c(i)-1)
            E4=0;
       else
            E4=nchoosek(Z-i_c(i)-1, N-j_c-1);
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

        deltaic(i)=i_c(i)/Z*((Z-i_c(i))/(Z-1))*(1/(1+exp(beta*(Pd-Pc)))-1/(1+exp(beta*(Pc-Pd))));
end

    plot(i_c/Z,deltaic,'*')
    hold on