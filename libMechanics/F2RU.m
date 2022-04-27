function [R U] =  F2RU(F)
            %calculate U => F=RU
           %calculate eigenvalues C and U
           I=[1 0 0; 0 1 0; 0 0 1];
           C=F'*F;
           [CVectors CValues]=eig(C);
           Uvalues=abs(sqrt(CValues));
           l1=Uvalues(1,1); l2 =Uvalues(2,2); l3=Uvalues(3,3);
           %I1=trace(CGR);       %I2=0.5*(I1^2-trace(CGR*CGR));       %I3=CGRdet;        
           i1=l1+l2+l3;i2=l1*l2+l1*l3+l2*l3;i3=l1*l2*l3;
           Dc=i1*i2-i3;
           %calculate U and R.
           U=(1/Dc)*(-(C*C)+(i1^2-i2)*C+i1*i3*I);
           %invU=inv(U);
           uU=unique(U);
%            uU
%            size(uU)
%            U
%            any(isnan(U(:)))
%            size(any(isnan(U(:))))
%            U=[0 0 0; 0 0 0; 0 0 0]
%                uU=unique(U);
           if size(uU,1)==1 && uU(1)==0 
               %'all zeros'
               R=I;
           elseif any(isnan(U(:)))
               R=I;
               U=I;
           elseif any(isinf(U(:)))
               R=I;
               U=I;
           else
               R=F/U;%matlab says it is faster
           end