 function[X]=rTV_2D(Pn,lambda)
%% This code is written by Metin Ertas, PhD on 01.07.2022. 
% Unlike the traditional TV which uses one neighboring pixel to be used in
% the gradient calculation for each partial derivative; in the proposed reinforced method rTV, another neighboring pixel adjacent to the gradient direction is considered in76
%the cost function and the gradient is still calculated on the central pixel where the derivation is taken on.

% The minimization problem is solved by using a gradient descent method.

%Inputs: Pn--- The image to be regularized.
%        lambda -- which tunes the impact of regularization.(The higher the lambda the higher the impact of denoising and vice versa)

%Output: Denoised Image.


[x,y]=size(Pn); %


 %% Map which will be used to get the real image while leaving expanded edges to "0"
map=zeros(4+x,4+y);
map(3:(x+2),3:(y+2))=1; 

%%%%%%%%%%%%%%%% Image Expansion %%%%%%%%%%%%%%%%%%%%%%%   
 X=extension(Pn,5);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Pn1=X;
error=10; % Random value to start the while loop
Grad=zeros(x+4,y+4);
ite=1;
if sum(sum(abs(X)))>1
 while error>0.0001&&ite<20

Y=X;
   Grad=2*(X-Pn1);
   for i=3:(x+2)
      for j=3:(y+2) 
          A=X(i,j)-X(i+1,j); B=X(i,j)-X(i+2,j); C=X(i,j)-X(i,j+1); D=X(i,j)-X(i,j+2);
          E=X(i-2,j)-X(i-1,j); F=X(i-2,j)-X(i,j); G=X(i-2,j)-X(i-2,j+1); H=X(i-2,j)-X(i-2,j+2);
          J=X(i-1,j)-X(i,j); K=X(i-1,j)-X(i+1,j); L=X(i-1,j)-X(i-1,j+1); M=X(i-1,j)-X(i-1,j+2);
          N=X(i,j-2)-X(i,j-1); O=X(i,j-2)-X(i,j); P=X(i,j-2)-X(i+1,j-2); R=X(i,j-2)-X(i+2,j-2);
          S=X(i,j-1)-X(i,j); T=X(i,j-1)-X(i,j+1); U=X(i,j-1)-X(i+1,j-1); Y=X(i,j-1)-X(i+2,j-1);
          
          Grad1=(2*(A+B)+2*(C+D))/ sqrt((A+B)^2+(C+D)^2+eps);
          Grad2= (E+F)/sqrt((E+F)^2+(G+H)^2+eps);
          Grad3= (J+K)/sqrt((J+K)^2+(L+M)^2+eps);
          Grad4= (N+O)/sqrt((N+O)^2+(P+R)^2+eps);
          Grad5= (S+T)/sqrt((S+T)^2+(U+Y)^2+eps);
          Grad(i,j)=Grad(i,j)+lambda*(Grad1-Grad2-Grad3-Grad4-Grad5);
                                     
      end
   end

   C_1=0+eps;
   C_2=0;
   nn=1;
   
   while C_1>C_2
       
       a=(X-Pn1).*map;
       C_1=sqrt(sum(sum((a).*(a))));
       for i=3:(x+2)
           for j=3:(y+2)
%              
                C_1=C_1+lambda*sqrt( (X(i,j)-X(i+1,j)+X(i,j)-X(i+2,j))^2+(X(i,j)-X(i,j+1)+X(i,j)-X(i,j+2))^2);
            
           end
       end
        
       X=X-10^(-6)*(2^(nn))*Grad;
       nn=nn+1;
       a=(X-Pn1).*map;
       C_2=sqrt(sum(sum((a).*(a))));
       for i=3:(x+2)
           for j=3:(y+2)
              
                C_2=C_2+lambda*sqrt( (X(i,j)-X(i+1,j)+X(i,j)-X(i+2,j))^2+(X(i,j)-X(i,j+1)+X(i,j)-X(i,j+2))^2);
             
           end
       end
   end
   
       X=X+10^(-6)*(2^(nn-1))*Grad;
       
   error=sum(abs(Y(:)-X(:)))/((x+2)*(y+2));
   ite=ite+1;
 end

end
 
 
 Pn1=X;
 X=zeros(x,y);

X(:,:)=Pn1(3:x+2,3:y+2);


