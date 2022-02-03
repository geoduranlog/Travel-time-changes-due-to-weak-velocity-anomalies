function [a,b,siga,yy,error,x]=regreslin(x,y)
%--test----
 %x=dist;
 %y=me'; 

 
 
 m=1;
 for i=1:length(x)
     if x(i)>1000
         
         largeoffset(m)=x(i);
    Y(m)=y(i); %me(i);
      
    m=m+1;
     end
  
 end
 
 
 %clear  largeoffset  Y
 
 x=largeoffset;
 y=Y;



n=length(x);
[p]=polyfit(x,y,1);
a=p(1); b=p(2);
yy=a*x+b; error=y-yy;
varx=std(x)*std(x);
sige=1/(n-2)*sum(error.*error);
siga=sqrt(sige/(n*varx));
%  
%  figure
%  plot(x,y,'bo',x,yy,'r-'); 
%  yyp=(a+siga)*x+b;
%  yym=(a-siga)*x+b;
%  hold on;
%  plot(x,yyp,'g-',x,yym,'k-');
 

end