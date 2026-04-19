function Output = g10(PopDec,n)
%  UNTITLED 此处显示有关此函数的摘要
%  此处显示详细说明
   %% 哦哦，终于理解他说不能借助相邻个体的最优是因为这个吗(相邻方向的最优解差异巨大)
   % 那这样感觉根本就是在随机乱找到了，不同方向之间的最优解完全不具备参考性。
   m=size(PopDec,2);
   Output=zeros(size(PopDec,1),n-m);
   j=2*pi*((m+1:n)-(m+1))/n;
%    Output=repmat(1/m*sum(PopDec,2),1,n-m);
   for i=m+1:n
      Output(:,i-m)=1/2*(1/(m.^3)*sum(sin((4*PopDec-2)*pi+j(i-m)),2).^3+1);
   end
end

