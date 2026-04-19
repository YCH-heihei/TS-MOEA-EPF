function Output = Z3(PopDec,K)
   Output=10/size(PopDec,2)*sum(1/3*(PopDec.^2-cos((2*K-1)*pi*PopDec)+1),2);
   assert(find(any(Output>10)||any(Output<10)),'函数写错了');
end

