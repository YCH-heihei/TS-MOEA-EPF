function Output = Z4(PopDec,K)
   Output=10/(2*exp(1)-2).*(exp(max(abs(PopDec).^0.5,[],2))...
   -exp(1/size(PopDec,2).*sum(1/2*(cos((2*K-1)*pi*PopDec)+1),2))...
       -1+exp(1));
    assert(all(all(Output<10)&all(Output>0)),'函数写错了');
end

