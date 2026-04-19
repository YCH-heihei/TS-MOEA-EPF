function Output = Z6(PopDec,K)
%    Output=-0.7*Z4(PopDec,K)+10*(1/size(PopDec,2)*sum(abs(PopDec),2).^0.002);
   Output=0.7*Z4(PopDec,K)+10*(1/size(PopDec,2)*sum(abs(PopDec),2).^0.002);
    assert(all(all(Output<10)&all(Output>0)),'函数写错了');
end

