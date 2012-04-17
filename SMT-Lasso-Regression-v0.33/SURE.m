function threshold=SURE(b,n)

ab=sort(abs(b),'ascend');
p=length(b);

for i=1:p
  t = ab(i);
  Rhat(i) = p/n - (2/n)*sum(ab<=t) + sum( (min(ab,t)).^2 );
end

[value,index] = min(Rhat);

if index==p
  threshold = ab(index);
else
  threshold = (ab(index) + ab(index+1))/2;
end

%figure,plot(ab)

