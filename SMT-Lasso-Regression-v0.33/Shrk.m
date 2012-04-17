function [bshrk]=Shrk(b,val,option)

if strcmp(option,'soft')
  bshrk=sign(b).*max(abs(b)-val,0);
elseif strcmp(option,'hard')
  bshrk=b;
  index=abs(b)<val;
  bshrk(index)=0;
else
  fprintf(1,'Wrong option');
  bshrk=0;
end
