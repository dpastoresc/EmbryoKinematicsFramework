function output=der1D(input)

if size(input,1)==1
    input=input';
end
output=zeros(size(input));
for i=1:size(input,1)-1
    output(i,:)=input(i+1,:)-input(i,:);
end
if size(input,2)==1
    output=output';
end