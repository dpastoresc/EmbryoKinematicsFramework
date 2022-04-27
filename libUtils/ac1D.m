function output=ac1D(input)

output=zeros(1,length(input))
for i=2:length(input)1
    output(i)=input(i-1)+input(i);
end