function[sorted] =  ksort(fdlist)
klength = length(fdlist);
for index = 1:2:klength
   sorted(ceil(index/2)) = fdlist(index);
end
return