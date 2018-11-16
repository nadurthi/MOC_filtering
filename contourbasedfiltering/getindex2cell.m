function g=getindex2cell(cellstrct,val)
    for i=1:length(cellstrct)
       if all(cellstrct{i}==val) || strcmp(cellstrct{i},val)
          g=i;
           return 
       end
    end
    
    g=NaN;
end
