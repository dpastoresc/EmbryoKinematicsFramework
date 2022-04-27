function cmap=loadCmap(name)


          x=load([ name '_cmap']) 
            x=struct2cell(x);
            x=x{1};
           % scale=[0:1:255];
            cmap=x/255;
      