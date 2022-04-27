function cmap=cmap2image(descriptorName)

    name=['colormaps/' descriptorName '_cmap']
    x=load(name);
    x=struct2cell(x);
    x=x{1};
    scale=[0:1:255];
    cmap=x/255;
    A=ind2rgb(scale,cmap);
    imwrite(A,['colormaps/' descriptorName '.png']);
