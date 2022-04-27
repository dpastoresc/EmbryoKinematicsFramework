function v2 = vector2longvector(v, factor)

      nn=sqrt(v(1)^2+v(2)^2+v(3)^2);
      v=v./nn; 
      v=int16(v.*factor);
      v2=double(v);