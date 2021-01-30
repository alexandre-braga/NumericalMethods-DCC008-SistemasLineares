#!/home/kuroneko/snap/octave -qf
1;
function dom = diagonal_dominante(a)
  dom = true(1);
  n = rows(a);
  m = columns(a);
  for i = 1:n
    somaLinha = 0;
    for j = 1:m
      if(i != j)
        somaLinha = somaLinha + a(i, j);
      endif
    endfor
    if(somaLinha > a(i, i))
      dom = false(1);
      return;
    endif
  endfor
  return;
endfunction

function analise(matriz)
  
  #a
  load("-mat",matriz);
  a = Problem.A
  n = rows(a);
  
  #b
  b = a * ones(n,1);
  
  #c
  dom = diagonal_dominante(a);
  printf("É diagonal dominante ? %d\n", dom);
  
  #save MetodosIterativos.text a b dom;

 endfunction

addpath(pwd);
args = argv();
printf("matriz = %s\n", args{1});
analise(args{1});
printf("Fim programa");

