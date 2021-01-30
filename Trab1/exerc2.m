#!/home/kuroneko/snap/octave -qf
1;
function [x,er,k] = jacobi(a, b, tol, kmax)
  n = rows(a);
  for i = 1:n
    x(i) = b(i)/a(i, i)
  endfor
  k = 0;
  while (k < kmax)
    for i = 1:n
      soma = 0;
      for j = 1:n
        if (i != j)
          soma = soma + a(i,j) * x(j);
        endif
      endfor
      aux(i) = (b(i) - soma)/a(i,i);
      er(i) =  (norm(x(i),inf) - norm(aux(i),inf))/norm(x(i), inf);
      if(abs(er(i)) < tol)
        return;
      endif
      x(i) = aux(i);
    endfor
    k = k + 1;
  endwhile
endfunction

function dom = diagonal_dominante(a)
  dom = true(1);
  n = rows(a);
  for i = 1:n
    somaLinha = 0;
    for j = 1:n
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
  
  #d
  [V lambda] = eig(a);
  raioEspec = max(abs(diag(lambda)));
  
  #e
  tol = input('Insira a tolerância: ');
  kmax = input('INsira o n máximo de iterações: ');
  [xJacobi,erJacobi,kJacobi] = jacobi(a, b, tol, kmax);
  save metodoJacobi.text xJacobi erJacobi kJacobi tol kmax b a;
  #solSeidel = 
  #solSOR = 

 endfunction

addpath(pwd);
args = argv();
printf("matriz = %s\n", args{1});
analise(args{1});
printf("Fim programa");

