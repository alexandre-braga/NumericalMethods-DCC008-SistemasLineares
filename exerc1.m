#!/home/kuroneko/snap/octave -qf
1;
function analise(matriz)
  
  #a
  load("-mat",matriz);
  a = Problem.A;
  n = rows(a);
  
  #b
  [L,U,P] = lu(a);
  
  #c
  #clf;
  #spy(a);
  #spy(L);
  #spy(U);
  
  #d
  preench = 100 - nnz(a)/(nnz(L)+nnz(U)) * 100;
  printf("preench = %f\n", preench);
  
  #e
  b = a * ones(n,1); 
  x = a\b; #Método de Gauss;
  
  #f
  difX = ones(n,1) - x;
  distX = norm(difX, inf) / norm(ones(n,1), inf);
  printf("distX = %f\n", distX);
  
  #g
  difA = a - (inv(P) * L * U);
  distA = norm(difA, inf) / norm(a, inf);
  printf("distA = %f\n", distA);
  
  #h
  difB = b - ((inv(P) * L * U) * x);
  distB = norm(difB, inf) / norm(b, inf);
  printf("distB = %f\n", distB);
  
  #i
  r = norm(b - a * x, inf);
  printf("r = %f\n", r);
  
  #j
  if(n < 10000)
    K = cond(a);
    save ultimaMatriz.text x distX distA distB r preench K b a;
  else 
    save ultimaMatriz.text x distX distA distB r preench b a;
  endif
  
 endfunction

addpath(pwd);
args = argv();
printf("matriz = %s\n", args{1});
analise(args{1});
printf("Fim programa");

