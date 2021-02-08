#!/home/kuroneko/snap/octave -qf
1;
function [x,er,k] = jacobi(a, b, tol, kmax)
  n = rows(a);
  x = zeros(n,1);
  aux = zeros(n,1);
  for i = 1:n,
    x(i) = b(i)/a(i, i);
  endfor
  for k = 1:kmax,
    for i = 1:n,
      soma = 0;
      for j = 1:n,
        if (i != j)
          soma = soma + a(i,j) * x(j);
        endif
      endfor
      aux(i) = x(i);
      x(i) = (b(i) - soma)/a(i,i);
    endfor
    er(k)  = (norm(x, inf) - norm(aux, inf))/norm(x, inf);
    if(abs(er(k)) < tol)
      return;
    endif
  endfor
  return;
endfunction

function [x,er,k] = sor(a, b, tol, kmax, w)
  n = rows(a);
  x = zeros(n,1);
  xAnt = zeros(n,1);
  for i = 1:n,
    x(i) = b(i)/a(i, i);
  endfor
  inf = tril(a, -1);
  sup = triu(a, 1);
  
  for k = 1:kmax,
    xAnt = x;
    for i = 1:n
      somainf = 0;
      somasup = 0;
      for j = 1:i,
        somainf = somainf + inf(i,j) * x(j);
      endfor
      for j = i:n,
        somasup = somasup + sup(i,j) * xAnt(j);
      endfor
      x(i) = (1-w) * xAnt(i) + w * (b(i) - somainf - somasup)/a(i,i);
    endfor
    er(k) = max(x - xAnt)/max(x);
    if(abs(er(k)) < tol)
      return;
    endif
  endfor
  return;
endfunction

function dom = diagonal_dominante(a)
  dom = true(1);
  n = rows(a);
  
  for i = 1:n
    somaLinha = 0;
    for j = 1:n
      if(i != j)
        if(a(i,j) != 0)
          somaLinha = somaLinha + a(i, j);
        endif
      endif
    endfor
    if(somaLinha >= a(i, i))
      dom = false(1);
      return;
    endif
  endfor
  return;
  
endfunction

function dom = diagonal_dominanteCholesky(a)
  try chol(a)
    dom = true(1);
  catch ME
    dom = false(1);
  end
  return
endfunction


function [BJ, BGS, BSOR] = fatora(a, w)
  L = tril(a,-1);
  U = triu(a, 1);
  D = diag(diag(a));
  BJ  = inv(D) * (-L -U);
  BGS = inv(L + D)*(-U);
  BSOR = inv(L + (1/w)*D )*( (1/w)*D -D -U);
  return;
endfunction

function lambdaB = raioEspec(B,n)
  [V lambda] = eig(B);
  lambdaB = max(abs(diag(lambda)));
  return;
endfunction  

function analise(matriz)
  
  #a
  load("-mat", matriz);
  a = Problem.A;
  n = rows(a);
  
  #b
  b = a * ones(n,1);
  
  #c
  if(n >= 10000)
    dom = diagonal_dominanteCholesky(a);
  else
    dom = diagonal_dominante(a);
  endif
  printf("É diagonal dominante ? %d\n", dom);
  
  #d
  
  w = input('Insira o parâmetro de relaxação W: ');
  [BJ,BGS,BSOR] = fatora(a, w);
  save fatoracoesB.text BJ BGS BSOR;
  
  #e
  if(n <= 10000)
    reJacobi = raioEspec(BJ,n);
    reSeidel = raioEspec(BGS,n);
    reSOR = raioEspec(BSOR,n);
    printf("reJacobi: %d\n reSeidel: %d\n reSOR: %d\n", reJacobi, reSeidel, reSOR);
    save raiosEspecW w reJacobi reSeidel reSOR;
  else
    if(dom)
      reJacobi = 0;
      reSeidel = 0;
      reSOR = 0;
     else
      reJacobi = 2;
      reSeidel = 2;
      reSOR = 2;
      endif
  endif
  
 
  
  if(reJacobi < 1)
    printf("Método Jacobi: \n");
    tol = input('Insira a tolerância: ');
    kmax = input('Insira o n máximo de iterações: ');
    [xJacobi,erJacobi,kJacobi] = jacobi(a, b, tol, kmax);
    save metodoJacobi.text erJacobi kJacobi tol kmax reJacobi xJacobi;
  else
    XJacobi = 0;
  endif
  
  if(reSeidel < 1)
    printf("Método Seidel: \n");
    tol = input('Insira a tolerância: ');
    kmax = input('Insira o n máximo de iterações: ');
    [xSeidel,erSeidel,kSeidel] = sor(a, b, tol, kmax, 1);
    save metodoSeidel.text erSeidel kSeidel tol kmax reSeidel xSeidel;
  else
    XSeidel = 0;
  endif
  
  if(reSOR < 1)
    printf("Método SOR: \n");
    tol = input('Insira a tolerância: ');
    kmax = input('Insira o n máximo de iterações: ');
    [xSOR,erSOR,kSOR] = sor(a, b, tol, kmax, w);
    save metodoSOR.text erSOR kSOR tol kmax w reSOR xSOR;
  else
    xSOR = 0;
  endif
  
  save metodosIterativos.text xJacobi xSeidel xSOR; 
  
  #f
  kj = 1:1:kJacobi;
  plot(kj,log(erJacobi));
  hold on;
  
  kgs = 1:1:kSeidel;
  plot(kgs,log(erSeidel));
  hold on;
  
  ks = 1:1:kSOR;
  plot(ks,log(erSOR));
  
  ylabel("log(er)");
  xlabel("k");
  title("Gráfico k x log(er)");
  
  endif
endfunction

addpath(pwd);
args = argv();
printf("matriz = %s\n", args{1});
analise(args{1});
printf("Fim programa\n");

