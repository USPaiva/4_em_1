from sympy import * # Necessário para calcular a deriva
import numpy as np
from numpy import log10 as log
from numpy import log as ln
import pylab as pl
import math 
#############################################
def f(x):
        return (x**3-5)
def f1(x):
        return (x**3) # Aqui definir quem será f1(x)
def f2(x):
        return (5*x**0) # Aqui definir quem será f1(x)
####################################################################################
a = float(input('Entre o valor mínimo do intervalo, a = '))
b = float(input('Entre o valor máximo do intervalo, b = '))
if f(a)*f(b) > 0:
 print('Não se cumpre o Teorema de Bolzano, precisa redefinir o intervalo [a, b]')
else:
 print('Teorema de Bolzano satisfeito, o programa continuará na procura da raiz')
j= float(input("insira o metodo de resolução que deseja, sendo:\n1=método de newton.\n2=bisseção.\n3=falsa posição.\n4=secante.\n"))

##########################  Resolução Newton    ###############################################################
if j== 1:
    y = Symbol('y') # Definir variável x como símbolo para efetuar a derivada
    f3 = f(y) # Em f3 fica a função f na notação tradicional
    f_prime = f3.diff(y) # Efetua a primeira derivada da função
    f_2prime = f_prime.diff(y) # Efetua a segunda derivada da função
    f3 = lambdify(y, f3)
    f_prime = lambdify(y, f_prime)
    f_2prime = lambdify(y, f_2prime)
    a1 = a
    b1 = b
    cond1 = f3(a1)*f_2prime(a1)
    cond2 = f3(b1)*f_2prime(b1)
    if (f3(a1)*f3(b1) < 0 and f_prime(a1)*f_prime(b1) >0 and f_2prime(a1)*f_2prime(b1) > 0):
        print('Haverá convergência porque as condições de Raphson e Fourier são satisfeitas')
        if cond1 > 0:
            X0 = a1
        else:
            X0 = b1
            erro = float(input('Informe o valor do erro epsilon:'))
            i = 0
            c = b1 - a1
        h=float(input("informe se deseja um erro:\n1-menor ou igual.\n2-somente igual.\n3-somente menor que o erro:\n"))
        if h==1:
            while (abs(c) > erro or math.fabs(f3(X0)) > erro):
                X = X0 - (f3(X0)/f_prime(X0))
                c = X - X0
                X0 = X
                print(X0, c)
                i = i + 1
                if i == 100:
                    break
                else:
                    print('Não haverá convergência do método nesse intervalo ({0}, {1})'.format(a1,b1))
        elif h==2:
            while (abs(c) == erro or math.fabs(f3(X0)) == erro):
                X = X0 - (f3(X0)/f_prime(X0))
                c = X - X0
                X0 = X
                print(X0, c)
                i = i + 1
                if i == 100:
                    break
                else:
                    print('Não haverá convergência do método nesse intervalo ({0}, {1})'.format(a1,b1))
        elif h==3:
            while (abs(c) >= erro or math.fabs(f3(X0)) >= erro):
                X = X0 - (f3(X0)/f_prime(X0))
                c = X - X0
                X0 = X
                print(X0, c)
                i = i + 1
                if i == 100:
                    break
                else:
                    print('Não haverá convergência do método nesse intervalo ({0}, {1})'.format(a1,b1))
    ErroRelativo = c / X * 100 
    print('raiz Xn = %f' %X0)
    print('f(Xn) = %f' %f3(X0))
    print('iterações n = %i' %(i))
    print('|xn - a| = {0} < erro ou |f(xn)| = {1} < erro '.format(round(c,5), round(abs(f3(X0)),5)))
    print('|xn - a|/xn = {0} %'.format(round(ErroRelativo,4)))

########################   RESOLUÇÃO BISSEÇÃO   #################################################
elif j == 2:
    a1 = a
    b1 = b
    erro = float(input('Informe o valor do erro epsilon:'))
    k = (log(b1-a1)-log(erro))/log(2) #número de iterações, ver slider 21 da aula 2:
    print('Estimativa do números de iterações é k = %i' %k)
    i=0
    c = b1 - a1
    x0 = (a1 + b1)/2
    while (c > erro or math.fabs(f(x0)) > erro):
        if (f(a1)*f(x0)) < 0:
            b1 = x0
        else:
            a1 = x0
        x = (a1 + b1)/2
        c = x - x0
        x0 = x
        i = i + 1
        if i == 100:
            break
    ErroRelativo = c / x * 100
    print('raiz Xn = %f' %x)
    print('f(Xn) = %f' %f(x))
    print('iterações n = %i' %(i))
    print('|xn - a| = {0} < erro ou |f(xn)| = {1} < erro '.format(round(c,5), round(abs(f(x)),5)))
    print('|xn - a|/xn = {0} %'.format(round(ErroRelativo,4)))

###################    FALSA POSIÇÃO   ##############################################################
elif j == 3:
    a1 = a
    b1 = b
    erro = float(input('Informe o valor do erro epsilon '))
    i = 0
    c = b1 - a1
    x0 = a1 - (((b1 - a1)*f(a1))/(f(b1)-f(a1)))
    while (abs(c) > erro or math.fabs(f(x0)) > erro):
        if (f(a1)*f(x0)) < 0.0:
            b1 = x0
        else:
            a1 = x0
        x = a1 - (((b1 - a1)*f(a1))/(f(b1)-f(a1)))
        c = x - x0
        x0 = x
        print(x0, c)
        i = i + 1
        if i == 100:
            break
    ErroRelativo = c / x * 100
    print('raiz Xn = %f' %x)
    print('f(Xn) = %f' %f(x))
    print('iterações n = %i' %(i))
    print('|xn - a| = {0} < erro ou |f(xn)| = {1} < erro '.format(round(c,5), round(abs(f(x)),5)))
    print('|xn - a|/xn = {0} %'.format(round(ErroRelativo,4)))

##################    SECANTE    ###################################################################
elif j==4:
    a1 = a
    b1 = b
    erro = float(input('Informe o valor do erro epsilon:'))
    c = b1 - a1
    x2 = b1
    x1 = a1
    i = 0
    while (abs(c) > erro or math.fabs(f(x1)) > erro):
        x2 = x1 - (((x1 - x0)*f(x1))/(f(x1)-f(x0)))
        c = x2 - x1
        x0 = x1
        x1 = x2
        print(x2, c)
        i = i + 1
        if i == 100:
            break
    ErroRelativo = c / x2 * 100
    print('raiz Xn = %f' %x2)
    print('f(Xn) = %f' %f(x2))
    print('iterações n = %i' %(i))
    print('|xn - a| = {0} < erro ou |f(xn)| = {1} < erro '.format(round(c,5), round(abs(f(x2)),5)))
    print('|xn - a|/xn = {0} %'.format(round(ErroRelativo,4)))

#######################    GRAFICOS    ###########################################
xl = np.linspace(a, b, 256, endpoint=True)
Cl= f(xl)
pl.plot(xl, Cl)
pl.show()
#############################################################
xl = np.linspace(a, b, 256, endpoint=True)
Cl, S = f1(xl), f2(xl) # Aqui escrever f1(x) e f2(x)
pl.plot(xl, Cl)
pl.plot(xl, S)
pl.show()