import matplotlib.pyplot as plt
import module_006 as my_module

latosa = my_module.allometric_params.latosa

print(latosa)

p1 = 1e-5 # kg/m2
p2 = 1e-5
p3 = 1e-5
p4 = 1e-5
p5 = 1e-5
p6 = 1e-5
p7 = 1e-5
p8 = 1e-5
p9 = 1e-5

sample = []
x = 0
tmp = 0


while True:
    
    if x > 3000: break
    if abs(tmp - p6) <= 1e-6:
        print(x)
        break
    
    # chama a funcao e guarda os outputs na variavel res
    tmp = p6
    res = my_module.allocation.allocate(p1,p2,p3,p4,p5,p6,p7,p8,p9)
    
    # tuple (lista imutável)
    # res tem os 11 elementos que saem da subrotina de fortran
    # nessa linha eu atualizo as variáveis que em teoria teriam que ser atualizadas.
    
    p1,p2,p3,p4,p5,p6,p7,p8,p9,a,b = res
    sample.append(p6)
    x += 1
        
plt.plot(sample)
plt.savefig("test_module006.png")
plt.clf() # limpa o grfico
 