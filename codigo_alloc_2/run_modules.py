import matplotlib.pyplot as plt
import module_006 as my_module


leaf_in = 0.0001
root_in = 0.0001
sap_in = 0.001
heart_in = 0.5
dens_in = 10
bminc_in = 2.5

res = my_module.allocation.alloc(leaf_in, root_in, sap_in, heart_in, bminc_in,dens_in)


#assessing a variable:
#t = my_module.program2

#print(t)

# #

# p1 = 100000 # kg/m2
# p2 = 300000


# sample = []
# x = 0
# tmp = 0


# while True:
    
#     if x > 3: break
    
#     # chama a funcao e guarda os outputs na variavel res
#     res = my_module.establishment.establish(p1,p2) #colocar só os inputs

#     r1 = res
#     sample.append(r1)
#     x += 1
#     p1 += 250
        
# plt.plot(sample)
# plt.savefig("test_module006.png")
# plt.clf() # limpa o grfico
 
