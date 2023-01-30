import module_alloc as alloc_module

woodens = alloc_module.traits.dwood
spec_la = alloc_module.traits.sla

print(woodens)

p1 = 2.4
p2 = 1.8
p3 = 19
p4 = 1
p5 = 3
p6 = 15
p7 = 3.4
# p8 = 2.4
# p9 = 2.4

def outputs(out):
    lst = ["lm2", "rm2", "cw2"]

    return dict(zip(lst, out))

for x in range(5):
    # chama a funcao e guarda os outputs na variavel res
    res = alloc_module.allocation_reform.alloc(p1,p2,p3,p4,p5,p6,p7,woodens,spec_la)
    
    # res tem os 11 elementos que saem da subrotina de fortran
    # nessa linha eu atualizo as variáveis que em teoria teriam que ser atualizadas.
    
    print(p1)

    
    resul = outputs(res)
    leaf_mass = resul['lm2'] #KgC/m2
    root_mass = resul['rm2'] #KgC/m2
    wood_mass = resul['cw2'] #KgC/m2

    print(leaf_mass, x)
    print(root_mass, x)
    print(wood_mass, x)




    

#     #ATUALIZANDO PARA PROXIMO STEP
#     # leaf_mass,root_mass,wood_mass,p4,p5,p6,p7,woodens,spec_la,a,b = res
