#/bin/sh
 
pyenv local 3.10.6


rm -rf *.mod


# Esse primeiro comando cria o arquivo interf.pyf, abram para olhar e tentem entender o que ele contem
python -m numpy.f2py -h interf.pyf constants.f90 allocation_reform.f90 -m module_alloc --overwrite-signature

# Este comando usa o interf.pyf e o codigo para gerar o módulo
python -m numpy.f2py interf.pyf -c constants.f90 allocation_reform.f90 

# esta ferramenta é do modulo numpy do python... se chama F2py. Leiam a documentacao!

# escrevi um pequeno script de python mostrando como usar o modulo criado
# test_module_006.py
