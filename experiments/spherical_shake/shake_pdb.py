# Este script crea las carpetas dependiendo de las pruebas que se quieran realizar 
# agitara los archivos pdb del archivo mainfile

import os 
import numpy as np
import random


#Funciones para crear archivos y hacer shake
def new_name(main,i,n):
    """ Crea un nombre en formato string para un archivo nuevo correspondiente 
    a las N pruebas del análisis de Montecarlo"""
    i=str(i)
    num=i
    if len(i)<len(str(n)):
        dif=len(str(n))-len(i)
        num='0'*dif+i
    return main.split('.')[0]+num+'.'+main.split('.')[1]

def randomX(x,R=1): # R indica el radio atómico en Angstroms 
    """Teniendo un vector x, se agita una esfera que lo contiene entregando una 
    nueva posición"""
    r=random.uniform(0,R)
    theta=random.uniform(0,2*np.pi)
    phi=random.uniform(0,np.pi)
    new_x=[0,0,0]
    new_x[0]=x[0]+r*np.cos(theta)*np.sin(phi)
    new_x[1]=x[1]+r*np.sin(theta)*np.sin(phi)
    new_x[2]=x[2]+r*np.cos(phi)
    return new_x

def dec3(x):
    """Redondeo de un número float a 3 decimales para correcta escritura en .pqr"""
    valor= str(round(x,3))
    if '.' in valor:
        decimal=valor.split('.')[-1]
        if len(decimal)<3:
            decimal=decimal+'000'
            decimal=decimal[:3]
            return valor.split('.')[0]+'.'+decimal
        else:
            return valor
    else:
        return valor+'.000'

if __name__=='__main__':

    #Opciones generales
    mainfile='1pgb.pdb'
    radios_prueba=[0.1,0.5,1,2,5]
    n_tests=100

    #Directorios de prueba
    for i in range(len(radios_prueba)):
        os.mkdir('tests-{}'.format(i))

    # Se guarda info archivo original

    # Lo que no es info relevante se guarda en remark, lo que va despues de atoms como ter, hetatm,master igual se queda
    remarks=""""""
    atoms=""""""
    end=""""""
    arch=open(mainfile)
    post_atoms=False
    for line in arch:
        data=line.split()
        if data[0]!='ATOM' and not post_atoms:
            remarks+=line
        elif data[0]=='ATOM':
            post_atoms=True
            atoms+=line
        else:
            end+=line
    arch.close()

    def new_file(main,i,n,path,R):
        file=open(path+new_name(main,i,n),'w')
        file.write(remarks)
        list_atoms=atoms.split('\n')
        for i in range(len(list_atoms)):
            linea=list_atoms[i].split()
            if len(linea)==0: continue
            _x=list(map(float,linea[6:9])) #5:8 pqr, 6:9 pdb
            ### ---------------------
            ### Aqui se cambia el radio de resolución en Angstroms
            ### como segundo parametro de la funcion con valores por omisión
            _x=randomX(_x,R)
            ### --------------------------------------------------
            x,y,z=list(map(dec3,_x))
            lineareal=list(list_atoms[i])
            c=0
            for j in z[::-1]+' ':
                lineareal[53-c]=j
                c+=1
            c=0
            for j in y[::-1]+' ':
                lineareal[45-c]=j
                c+=1
            c=0
            for j in x[::-1]+' ':
                lineareal[37-c]=j
                c+=1
            file.write(''.join(lineareal)+'\n')
        file.write(end)
        file.close()

    # Ahora se crean los directorios con los archivos de prueba

    for i_file in range(len(radios_prueba)):
        path='tests-'+str(i_file)+'/'
        for i in range(n_tests):
            new_file(mainfile,i,n_tests,path,radios_prueba[i_file])


