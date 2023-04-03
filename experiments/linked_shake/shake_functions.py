import numpy as np
import os
import random

def importar_atomos(pqr_file):
    """Se toma la información de un archivo pqr para traspasarla a una lista pqr_list que 
    contiene la información textual en formato lista de strings de cada linea del pqr"""
    pqr_list=[]
    file=open(pqr_file)
    for linea in file: 
        if 'ATOM' not in linea: continue
        pqr_list.append(linea.strip().split())
    file.close()
    return pqr_list

def xyz_atomo(atomo,pqr_file):
    """Retorna una lista de flotantes correspondientes a las coordenadas cartesianas de cierto átomo"""
    atomos=importar_atomos(pqr_file)
    for _,_ID,a,_,_,x,y,z,_,_ in atomos:
        if atomo==a:
            return list(map(float,(x,y,z)))
    return None


def dec3(x):
    """Redondeo de un número float a 3 decimales para correcta escritura en .pqr"""
    return str(round(x,3))

def agitar_m2(pqr_file,mol2_file,dic_IDatoms,alpha=np.pi/4,p_r=0.2):
    """Retorna una lista de elementos fácilmente traspasable a un archivo pqr, 
    cada elemento es una lista correspondiente a una línea de un archivo pqr.
    
    Usa el metodo de agitacion para un porcentaje de distancia entre atomos y 
    un angulo de probabilidad."""
    pqr_list=importar_atomos(pqr_file)
    bonds_list=encontrar_conectividad(mol2_file)
    for _,origin,target,_ in bonds_list:
        atomo=dic_IDatoms[origin]
        x1,y1,z1=xyz_atomo(atomo,pqr_file)
        atomo2=dic_IDatoms[target]
        x2,y2,z2=xyz_atomo(atomo2,pqr_file)
        X,Y,Z=x2-x1,y2-y1,z2-z1
        # Se definen coordenadas cilindricas
        
        R=(X**2+Y**2+Z**2)**0.5
        
        if X>0:
            theta=np.arctan(Y/X)
        elif X==0 and Y!=0:
            theta=np.pi/2*Y/abs(Y)
        elif X<0:
            theta=np.pi + np.arctan(Y/X)
        elif X==0 and Y==0:
            theta=np.uniform(0,2*np.pi)

        phi=np.arccos(Z/R)

        n_R=random.uniform(R*(1-p_r),R*(1+p_r))
        n_theta=random.uniform(theta-alpha,theta+alpha)
        n_phi=random.uniform(phi-alpha,phi+alpha)

        #Se vuelve coordenadas cartesianas
        nX= X + n_R*np.cos(n_theta)*np.sin(n_phi)
        nY= Y + n_R*np.sin(n_theta)*np.sin(n_phi)
        nZ= Z + n_R*np.cos(n_phi)

        for i in range(len(pqr_list)):
            if pqr_list[i][2]==atomo2:
                pqr_list[i][5]=dec3(nX+x1)
                pqr_list[i][6]=dec3(nY+y1)
                pqr_list[i][7]=dec3(nZ+z1)
    return pqr_list

def new_name(main,i,n):
    """ Crea un nombre en formato string para un archivo nuevo correspondiente 
    a las N pruebas del análisis de Montecarlo"""
    i=str(i)
    num=i
    if len(i)<len(str(n)):
        dif=len(str(n))-len(i)
        num='0'*dif+i
    return main.split('.')[0]+num+'.'+main.split('.')[1]

def lista2pqr(pqr_list,nombre):
    """Teniendo una lista de datos se escribe un archivo pqr (explicitar la extensión)"""
    nuevo=open(nombre,'w')
    for linea in pqr_list:
        nuevo.write('\t'.join(linea)+'\n')
    nuevo.close()
    return None

#Funciones para conectividad de átomos

def encontrar_conectividad(mol2_file):
    """A partir de un archivo mol2 se retorna una lista con la información de la conectividad de los átomos, la estructura de cada sublista es:
        [ [bond_id origin_atom_id target_atom_id bond_type [status_bits] ], ... ]
        """    
    archivo=open(mol2_file)
    is_bond_info=False
    bonds=[]
    for linea in archivo:
        if '@<TRIPOS>BOND' not in linea and not is_bond_info: continue
        elif '@<TRIPOS>BOND' in linea:
            is_bond_info=True
        elif is_bond_info and '@<TRIPOS>' not in linea and len(linea.split())>=4:
            data=linea.strip().split()
            bonds.append(data)
        elif is_bond_info and '@<TRIPOS>' in linea:
            is_bond_info=False
            archivo.close()
            return bonds
    #En el caso de que la información relacionada a conectividad sea la última del archivo
    archivo.close()
    return bonds


def diccionario_ID_nombre(mol2_file):
    """ Retorna un diccionario cuya llave es el ID de un atomo y su valor el nombre de este (por ejemplo {1:C1,6:H3})"""
    # Se lee informacion de TRIPOS con respecto a los atomos
    archivo=open(mol2_file)
    is_atom_info=False
    dic={}
    for linea in archivo:
        if '@<TRIPOS>ATOM' not in linea and not is_atom_info: continue
        elif '@<TRIPOS>ATOM' in linea:
            is_atom_info=True
        elif is_atom_info and '@<TRIPOS>' not in linea:
            data=linea.strip().split()
            dic[data[0]]=data[1]
        elif is_atom_info and '@<TRIPOS>' in linea:
            is_atom_info=False
            archivo.close()
            return dic
    archivo.close()
    return dic