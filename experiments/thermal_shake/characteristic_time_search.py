# Characteristic time search

# Programa para encontrar el tiempo característico de una molécula pequeña en base a método de Montecarlo

# Librerias a utilizar
import os
import numpy as np
from time import time

# Comienza programa para ejecutar N simulaciones de Montecarlo

START_TIME=time()

# MSMS es un archivo ejecutable, en caso de estar trabajando en linux
# os.system('chmod +x bem_electrostatics/ExternalSoftware/MSMS/msms')

# Opciones de Montecarlo y directorios
##### Modificable ############
n_tests=100
path='tests'
pqr_path='mobley_test_pqr'
#############################

if path not in os.listdir(): os.mkdir(path)

# Diccionario de parámetros individuales para cada test y función para completar archivo de parámetros

def completar_archivo_entrada(dic_parametros):
    """En base a los parametros de la funcion se completa el archivo de parametros de entrada para el solver externo"""
    arch=open('input_config.prm') # Este archivo siempre se llamara asi
    temp=open('temp.txt','w')
    for linea in arch:
        if '=' in linea:
            par,_=linea.strip().split('=')
            temp.write('{}={}\n'.format(par,dic_parametros[par]))
        else:
            temp.write(linea)
    arch.close()
    temp.close()

    # Se elimina archivo original y renombra archivo temporal
    os.remove('input_config.prm')
    os.rename('temp.txt','input_config.prm')
    return None

def calcular_std_desv(archivo):
    """En base a un archivo csv con los datos de una prueba a tc especifico, retorna un flotante asociado a la desviacion estandar de las pruebas"""
    arch=open(archivo)
    arch.readline() # Skip header
    data=[]
    for linea in arch:
        G=float(linea.strip().split(',')[1])
        if not np.isnan(G): 
            data.append(G)
    data=np.array(data)
    # Aqui seria interesante agregar alguna forma de corroborar que existan una cantidad de datos minima para calcular 
    # una desviacion estandar, en caso de que hayan errores o numeros nan
    return np.std(data)

def agregar_calculo(arch_historial,linea):
    """Se agrega la linea al archivo historial de pruebas"""
    arch = open(arch_historial,'a')
    arch.write(linea)
    arch.close()
    return None

parametros = {
    'NAME':'',
    'TESTPATH':'',
    'PQRFILE':'',
    'N_TESTS':n_tests,
    'TIEMPO':0,
    'SIGMA_MOBLEY':0,
    'EP_IN':1,
    'KAPPA':0,
    'DENSIDAD':3,
    'OUTPUT_FILE':''
}

# Masas de átomos y sigmas electrostáticos
from utilities import masas, info_mobley, cantidad_atomos, moleculas #,epsilon_in # Pendiente de desarrollar

#En caso de pausar ejecución, se guarda un archivo .txt con el historial de montecarlos realizados para continuar desde el ultimo punto
hist_file=open('historial.txt','a') # si no esta creado, se crea
# Cada linea tendra la forma molecula-{N_tc} (N_tc puede ser 1 como limite inferior, 2 limite superior, 3 tc candidato)
hist_file.close()
tiempos_prueba          = [1e-9,1e-8]
nombres_tiempos_prueba  = ['T1','T2','T3']
formato_linea           = '{}-{}\n'
formato_dir             = '{}-{}'
formato_csv             = '{}-{}.csv'
path_results            = 'CSV_data'
arch_resumen            = 'resumen_tiempos.csv'

#Ejecucion programa principal

# Directorio de resultados Montecarlo
if path_results not in os.listdir(): os.mkdir(path_results)

# CSV con info de tiempos caracteristicos y parametros ecuacion
if arch_resumen not in os.listdir(): 
    resumen = open('resumen_tiempos.csv','w')
    resumen.write('Molecula,alpha,beta,N atomos,tiempo caracteristico\n')
    resumen.close()

for molecula in moleculas:
    # Actualizan parametros en diccionario
    parametros['NAME']          = molecula
    parametros['PQRFILE']       = os.path.join(pqr_path, molecula) + '.pqr' # Directorio pqr original
    parametros['SIGMA_MOBLEY']  = info_mobley[molecula][1]       # Valor de sigma objetivo
    #parametros['EP_IN'] = epsilon_in[molecula]                 # Pendiente de desarrollar, por el momento es ep=1

    # Se crean directorios para las pruebas en caso de no existir
    for nombre_t in nombres_tiempos_prueba:
        if formato_dir.format(molecula,nombre_t) not in os.listdir(path): 
            os.mkdir(os.path.join(path,formato_dir.format(molecula,nombre_t)))

    # Comprobacion en caso de que ya se hayan ejecutado pruebas de Montecarlo
    skipT1 = False
    skipT2 = False
    skipT3 = False
    hist_file=open('historial.txt')

    for linea in hist_file:
        if formato_linea.format(molecula,'T1')==linea: skipT1=True
        if formato_linea.format(molecula,'T2')==linea: skipT2=True
        if formato_linea.format(molecula,'T3')==linea: skipT3=True
    hist_file.close()
    
    i=0
    for tc in tiempos_prueba:
        # Se saltan pruebas ya realizadas
        if skipT1 and skipT2: 
            break
        elif skipT1 and not skipT2 and i==0: 
            i+=1
            continue
        # Se asume que no puede estar calculado T2 sin T1
        
        # Se ajustan parametros en diccionario
        parametros['TIEMPO']      = tc
        parametros['TESTPATH']   = os.path.join(path,formato_dir.format(molecula,nombres_tiempos_prueba[i]))
        parametros['OUTPUT_FILE'] = os.path.join(path_results,formato_csv.format(molecula,nombres_tiempos_prueba[i]))
        completar_archivo_entrada(parametros)
        # Se calculan propiedades de interés usando solver externo
        os.system('python CTS_Solver.py')
        if 'error_file.txt' in os.listdir(): break # En caso de parar la ejecucion
        agregar_calculo('historial.txt',formato_linea.format(molecula,nombres_tiempos_prueba[i]))
        i+=1
    if skipT3: continue
    
    # Calculo de desviaciones estandar para T1 y T2
    std_desv=[0,0]
    for i in range(2):
        archivo     = os.path.join(path_results, formato_csv.format(molecula,nombres_tiempos_prueba[i]))
        std_desv[i] = calcular_std_desv(archivo)

    # Calcular tc*
    alpha,beta  = np.polyfit(np.log(tiempos_prueba),np.log(std_desv),1)
    sigma_e     = info_mobley[molecula][1]
    tc          = np.exp( (np.log(sigma_e) -beta) / alpha )
    
    # Escribir datos en resumen para mantener un historial
    resumen = open(arch_resumen,'a')
    resumen.write('{},{},{},{},{}\n'.format(molecula,alpha,beta,cantidad_atomos[molecula],tc))
    resumen.close()

    # Calculos a tiempo caracteristico objetivo para lograr sigma de Mobley
    i=2
    parametros['TIEMPO']      = tc
    parametros['TESTPATH']    = os.path.join(path,formato_dir.format(molecula,nombres_tiempos_prueba[i]))
    parametros['OUTPUT_FILE'] = os.path.join(path_results,formato_csv.format(molecula,nombres_tiempos_prueba[i]))
    completar_archivo_entrada(parametros)
    os.system('python CTS_Solver.py')
    if 'error_file.txt' in os.listdir(): break # En caso de parar la ejecucion
    agregar_calculo('historial.txt',formato_linea.format(molecula,nombres_tiempos_prueba[i]))


if 'error_file.txt' in os.listdir():
        os.remove('error_file.txt')
