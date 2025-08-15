import random
import numpy as np

AMINOACIDOS = {
    'A': 0.0825,  # Alanina
    'R': 0.0553,  # Arginina
    'N': 0.0406,  # Asparagina
    'D': 0.0546,  # Ácido aspártico
    'C': 0.0137,  # Cisteína
    'Q': 0.0393,  # Glutamina
    'E': 0.0675,  # Ácido glutámico
    'G': 0.0707,  # Glicina
    'H': 0.0227,  # Histidina
    'I': 0.0596,  # Isoleucina
    'L': 0.0968,  # Leucina (ajustado +0.0002)
    'K': 0.0584,  # Lisina
    'M': 0.0242,  # Metionina
    'F': 0.0386,  # Fenilalanina
    'P': 0.0470,  # Prolina
    'S': 0.0656,  # Serina
    'T': 0.0534,  # Treonina
    'W': 0.0108,  # Triptófano
    'Y': 0.0292,  # Tirosina
    'V': 0.0695   # Valina 
}

# 1A) Proteínas aleatorias
#   I) Escriba un código que genere secuencias de proteínas al azar, de un largo preestablecido por el
#    usuario. 
#   II) Modifique el mismo tal que sea posible fijar la frecuencia o probabilidad de ocurrencia de cada 
#   uno de los 20 aminoácidos.
#   Ayuda: Puede definir un diccionario con frecuencias y usar numpy.random.choice().

def random_prot():
    largo_prot = int(input('Seleccione el largo de la proteína de salida:\n'))
    prot = ''
    for i in range(largo_prot):
        amino = random.choice(list(AMINOACIDOS.keys()))
        prot += amino
    
    print('Proteína aleatoria:\n')
    print(prot)

def random_prot_frecuencia():
    largo_prot = int(input('Seleccione el largo de la proteína de salida:\n'))
    
    aminoacidos = list(AMINOACIDOS.keys())
    frecuencias = list(AMINOACIDOS.values())

    prot = ''
    for i in range(largo_prot):
        aminoacido_aleatorio = np.random.choice(aminoacidos, p=frecuencias)
        prot += aminoacido_aleatorio
    
    print('Proteína aleatoria:\n')
    print(prot)

def main():
    comando = input('Seleccione R para proteina random y F para una proteína con frecuencias tenidas en cuenta:\n')
    if comando == 'R': random_prot()
    elif comando == 'F': random_prot_frecuencia()
    else: print('Comando incorrecto.') 

main()

# 1B) Frecuencias en secuencias
#   I) Determine la frecuencia de cada uno de los 20 aminoácidos para un set de secuencias de 
#       proteínas y para una secuencia de ADN, ambas al azar, a partir de los codones correspondientes. 
#   II) Que el código sea capaz de mostrar los resultados en un gráfico de barras.
#
#   III) Ahora analice y discuta los resultados obtenidos, ¿Cómo se comparan la distribución de 
#   aminoácidos respecto las generadas para proteínas al azar?

