import math
import matplotlib

matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import numpy as np


def precisao(x, indice_max, extradorso, intradorso, coefextra, coefintra):
    j = 0
    while x[indice_max] > float(extradorso[j][0]):
        j += 1
    x_preciso = np.linspace(float(extradorso[j - 1][0]), float(extradorso[j + 1][0]), num=500)
    yextra_preciso = []
    yintra_preciso = []
    for x_i in x_preciso:
        yextra_preciso.append(calculay(coefextra, j, extradorso, x_i))
        yintra_preciso.append(calculay(coefintra, j, intradorso, x_i))
    indice, novaespessura = espessura(yextra_preciso, yintra_preciso)
    return indice, novaespessura


def espessura(extradorso, intradorso):
    # Encontra a espessura máxima e a sua posição
    # Como a espessura máxima se encontra entre 10% e 20% do perfil, procuraremos somente para x entre 0.1 e 0.2
    # equivalente a procurar somente para x com indíces de 100 a 150
    max_espessura = 0
    for i in range(len(extradorso)):
        dif = extradorso[i] - intradorso[i]
        if dif > max_espessura:
            max_espessura = dif
            max_espessura_indice = i
    return max_espessura_indice, max_espessura


def coeficientes(n, coordenadas):
    # input
    a_i = []  # for i = 0, 1, ..., n
    for i in range(n + 1):
        a_i.append(float(coordenadas[i][1]))

    # step 1
    h_i = []  # for i = 0, 1, ..., n-1
    for i in range(n):
        h_i.append(float(coordenadas[i + 1][0]) - float(coordenadas[i][0]))

    # step 2
    alpha_i = ["null"] * (n)  # for i = 1, 2, ..., n-1
    for i in range(1, n):
        alpha = (3 / h_i[i]) * (a_i[i + 1] - a_i[i]) - (3 / h_i[i - 1]) * (a_i[i] - a_i[i - 1])
        alpha_i[i] = alpha

    # step 3
    l_i = [0] * (n + 1)  # for i = 0, 1, ..., n
    mi_i = [0] * n  # for i = 0, 1, ..., n-1
    z_i = [0] * (n + 1)  # for i = 0, 1, ..., n
    l_i[0] = 1
    mi_i[0] = 0
    z_i[0] = 0

    # step 4
    for i in range(1, n):
        l_i[i] = 2 * (float(coordenadas[i + 1][0]) - float(coordenadas[i - 1][0])) - h_i[i - 1] * mi_i[i - 1]
        mi_i[i] = h_i[i] / l_i[i]
        z_i[i] = (alpha_i[i] - h_i[i - 1] * z_i[i - 1]) / l_i[i]

    # step 5
    l_i[n] = 1
    z_i[n] = 0
    c_i = [0] * (n + 1)  # for i = 0, 1, ..., n
    c_i[n] = 0

    # step 6
    b_i = [0] * (n + 1)  # for i = 0, 1, ..., n
    d_i = [0] * (n + 1)  # for i = 0, 1, ..., n
    for j in range(n - 1, -1, -1):
        c_i[j] = z_i[j] - mi_i[j] * c_i[j + 1]
        b_i[j] = (a_i[j + 1] - a_i[j]) / h_i[j] - h_i[j] * (c_i[j + 1] + 2 * c_i[j]) / 3
        d_i[j] = (c_i[j + 1] - c_i[j]) / (3 * h_i[j])

    # step 7
    output = [a_i, b_i, c_i, d_i]
    return output


def calculay(coeficientes, indice, tabela, ponto):
    # indice - 1 <= ponto <= indice
    j = indice - 1
    # j <= ponto <= j + 1
    a = coeficientes[0][j]
    b = coeficientes[1][j]
    c = coeficientes[2][j]
    d = coeficientes[3][j]
    dif = ponto - float(tabela[j][0])
    y = a + b * (dif) + c * (dif) ** 2 + d * (dif) ** 3
    return y


def interpola(coeficientes, tabela, pontos):
    x_tabelado = []
    y_calculado = []
    for coordenada in tabela:
        x_tabelado.append(float(coordenada[0]))
    x_tabelado.sort()
    for ponto in pontos:
        k = 0
        while ponto > x_tabelado[k]:
            k += 1
        y_calculado.append(calculay(coeficientes, k, tabela, ponto))
    return y_calculado


def plota(x, y1, y2):
    plt.plot(x, y1)
    plt.plot(x, y2)

    # Define os limites dos eixos x e y
    plt.xlim(0, 1)
    plt.ylim(-0.3, 0.3)

    # Define os intervalos entre x e y
    plt.xticks(np.arange(0, 1.1, 0.1))
    plt.yticks(np.arange(-0.3, 0.4, 0.1))

    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.grid(True)
    return plt


def main():
    print("Bem-vindo ao primeiro trabalho computacional de MAP2220")
    print("Feito por Matheus Oliveira da Silva, nºUSP: 13696262")
    file_name = input("Insira o nome do arquivo que deseja ler: ")
    try:
        with open(file_name, 'r') as file:
            text = file.readlines()
    except FileNotFoundError:
        print(f"Arquivo '{file_name}' não encontrado.")
    except Exception as e:
        print(f"Um erro ocorreu: {e}")

    nome = text[0]
    numpontos = text[1].split()
    # Número de pontos em cada superficie
    np_extradorso = int(numpontos[0].split(".")[0])
    np_intradorso = int(numpontos[1].split(".")[0])
    print(f"Você abriu o arquivo de {nome}")
    print(f"Com {np_extradorso} pontos na superfície extradorsal")
    print(f"e {np_intradorso} pontos na superfície intradorsal")

    intradorso = []
    extradorso = []
    # A lista de coordenadas de cada modelo possui um elemento vazio que divide entre coordenadas do extradorso e
    # coordenadas do intradorso
    for i in range(np_extradorso):
        extradorso.append(text[i + 3].split())
    for j in range(4 + np_extradorso, len(text)):
        intradorso.append(text[j].split())
    thetas = np.linspace(0, math.pi, num=500)
    x = (1 / 2) * (1 - np.cos(thetas))
    coefextra = coeficientes(np_extradorso - 1, extradorso)
    coefintra = coeficientes(np_intradorso - 1, intradorso)
    yextra = interpola(coefextra, extradorso, x)
    yintra = interpola(coefintra, intradorso, x)


    indice_max, espessura_max = espessura(yextra, yintra)  # indice da maior espessura entre os 500 pontos
    j = 0
    while x[indice_max] > float(extradorso[j][0]):
        j += 1
    n_prec = input("Insira o número de novos pontos para aumentar a precisão (escolha um número menor que 1000 para "
                   "máquinas menos potentes): ")
    x_preciso = np.linspace(float(extradorso[j - 1][0]), float(extradorso[j + 1][0]), num=int(n_prec))
    yextra_preciso = []
    yintra_preciso = []
    for x_i in x_preciso:
        yextra_preciso.append(calculay(coefextra, j, extradorso, x_i))
        yintra_preciso.append(calculay(coefintra, j, intradorso, x_i))
    indce, novaespessura = espessura(yextra_preciso, yintra_preciso)
    print(espessura_max)
    print(x[indice_max])
    print(f"A maior espessura é de {novaespessura}, no ponto x = {x_preciso[indce]}")
    grafico = plota(x, yextra, yintra)
    grafico.show()

main()
