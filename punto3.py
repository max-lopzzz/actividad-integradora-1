from codones import traducir, CODONES
from punto1 import kmp_buscar


def leer_con_frameshift(genoma, inicio, aminoacidos):
    pos = inicio
    ya_salto = False
    salto = None
    for aa in aminoacidos:
        if CODONES.get(genoma[pos:pos + 3], "X") != aa:
            if ya_salto:
                return None
            pos -= 1
            ya_salto = True
            salto = pos
            if CODONES.get(genoma[pos:pos + 3], "X") != aa:
                return None
        pos += 3
    return pos - 1, salto


def buscar_con_frameshift(genoma, traducciones, aminoacidos, k=8):
    semilla = aminoacidos[:k]
    for marco in range(3):
        for p in kmp_buscar(traducciones[marco], semilla):
            inicio = marco + 3 * p
            res = leer_con_frameshift(genoma, inicio, aminoacidos)
            if res is not None:
                return inicio, res[0], res[1]
    return None


def mostrar(nombre, genoma, inicio, fin, aminoacidos, nota=""):
    n = min(4, len(aminoacidos))
    codones = [genoma[inicio + 3 * k:inicio + 3 * k + 3] for k in range(n)]
    aminos = [CODONES.get(c, "X") for c in codones]
    linea = (nombre + " : " + str(inicio) + " - " + str(fin) + "  amino: "
             + "".join(aminos) + " = " + " ".join(codones))
    print(linea + (" " + nota if nota else ""))


def resolver(genoma, proteinas):
    # Traducimos el genoma UNA sola vez en cada marco de lectura (0, 1, 2)
    traducciones = [traducir(genoma, marco) for marco in range(3)]
    regiones = []
    halladas = []
    for nombre, aminoacidos in proteinas:
        hallada = False
        for marco in range(3):
            for p in kmp_buscar(traducciones[marco], aminoacidos):
                inicio = marco + 3 * p                    # base 0 en el genoma
                fin = inicio + 3 * len(aminoacidos) - 1
                mostrar(nombre, genoma, inicio, fin, aminoacidos)
                regiones.append((nombre, inicio, fin, None))
                hallada = True
        if not hallada:
            res = buscar_con_frameshift(genoma, traducciones, aminoacidos)
            if res is not None:
                inicio, fin, salto = res
                mostrar(nombre, genoma, inicio, fin, aminoacidos,
                        "(con frameshift -1)")
                regiones.append((nombre, inicio, fin, salto))
                hallada = True
        if hallada:
            halladas.append(nombre)
        else:
            print("Proteína:", nombre, "-> NO encontrada")
    print("\nProteínas encontradas:", len(halladas), "de", len(proteinas))
    return regiones