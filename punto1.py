import os
import sys

ARCHIVO_GENOMA = "SARS-COV-2-MN908947.3.txt"
GENES = [
    ("M", "gen-M.txt"),
    ("S", "gen-S.txt"),
    ("ORF1AB", "gen-ORF1AB.txt"),
]


def leer_secuencia(ruta):
    """
    Lee un archivo de secuencia y regresa solo las letras en mayúscula.
    Ignora encabezados FASTA (líneas que empiezan con '>'), saltos de línea,
    espacios, '\\r' y números.
    Complejidad: O(L), L = tamaño del archivo.
    """
    partes = []
    with open(ruta, "r") as archivo:
        for linea in archivo:
            if linea.startswith(">"):
                continue
            partes.append("".join(c for c in linea.upper() if "A" <= c <= "Z"))
    return "".join(partes)


def calcular_lps(patron):
    """
    Tabla LPS de KMP: lps[i] = longitud del prefijo propio más largo de
    patron[0..i] que también es sufijo de patron[0..i].
    Complejidad: O(m) tiempo y espacio, m = len(patron).
    """
    m = len(patron)
    lps = [0] * m
    longitud = 0  # longitud del prefijo-sufijo actual
    i = 1
    while i < m:
        if patron[i] == patron[longitud]:
            longitud += 1
            lps[i] = longitud
            i += 1
        elif longitud != 0:
            longitud = lps[longitud - 1]  # retrocede con lo ya calculado, sin avanzar i
        else:
            lps[i] = 0
            i += 1
    return lps


def kmp_buscar(texto, patron):
    """
    Regresa una lista con TODOS los índices (base 0) donde inicia `patron`
    dentro de `texto`, incluyendo apariciones traslapadas.
    Nunca retrocede en el texto: ante un fallo usa la tabla LPS.
    Complejidad: O(n + m) tiempo, O(m) espacio.
    """
    n, m = len(texto), len(patron)
    if m == 0 or m > n:
        return []

    lps = calcular_lps(patron)
    apariciones = []
    i = j = 0  # i recorre el texto, j el patrón
    while i < n:
        if texto[i] == patron[j]:
            i += 1
            j += 1
            if j == m:  # coincidencia completa
                apariciones.append(i - m)
                j = lps[j - 1]
        elif j != 0:
            j = lps[j - 1]
        else:
            i += 1
    return apariciones


def punto1(carpeta="."):
    """
    Abre el genoma de Wuhan y los tres genes, busca cada gen con KMP e
    imprime nombre, índices de aparición y primeros 12 caracteres.
    Regresa un diccionario {nombre_gen: [(inicio, fin), ...]} por si
    otro punto lo necesita.
    """
    genoma = leer_secuencia(os.path.join(carpeta, ARCHIVO_GENOMA))

    print("=" * 70)
    print(" PUNTO 1: Apariciones de los genes en SARS-COV-2-MN908947.3 (Wuhan)")
    print("=" * 70)
    print(f"Longitud del genoma: {len(genoma)} nucleotidos")
    print("(Indices en base 0; entre parentesis la posicion en base 1, como en NCBI)\n")

    resultados = {}
    for nombre, archivo in GENES:
        ruta = os.path.join(carpeta, archivo)
        print(f"Gen: {nombre}  ({archivo})")
        try:
            gen = leer_secuencia(ruta)
        except OSError:
            print("  ERROR: no se pudo abrir el archivo.\n")
            continue

        print(f"  Longitud del gen:     {len(gen)} nucleotidos")
        print(f"  Primeros 12 chars:    {gen[:12]}")

        indices = kmp_buscar(genoma, gen)
        resultados[nombre] = [(ini, ini + len(gen) - 1) for ini in indices]

        if not indices:
            print("  Indices en el genoma: no se encontro en el genoma.\n")
            continue

        print(f"  Apariciones:          {len(indices)}")
        for ini, fin in resultados[nombre]:
            print(f"  Indices en el genoma: {ini} - {fin}   ({ini + 1} - {fin + 1})")
        print()

    return resultados


if __name__ == "__main__":
    punto1(sys.argv[1] if len(sys.argv) > 1 else ".")
