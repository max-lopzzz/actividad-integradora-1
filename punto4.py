import os
from utils import leer_secuencia
from codones import CODONES


ARCHIVO_WUHAN = "SARS-COV-2-MN908947.3.txt"
ARCHIVO_TEXAS = "SARS-COV-2-MT106054.1.txt"

REGIONES = [
    ("ORF1ab", 265, 21555),
    ("S", 21562, 25384),
    ("ORF3a", 25392, 26220),
    ("E", 26244, 26472),
    ("M", 26522, 27191),
    ("ORF6", 27201, 27387),
    ("ORF7a", 27393, 27759),
    ("ORF7b", 27755, 27887),
    ("ORF8", 27893, 28259),
    ("N", 28273, 29533),
    ("ORF10", 29557, 29674),
]


def comparar_genomas(wuhan, texas):
    # Compara posicion por posicion, solo hasta donde llega el mas corto
    longitud = min(len(wuhan), len(texas))
    diferencias = []

    racha_actual = 0
    racha_maxima = 0
    inicio_actual = 0
    inicio_maximo = -1

    for i in range(longitud):
        if wuhan[i] == texas[i]:
            # Si son iguales, la racha sigue creciendo
            if racha_actual == 0:
                inicio_actual = i

            racha_actual += 1

            if racha_actual > racha_maxima:
                racha_maxima = racha_actual
                inicio_maximo = inicio_actual
        else:
            # Si son distintas, guardamos la diferencia y la racha empieza de nuevo
            diferencias.append((i, wuhan[i], texas[i]))
            racha_actual = 0

    if racha_maxima == 0:
        return diferencias, -1, -1, 0

    fin_maximo = inicio_maximo + racha_maxima - 1
    return diferencias, inicio_maximo, fin_maximo, racha_maxima


def traducir_codon(codon):
    # Busca el aminoacido en la tabla del equipo (X si no se encuentra)
    return CODONES.get(codon, "X")


def analizar_codon(indice, wuhan, texas):
    # Busca en que gen cae la diferencia y saca su codon completo
    for gen, inicio, fin in REGIONES:
        if inicio <= indice < fin:
            inicio_marco = inicio
            codones_antes = 0

            # ORF1ab tiene un salto de marco en 13467, ahi se vuelve a contar de tres en tres
            # El primer tramo tiene 13203 letras, o sea 4401 codones
            if gen == "ORF1ab" and indice >= 13467:
                inicio_marco = 13467
                codones_antes = 4401

            # El codon empieza en la ultima posicion que cae en el marco de lectura
            inicio_codon = inicio_marco + 3 * ((indice - inicio_marco) // 3)

            # Si el codon se sale del gen, no se analiza
            if inicio_codon + 3 > fin:
                break

            num_codon = codones_antes + (inicio_codon - inicio_marco) // 3 + 1
            codon_wuhan = wuhan[inicio_codon:inicio_codon + 3]
            codon_texas = texas[inicio_codon:inicio_codon + 3]
            aa_wuhan = traducir_codon(codon_wuhan)
            aa_texas = traducir_codon(codon_texas)

            if aa_wuhan != aa_texas:
                cambio = "Si"
            else:
                cambio = "No"

            return (gen, num_codon, codon_wuhan, codon_texas,
                    aa_wuhan, aa_texas, cambio)

    # Si no cae en ningun gen, no se calcula el aminoacido
    return None


def guardar_racha(wuhan, inicio, fin, nombre_archivo="racha_mas_larga.txt"):
    # Guarda la cadena mas larga sin diferencias en un archivo
    if inicio == -1:
        return

    with open(nombre_archivo, "w", encoding="utf-8") as archivo:
        archivo.write(">Racha sin diferencias, indices " + str(inicio) + " a " + str(fin) + "\n")

        for i in range(inicio, fin + 1, 70):
            archivo.write(wuhan[i:min(i + 70, fin + 1)] + "\n")


def punto4(carpeta="Archivos"):
    # Lee los dos genomas y hace la comparacion
    wuhan = leer_secuencia(os.path.join(carpeta, ARCHIVO_WUHAN))
    texas = leer_secuencia(os.path.join(carpeta, ARCHIVO_TEXAS))

    diferencias, inicio, fin, longitud_racha = comparar_genomas(wuhan, texas)
    longitud_comparada = min(len(wuhan), len(texas))
    iguales = longitud_comparada - len(diferencias)

    # Armamos el resumen para mostrarlo y guardarlo
    lineas = [
        "Longitud Wuhan: " + str(len(wuhan)),
        "Longitud Texas: " + str(len(texas)),
        "Posiciones comparadas: " + str(longitud_comparada),
        "Posiciones iguales: " + str(iguales),
        "Diferencias encontradas: " + str(len(diferencias)),
        "",
        "Indice (desde 0) | Wuhan | Texas"
    ]

    for indice, letra_wuhan, letra_texas in diferencias:
        lineas.append(str(indice) + " (" + str(indice + 1) + ") | " + letra_wuhan + " | " + letra_texas)

    lineas += [
        "",
        "Racha mas larga sin diferencias:",
        "Inicio: " + str(inicio),
        "Fin (incluido): " + str(fin),
        "Longitud: " + str(longitud_racha)
    ]

    # Si un genoma es mas largo, las letras que sobran se cuentan aparte
    if len(wuhan) > longitud_comparada:
        lineas.append("Letras extra al final de Wuhan: " + str(len(wuhan) - longitud_comparada))
        lineas.append("Contenido extra: " + wuhan[longitud_comparada:])
    elif len(texas) > longitud_comparada:
        lineas.append("Letras extra al final de Texas: " + str(len(texas) - longitud_comparada))
        lineas.append("Contenido extra: " + texas[longitud_comparada:])

    # Codones y aminoacidos de cada diferencia
    lineas += ["", "Analisis de codones:"]

    for indice, letra_wuhan, letra_texas in diferencias:
        lineas.append("")
        lineas.append("Indice " + str(indice) + " (" + str(indice + 1) + "): " + letra_wuhan + " -> " + letra_texas)
        resultado = analizar_codon(indice, wuhan, texas)

        if resultado is None:
            lineas.append("No cae en ningun gen conocido")
        else:
            gen, num, codon_w, codon_t, aa_w, aa_t, cambio = resultado
            lineas.append("Gen: " + gen + " (codon " + str(num) + ")")
            lineas.append("Codon Wuhan: " + codon_w + "  Codon Texas: " + codon_t)
            lineas.append("Aminoacido Wuhan: " + aa_w + "  Aminoacido Texas: " + aa_t)
            lineas.append("Cambia el aminoacido: " + cambio)

    resumen = "\n".join(lineas)
    print(resumen)

    with open("resultados_punto4.txt", "w", encoding="utf-8") as archivo:
        archivo.write(resumen + "\n")

    guardar_racha(wuhan, inicio, fin)
    return diferencias, inicio, fin, longitud_racha


if __name__ == "__main__":
    if os.path.isdir("Archivos"):
        punto4("Archivos")
    else:
        punto4(".")
