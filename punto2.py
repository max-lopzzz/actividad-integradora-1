# Mancher propio, longitud de palindromo mas largo por gen guardado en un archivo de salida

from utils import leer_secuencia


def manacher(s):
    """Devuelve (inicio, longitud) del palíndromo más largo de s en O(n).

    Se inserta un separador '#' entre cada carácter (y en los extremos)
    para tratar igual los palíndromos de longitud par e impar.
    Ej: "ABBA" -> "#A#B#B#A#"
    """
    t = "#" + "#".join(s) + "#"
    n = len(t)
    p = [0] * n          # p[i] = radio del palíndromo centrado en t[i]
    centro = 0           # centro del palíndromo más a la derecha conocido
    derecha = 0          # borde derecho de ese palíndromo

    for i in range(n):
        if i < derecha:
            # Reutiliza lo ya calculado usando el espejo de i respecto al centro
            p[i] = min(derecha - i, p[2 * centro - i])
        # Intenta expandir alrededor de i
        while (i + p[i] + 1 < n and i - p[i] - 1 >= 0
               and t[i + p[i] + 1] == t[i - p[i] - 1]):
            p[i] += 1
        # Actualiza el palíndromo más a la derecha
        if i + p[i] > derecha:
            centro = i
            derecha = i + p[i]

    # El mayor radio en t es igual a la longitud del palíndromo en s
    radio_max = max(p)
    centro_max = p.index(radio_max)
    inicio = (centro_max - radio_max) // 2
    return inicio, radio_max


def punto2(genes, archivo_salida="palindromos.txt"):
    """genes: dict {nombre: ruta}. Muestra y guarda el palíndromo más largo de cada gen."""
    lineas = []
    for nombre, ruta in genes.items():
        secuencia = leer_secuencia(ruta)
        inicio, longitud = manacher(secuencia)
        palindromo = secuencia[inicio:inicio + longitud]
        lineas.append(
            f"Gen {nombre}: longitud del palíndromo más largo = {longitud} "
            f"(índices {inicio}-{inicio + longitud - 1}) -> {palindromo}"
        )
    with open(archivo_salida, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")
    for l in lineas:
        print(l)


if __name__ == "__main__":
    punto2({
        "M": "Archivos/gen-M.txt",
        "S": "Archivos/gen-S.txt",
        "ORF1AB": "Archivos/gen-ORF1AB.txt",
    })