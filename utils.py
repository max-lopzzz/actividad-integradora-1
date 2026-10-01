def leer_secuencia(ruta):
    #Lee un archivo y devuelve la secuencia como un solo string en mayúsculas sin saltos de línea
    with open(ruta, "r", encoding="utf-8") as f:
        lineas = f.read().splitlines()
    partes = []
    for linea in lineas:
        linea = linea.strip()
        if linea and not linea.startswith(">"):
            partes.append(linea)
    return "".join(partes).upper()


def leer_proteinas(ruta):
    #Devuelve una lista de (nombre, aminoacidos) leyendo el archivo FASTA
    proteinas = []
    nombre, partes = None, []
    with open(ruta, "r", encoding="utf-8") as f:
        for linea in f:
            linea = linea.strip()
            if not linea:
                continue
            if linea.startswith(">"):
                if nombre is not None:
                    proteinas.append((nombre, "".join(partes).upper()))
                nombre, partes = linea[1:], []
            else:
                partes.append(linea)
    if nombre is not None:
        proteinas.append((nombre, "".join(partes).upper()))
    return proteinas