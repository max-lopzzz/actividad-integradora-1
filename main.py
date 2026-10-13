import os
from utils import leer_secuencia, leer_proteinas

BASE = os.path.dirname(os.path.abspath(__file__))
CARPETA = os.path.join(BASE, "Archivos")


def ruta(nombre):
    return os.path.join(CARPETA, nombre)


def usar(titulo, modulo, funcion, *args):
    print("\n" + "-" * 70)
    print(titulo)
    print("-" * 70)
    try:
        mod = __import__(modulo)
        return getattr(mod, funcion)(*args)
    except Exception as e:
        print("(no se pudo correr " + modulo + "." + funcion + ":", repr(e) + ")")
        return None


def main():
    # Lectura de archivos
    wuhan = leer_secuencia(ruta("SARS-COV-2-MN908947.3.txt"))
    texas = leer_secuencia(ruta("SARS-COV-2-MT106054.1.txt"))
    proteinas = leer_proteinas(ruta("seq-proteins.txt"))
    genes = {
        "M": ruta("gen-M.txt"),
        "S": ruta("gen-S.txt"),
        "ORF1AB": ruta("gen-ORF1AB.txt"),
    }

    usar("PUNTO 1: genes en el genoma", "punto1", "punto1", CARPETA)
    usar("PUNTO 2: palíndromos más largos", "punto2", "punto2", genes)
    usar("PUNTO 3: proteínas en el genoma", "punto3", "resolver",
         wuhan, proteinas)
    usar("PUNTO 4: Wuhan vs Texas", "punto4", "punto4", CARPETA)


if __name__ == "__main__":
    main()