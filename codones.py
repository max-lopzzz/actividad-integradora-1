BASES = "TCAG"
AMINOS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONES = {}
_i = 0
for a in BASES:
    for b in BASES:
        for c in BASES:
            CODONES[a + b + c] = AMINOS[_i]
            _i += 1
def traducir(secuencia, marco=0):
    aminoacidos = []
    for i in range(marco, len(secuencia) - 2, 3):
        aminoacidos.append(CODONES.get(secuencia[i:i + 3], "X"))
    return "".join(aminoacidos)