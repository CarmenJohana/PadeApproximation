import csv
import matplotlib.pyplot as plt
from fractions import Fraction

archivos = ["datos_prob3.csv", "datos_prob4.csv"]

datos = {}

for archivo in archivos:
    n = int(archivo.split("datos_prob")[1].split(".csv")[0])
    qs = []
    probs = []
    with open(archivo, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            qs.append(int(row["q"]))
            # Convertir fracción a float
            probs.append(float(Fraction(row["prob"])))
    datos[n] = (qs, probs)

# --- Graficar ---
plt.figure(figsize=(8,5))

for n, (qs, probs) in sorted(datos.items()):
    plt.plot(qs, probs, marker='o', label=f"n={n}")

plt.xlabel("q (número primo)")
plt.ylabel("Probabilidad")
plt.title("Probabilidad de recuperación del polinomio minimal")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
