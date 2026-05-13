import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os

NAZWA_PLIKU = 'orbity_gwiazd.txt'

if not os.path.exists(NAZWA_PLIKU):
    print(f"Błąd: Nie znaleziono pliku '{NAZWA_PLIKU}'.")
    exit()

print("Wczytywanie danych (to może chwilę potrwać przy dużej liczbie gwiazd)...")
# Wczytujemy całą macierz. Każdy wiersz to jeden krok czasowy.
data = np.loadtxt(NAZWA_PLIKU)

# Obliczamy liczbę gwiazd na podstawie liczby kolumn (2 kolumny = 1 gwiazda)
LICZBA_KLATEK, KOLUMNY = data.shape
LICZBA_GWIAZD = KOLUMNY // 2
print(f"Wczytano {LICZBA_KLATEK} klatek dla {LICZBA_GWIAZD} gwiazd.")

# --- USTAWIENIA WYKRESU ---
fig, ax = plt.subplots(figsize=(9, 9), facecolor='black')
ax.set_facecolor('black')
ax.set_xlim(0, 301)
ax.set_ylim(0, 301)

# Stylizacja w kosmicznym klimacie
ax.tick_params(colors='gray')
for spine in ax.spines.values():
    spine.set_edgecolor('gray')
ax.grid(True, linestyle=':', alpha=0.2, color='gray')

# Rysujemy Jądro Galaktyki (Czarną Dziurę) na środku
ax.plot(150.5, 150.5, marker='*', color='white', markersize=12)

# Tworzymy obiekt "scatter" - na początku pusty. 
# s=1 to rozmiar punktu (gwiazdy). alpha=0.6 daje efekt świecenia/przezroczystości
scatter = ax.scatter([], [], s=1.5, color='cyan', alpha=0.6)

# Parametr przyspieszenia (jeśli animacja jest za wolna, zwiększ np. na 2, 3, 5)
SKOK = 1 

# --- FUNKCJA AKTUALIZUJĄCA KLATKI ---
def update(frame):
    real_frame = min(frame * SKOK, LICZBA_KLATEK - 1)
    
    # Magia Numpy: wyciągamy wszystkie X (kolumny parzyste) i Y (kolumny nieparzyste)
    # z bieżącego wiersza w ułamku sekundy
    x = data[real_frame, 0::2]
    y = data[real_frame, 1::2]
    
    # Aktualizujemy pozycje wszystkich 1500 gwiazd jednym poleceniem
    scatter.set_offsets(np.c_[x, y])
    
    # Aktualizujemy tytuł okna, żeby pokazywał upływ czasu (krok symulacji)
    ax.set_title(f"Ewolucja Galaktyki | Krok: {real_frame * 5}", color='white', fontsize=14)
    
    return scatter,

# --- URUCHOMIENIE ANIMACJI ---
print("Generowanie animacji...")
# interval=16 daje ok. 60 FPS (1000 ms / 16 ms)
ani = animation.FuncAnimation(
    fig, update, frames=LICZBA_KLATEK // SKOK, interval=16, blit=False, repeat=False
)

plt.show()