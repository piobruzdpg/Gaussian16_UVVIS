# Gaussian UV-Vis Spectrum Generator

Narzędzie typu desktop GUI (napisane w Pythonie) służące do analizy i wizualizacji widm elektronowych (UV-Vis) na podstawie plików wynikowych programu Gaussian (TD-DFT).

Aplikacja została stworzona z myślą o pracy naukowej – pozwala na szybką weryfikację przejść, kontrolę kontaminacji spinu oraz generowanie wykresów gotowych do publikacji.

## Główne funkcjonalności

* **Wczytywanie danych**: Obsługa standardowych plików logów Gaussiana (`.log`, `.out`).
* **Wizualizacja**:
    * Interaktywny wykres widma (splot gaussowski) oraz przejść (sticks).
    * Estetyka "publication-ready" (białe tło, wyraźne osie).
    * Możliwość definiowania szerokości połówkowej (FWHM/$\sigma$).
* **Analiza danych**:
    * Tabela przejść z filtrowaniem po sile oscylatora ($f$).
    * **Weryfikacja spinu**: Automatyczne porównanie obliczonego $<S^2>$ z wartością idealną i oznaczanie stanów z dużą kontaminacją spinu.
* **Eksport**:
    * **Excel (.xlsx)**: Zapis danych do jednego pliku z podziałem na arkusze: widmo (punkty), surowe przejścia oraz przejścia przefiltrowane.
    * **Grafika**: Zapis wykresów do formatów PNG, PDF lub SVG (domyślnie 600 DPI).

## Wymagania

Skrypt wymaga Pythona 3 oraz następujących bibliotek:

* `customtkinter` (GUI)
* `pandas` & `numpy` (analiza danych)
* `matplotlib` (wykresy)
* `openpyxl` (obsługa zapisu do Excela)

## Instalacja i uruchomienie

1.  Zainstaluj wymagane biblioteki:
    ```bash
    pip install customtkinter pandas matplotlib numpy openpyxl
    ```
2.  Uruchom skrypt:
    ```bash
    python main_ctk.py
    ```
