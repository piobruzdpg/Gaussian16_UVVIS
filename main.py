import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import re
import io
import base64  # Do generowania linków do pobrania plików

# --- Stałe fizyczne i konwersje ---
H_PLANCK = 6.62607015e-34  # J*s
C_LIGHT = 299792458  # m/s
C_LIGHT_CM = C_LIGHT * 100  # cm/s
EV_TO_J = 1.602176634e-19  # J/eV
EV_TO_CM1 = 8065.544  # cm^-1 / eV (1 / (H_PLANCK * C_LIGHT_CM / EV_TO_J))
NM_TO_CM = 1e-7  # cm/nm
AVOGADRO = 6.02214076e23  # mol^-1
LN_10 = np.log(10)
# Stała z równania 5 (przeliczona dla wygody - uwzględnia stałe fizyczne i jednostki)
# Oryginał: 1.3062974e8 * f / sigma_cm1 * exp(...) -> epsilon w L mol^-1 cm^-1
EPSILON_PREFACTOR = 1.3062974e8


# --- Funkcje pomocnicze ---

def parse_gaussian_log(uploaded_file):
    """
    Parsuje plik log Gaussian w poszukiwaniu danych o stanach wzbudzonych (TD-DFT/EOM).
    """
    excited_states = []
    # Upewniamy się, że czytamy od początku pliku
    uploaded_file.seek(0)
    content = uploaded_file.read().decode("utf-8")
    lines = content.splitlines()

    # Wzór na główną linię stanu wzbudzonego
    # Zmiana z ([\w-]+) na (.+?) oraz doprecyzowanie liczb [-+]?\d*\.\d+
    excited_state_pattern = re.compile(
        r"^\s*Excited State\s+(\d+):\s+(.+?)\s+([-+]?\d*\.\d+)\s+eV\s+([-+]?\d*\.\d+)\s+nm\s+f=([-+]?\d*\.\d+)\s+<S\*\*2>=([-+]?\d*\.\d+)"
    )

    # Wzór na linie przejść orbitalnych
    # Zmiana z -> na (?:->|<-). (?:...) to grupa niechwytająca
    transition_pattern = re.compile(r"^\s*(\d+[AB]?)\s+(?:->|<-)\s+(\d+[AB]?)\s+([-+]?\d*\.\d+)")

    current_state_transitions = []
    current_state_data = {}

    for i, line in enumerate(lines):
        match_state = excited_state_pattern.match(line)
        if match_state:
            # Zapisz poprzedni stan (jeśli istniał)
            if current_state_data:
                current_state_data['Transitions'] = " | ".join(current_state_transitions)
                excited_states.append(current_state_data)

            # Rozpocznij nowy stan
            state_num, symmetry, energy_ev, wavelength_nm, osc_strength, s2 = match_state.groups()
            current_state_data = {
                "State": int(state_num),
                "Symmetry": symmetry,
                "Energy (eV)": float(energy_ev),
                "Wavelength (nm)": float(wavelength_nm),
                "Oscillator Strength (f)": float(osc_strength),
                "<S**2>": float(s2),
                "Transitions": ""  # Placeholder
            }
            current_state_transitions = []

            # Sprawdź kolejne linie w poszukiwaniu przejść dla tego stanu
            j = i + 1
            while j < len(lines):
                line_trans = lines[j].strip()
                # Sprawdź, czy linia pasuje do wzorca przejścia
                match_trans = transition_pattern.match(line_trans)
                if match_trans:
                    orb_from, orb_to, contrib = match_trans.groups()
                    # Formatuj opis przejścia
                    current_state_transitions.append(f"{orb_from} -> {orb_to} ({contrib})")
                # Jeśli linia jest pusta lub zaczyna się od czegoś innego niż liczba/spacja, zakończ szukanie przejść
                elif not line_trans or not (line_trans[0].isdigit() or line_trans.startswith(' ')):
                    break
                # Heurystyka: Jeśli napotkamy kolejny 'Excited State', to też kończymy
                elif "Excited State" in lines[j]:
                    break
                j += 1

        # Jeśli jesteśmy na końcu pliku, dodaj ostatni znaleziony stan
        if i == len(lines) - 1 and current_state_data:
            current_state_data['Transitions'] = " | ".join(current_state_transitions)
            excited_states.append(current_state_data)

    if not excited_states:
        return None  # Nie znaleziono stanów wzbudzonych

    df = pd.DataFrame(excited_states)
    # Oblicz energię w cm^-1
    df['Energy (cm^-1)'] = df['Energy (eV)'] * EV_TO_CM1
    return df


def gaussian_broadening(nu_cm, nu_i_cm, f_i, sigma_cm):
    """
    Oblicza wkład pojedynczego przejścia do widma (Równanie 5 z dokumentu).
    """
    if sigma_cm <= 0:  # Unikaj dzielenia przez zero
        return np.zeros_like(nu_cm)
    exponent = -((nu_cm - nu_i_cm) / sigma_cm) ** 2
    # Zapobiegaj underflow/overflow dla dużych/małych wykładników
    exponent = np.clip(exponent, -700, 700)
    epsilon_i = EPSILON_PREFACTOR * (f_i / sigma_cm) * np.exp(exponent)
    return epsilon_i


def generate_spectrum(excitations_df, lambda_min, lambda_max, lambda_step, sigma_ev):
    """
    Generuje widmo UV-Vis sumując broadened pasma Gaussa.
    """
    if excitations_df is None or excitations_df.empty:
        return None, None

    # Konwertuj sigma z eV na cm^-1
    sigma_cm = sigma_ev * EV_TO_CM1

    # Utwórz siatkę długości fal (nm) i liczb falowych (cm^-1)
    wavelength_nm_plot = np.arange(lambda_min, lambda_max + lambda_step, lambda_step)
    # Unikaj dzielenia przez zero dla lambda = 0
    wavelength_nm_plot = wavelength_nm_plot[wavelength_nm_plot > 0]
    wavenumber_cm_plot = 1.0 / (wavelength_nm_plot * NM_TO_CM)

    # Inicjalizuj całkowite epsilon
    total_epsilon = np.zeros_like(wavenumber_cm_plot)

    # Iteruj po każdym stanie wzbudzonym i dodaj jego wkład
    for _, row in excitations_df.iterrows():
        nu_i_cm = row['Energy (cm^-1)']
        f_i = row['Oscillator Strength (f)']
        if f_i > 0:
            total_epsilon += gaussian_broadening(wavenumber_cm_plot, nu_i_cm, f_i, sigma_cm)

    spectrum_df = pd.DataFrame({
        'Wavelength (nm)': wavelength_nm_plot,
        'Wavenumber (cm^-1)': wavenumber_cm_plot,
        'Molar Absorptivity (L mol^-1 cm^-1)': total_epsilon
    })

    return spectrum_df, sigma_cm


def create_download_link(df, filename, link_text):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">{link_text}</a>'
    return href


def create_figure_download_link(fig, filename, format, link_text):
    img_bytes = fig.to_image(format=format)
    b64 = base64.b64encode(img_bytes).decode()
    href = f'<a href="data:image/{format};base64,{b64}" download="{filename}">{link_text}</a>'
    return href


# --- Interfejs Streamlit ---
st.set_page_config(layout="wide")
st.title("Generator Widm UV-Vis z plików log Gaussian")

st.markdown("""
Ta aplikacja wczytuje plik log z obliczeń stanów wzbudzonych programu Gaussian (np. TD-DFT, EOM-CCSD)
i generuje symulowane widmo UV-Vis, stosując rozmycie Gaussowskie dla każdego przejścia.
""")

# --- Sidebar z opcjami ---
st.sidebar.header("Opcje Generowania Widma")
uploaded_file = st.sidebar.file_uploader("1. Wgraj plik log Gaussian (.log lub .out)", type=["log", "out"])

default_sigma_ev = 0.4
default_lambda_min = 200.0
default_lambda_max = 700.0
default_lambda_step = 1.0
default_f_threshold = 0.01

sigma_ev = st.sidebar.number_input(
    "2. Szerokość połówkowa Gaussa σ (eV)",
    min_value=0.01,
    value=default_sigma_ev,
    step=0.05,
    format="%.2f",
    help=f"Parametr σ kontroluje szerokość pików w widmie. Wartość domyślna w GaussView to {default_sigma_ev} eV."
)

lambda_min = st.sidebar.number_input("3. Minimalna długość fali (nm)", min_value=10.0, value=default_lambda_min,
                                     step=10.0)
lambda_max = st.sidebar.number_input("4. Maksymalna długość fali (nm)", min_value=lambda_min + 10.0,
                                     value=default_lambda_max, step=10.0)
lambda_step = st.sidebar.number_input("5. Krok długości fali (nm)", min_value=0.1, max_value=10.0,
                                      value=default_lambda_step, step=0.1)

f_threshold = st.sidebar.slider(
    "6. Próg siły oscylatora (f) do wyświetlenia",
    min_value=0.0,
    max_value=1.0,
    value=default_f_threshold,
    step=0.001,
    format="%.4f",
    help="Filtruje listę wyświetlanych przejść."
)

# --- Główna część aplikacji ---
if uploaded_file is not None:
    try:
        excitations_df_raw = parse_gaussian_log(uploaded_file)

        if excitations_df_raw is None or excitations_df_raw.empty:
            st.error("Nie znaleziono danych o stanach wzbudzonych w podanym pliku.")
        else:
            st.success(f"Znaleziono {len(excitations_df_raw)} stanów wzbudzonych w pliku '{uploaded_file.name}'.")

            # --- Wyświetlanie danych o wzbudzeniach ---
            st.subheader("Dane o Obliczonych Stanach Wzbudzonych")
            st.markdown(f"Wyświetlane są stany z siłą oscylatora (f) ≥ {f_threshold:.4f}")

            excitations_df_filtered = excitations_df_raw[
                excitations_df_raw['Oscillator Strength (f)'] >= f_threshold].copy()
            excitations_df_display = excitations_df_filtered[[
                "State", "Symmetry", "Energy (eV)", "Wavelength (nm)", "Oscillator Strength (f)", "Transitions",
                "<S**2>"
            ]].round(4)

            st.dataframe(excitations_df_display)

            # --- Generowanie i wyświetlanie widma ---
            st.subheader("Wygenerowane Widmo UV-Vis")
            spectrum_df, sigma_cm_used = generate_spectrum(
                excitations_df_raw,
                lambda_min,
                lambda_max,
                lambda_step,
                sigma_ev
            )

            if spectrum_df is not None:
                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=spectrum_df['Wavelength (nm)'],
                    y=spectrum_df['Molar Absorptivity (L mol^-1 cm^-1)'],
                    mode='lines',
                    name=f'Widmo (σ={sigma_ev} eV)'
                ))

                # --- POPRAWKA: Filtrowanie markerów i ustawianie zakresu ---

                # Wybieramy tylko te piki, które mieszczą się w wybranym zakresie
                peaks_df = excitations_df_raw[
                    (excitations_df_raw['Oscillator Strength (f)'] > 0) &
                    (excitations_df_raw['Wavelength (nm)'] >= lambda_min) &
                    (excitations_df_raw['Wavelength (nm)'] <= lambda_max)
                    ]

                # Obliczanie pozycji Y dla markerów
                peak_epsilon_values = []
                for _, peak in peaks_df.iterrows():
                    idx = (np.abs(spectrum_df['Wavelength (nm)'] - peak['Wavelength (nm)'])).argmin()
                    if idx < len(spectrum_df):
                        peak_epsilon_values.append(spectrum_df['Molar Absorptivity (L mol^-1 cm^-1)'].iloc[idx])
                    else:
                        peak_epsilon_values.append(0)

                fig.add_trace(go.Scatter(
                    x=peaks_df['Wavelength (nm)'],
                    y=peak_epsilon_values,
                    mode='markers',
                    marker=dict(symbol='x', size=8, color='red'),
                    name='Położenie obliczonych przejść (eV)',
                    hoverinfo='text',
                    hovertext=[
                        f"State: {row['State']}<br>λ: {row['Wavelength (nm)']:.2f} nm<br>E: {row['Energy (eV)']:.3f} eV<br>f: {row['Oscillator Strength (f)']:.4f}"
                        for _, row in peaks_df.iterrows()
                    ]
                ))

                fig.update_layout(
                    title=f"Symulowane widmo UV-Vis (σ = {sigma_ev:.2f} eV ≈ {sigma_cm_used:.1f} cm⁻¹)",
                    xaxis_title="Długość fali λ (nm)",
                    yaxis_title="Molowa absorpcja ε (L mol⁻¹ cm⁻¹)",
                    hovermode="x unified",
                    xaxis=dict(range=[lambda_min, lambda_max])  # Wymuszenie zakresu osi X
                )
                # --- KONIEC POPRAWKI ---

                st.plotly_chart(fig, use_container_width=True)

                # --- Opcje pobierania ---
                st.subheader("Zapisz wyniki")
                col1, col2 = st.columns(2)

                with col1:
                    st.markdown("Pobierz dane widma:")
                    csv_filename = f"spectrum_{uploaded_file.name.split('.')[0]}_sigma{sigma_ev:.2f}eV.csv"
                    st.markdown(create_download_link(spectrum_df, csv_filename, 'Pobierz jako CSV'),
                                unsafe_allow_html=True)

                with col2:
                    st.markdown("Pobierz figurę widma:")
                    png_filename = f"spectrum_{uploaded_file.name.split('.')[0]}_sigma{sigma_ev:.2f}eV.png"
                    svg_filename = f"spectrum_{uploaded_file.name.split('.')[0]}_sigma{sigma_ev:.2f}eV.svg"
                    try:
                        st.markdown(create_figure_download_link(fig, png_filename, "png", 'Pobierz jako PNG'),
                                    unsafe_allow_html=True)
                        st.markdown(create_figure_download_link(fig, svg_filename, "svg", 'Pobierz jako SVG'),
                                    unsafe_allow_html=True)
                    except Exception as e:
                        st.warning(f"Nie można wygenerować linków do pobrania figury. Błąd: {e}")

            else:
                st.warning("Nie można wygenerować widma. Sprawdź dane wejściowe.")

    except Exception as e:
        st.error(f"Wystąpił błąd podczas przetwarzania pliku: {e}")
        st.exception(e)