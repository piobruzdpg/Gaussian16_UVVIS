import customtkinter as ctk
import tkinter as tk
from tkinter import filedialog, messagebox
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import re
import os

# --- Constants & Physics ---
H_PLANCK = 6.62607015e-34
C_LIGHT = 299792458
C_LIGHT_CM = C_LIGHT * 100
EV_TO_J = 1.602176634e-19
EV_TO_CM1 = 8065.544
NM_TO_CM = 1e-7
EPSILON_PREFACTOR = 1.3062974e8


# --- Math Logic (Ported from your script) ---

def parse_gaussian_log(filepath):
    """Parses Gaussian log file for Excited States."""
    excited_states = []

    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        return None, str(e)

    lines = content.splitlines()

    # Regex patterns
    excited_state_pattern = re.compile(
        r"^\s*Excited State\s+(\d+):\s+(.+?)\s+([-+]?\d*\.\d+)\s+eV\s+([-+]?\d*\.\d+)\s+nm\s+f=([-+]?\d*\.\d+)\s+<S\*\*2>=([-+]?\d*\.\d+)"
    )
    transition_pattern = re.compile(r"^\s*(\d+[AB]?)\s+(?:->|<-)\s+(\d+[AB]?)\s+([-+]?\d*\.\d+)")

    current_state_transitions = []
    current_state_data = {}

    for i, line in enumerate(lines):
        match_state = excited_state_pattern.match(line)
        if match_state:
            if current_state_data:
                current_state_data['Transitions'] = " | ".join(current_state_transitions)
                excited_states.append(current_state_data)

            state_num, symmetry, energy_ev, wavelength_nm, osc_strength, s2 = match_state.groups()
            current_state_data = {
                "State": int(state_num),
                "Symmetry": symmetry,
                "Energy (eV)": float(energy_ev),
                "Wavelength (nm)": float(wavelength_nm),
                "Oscillator Strength (f)": float(osc_strength),
                "<S**2>": float(s2),
                "Transitions": ""
            }
            current_state_transitions = []

            j = i + 1
            while j < len(lines):
                line_trans = lines[j].strip()
                match_trans = transition_pattern.match(line_trans)
                if match_trans:
                    orb_from, orb_to, contrib = match_trans.groups()
                    current_state_transitions.append(f"{orb_from} -> {orb_to} ({contrib})")
                elif not line_trans or not (line_trans[0].isdigit() or line_trans.startswith(' ')):
                    break
                elif "Excited State" in lines[j]:
                    break
                j += 1

    if current_state_data:
        current_state_data['Transitions'] = " | ".join(current_state_transitions)
        excited_states.append(current_state_data)

    if not excited_states:
        return None, "No excited states found."

    df = pd.DataFrame(excited_states)
    df['Energy (cm^-1)'] = df['Energy (eV)'] * EV_TO_CM1
    return df, None


def gaussian_broadening(nu_cm, nu_i_cm, f_i, sigma_cm):
    if sigma_cm <= 0: return np.zeros_like(nu_cm)
    exponent = -((nu_cm - nu_i_cm) / sigma_cm) ** 2
    exponent = np.clip(exponent, -700, 700)
    epsilon_i = EPSILON_PREFACTOR * (f_i / sigma_cm) * np.exp(exponent)
    return epsilon_i


def generate_spectrum(excitations_df, lambda_min, lambda_max, lambda_step, sigma_ev):
    if excitations_df is None or excitations_df.empty:
        return None

    sigma_cm = sigma_ev * EV_TO_CM1

    # Avoid zero division if lambda starts at 0 (unlikely but safe)
    start = max(0.1, lambda_min)
    wavelength_nm_plot = np.arange(start, lambda_max + lambda_step, lambda_step)
    wavenumber_cm_plot = 1.0 / (wavelength_nm_plot * NM_TO_CM)

    total_epsilon = np.zeros_like(wavenumber_cm_plot)

    for _, row in excitations_df.iterrows():
        nu_i_cm = row['Energy (cm^-1)']
        f_i = row['Oscillator Strength (f)']
        if f_i > 0:
            total_epsilon += gaussian_broadening(wavenumber_cm_plot, nu_i_cm, f_i, sigma_cm)

    spectrum_df = pd.DataFrame({
        'Wavelength (nm)': wavelength_nm_plot,
        'Molar Absorptivity': total_epsilon
    })
    return spectrum_df


# --- GUI Application ---

ctk.set_appearance_mode("Dark")  # Modes: "System" (standard), "Dark", "Light"
ctk.set_default_color_theme("blue")  # Themes: "blue" (standard), "green", "dark-blue"


class UVVisApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("UV-Vis Spectrum Generator (Gaussian)")
        self.geometry("1200x800")

        # Data placeholders
        self.df_raw = None
        self.df_spectrum = None
        self.current_filepath = None

        # --- Grid Layout Configuration ---
        self.grid_columnconfigure(0, weight=0)  # Left Panel (Fixed width)
        self.grid_columnconfigure(1, weight=1)  # Right Area (Expandable)
        self.grid_rowconfigure(0, weight=1)  # Full height

        # --- Left Panel (Controls) ---
        self.left_frame = ctk.CTkFrame(self, width=250, corner_radius=0)
        self.left_frame.grid(row=0, column=0, sticky="nsew", padx=0, pady=0)

        self.create_controls()

        # --- Right Panel (Plot + Table) ---
        self.right_frame = ctk.CTkFrame(self, corner_radius=0, fg_color="transparent")
        self.right_frame.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)
        self.right_frame.grid_rowconfigure(0, weight=3)  # Plot area dominates
        self.right_frame.grid_rowconfigure(1, weight=1)  # Table area smaller
        self.right_frame.grid_columnconfigure(0, weight=1)

        # Plot Area
        self.plot_frame = ctk.CTkFrame(self.right_frame, fg_color="white")  # White for publication look
        self.plot_frame.grid(row=0, column=0, sticky="nsew", pady=(0, 10))

        # Initialize Matplotlib Figure
        self.init_plot()

        # Table Area
        self.table_frame = ctk.CTkFrame(self.right_frame)
        self.table_frame.grid(row=1, column=0, sticky="nsew")

        self.table_label = ctk.CTkLabel(self.table_frame, text="Transition Table (Filtered)",
                                        font=ctk.CTkFont(size=14, weight="bold"))
        self.table_label.pack(pady=5)

        # ZMIANA: wrap="none" wyłącza zawijanie, font Courier
        self.table_text = ctk.CTkTextbox(self.table_frame, font=("Courier", 12), wrap="none")
        self.table_text.pack(expand=True, fill="both", padx=5, pady=5)
        self.table_text.configure(state="disabled")

    def create_controls(self):
        # Title
        lbl_title = ctk.CTkLabel(self.left_frame, text="Parameters", font=ctk.CTkFont(size=20, weight="bold"))
        lbl_title.grid(row=0, column=0, padx=20, pady=(20, 10))

        # File Loader
        self.btn_load = ctk.CTkButton(self.left_frame, text="Open Log File", command=self.load_file)
        self.btn_load.grid(row=1, column=0, padx=20, pady=10)

        # Parameters (Teraz podajemy nazwę zmiennej jawnie jako ostatni argument)
        self.add_param_input(2, "FWHM (σ) [eV]:", "0.4", "entry_fwhm")
        self.add_param_input(3, "Min Wavelength [nm]:", "200", "entry_min")
        self.add_param_input(4, "Max Wavelength [nm]:", "700", "entry_max")
        self.add_param_input(5, "Step [nm]:", "1.0", "entry_step")
        self.add_param_input(6, "Osc. Strength Threshold:", "0.005", "entry_osc")  # Tu był błąd

        # NOWE POLA - WYMIARY WYKRESU
        ctk.CTkLabel(self.left_frame, text="Figure Size", font=ctk.CTkFont(size=14, weight="bold")).grid(row=7,
                                                                                                         column=0,
                                                                                                         pady=(15, 5))
        self.add_param_input(8, "Width [cm]:", "15.0", "entry_width")
        self.add_param_input(9, "Height [cm]:", "10.0", "entry_height")

        # Update Button (przesuwamy go niżej)
        self.btn_update = ctk.CTkButton(self.left_frame, text="Update Plot & Table", command=self.update_view,
                                        fg_color="green")
        self.btn_update.grid(row=10, column=0, padx=20, pady=20)

        # Export Section (również przesuwamy niżej)
        ctk.CTkLabel(self.left_frame, text="Export", font=ctk.CTkFont(size=16, weight="bold")).grid(row=11, column=0,
                                                                                                    pady=(10, 5))

        self.btn_save_csv = ctk.CTkButton(self.left_frame, text="Save Data (CSV)", command=self.save_csv)
        self.btn_save_csv.grid(row=12, column=0, padx=20, pady=5)

        self.btn_save_fig = ctk.CTkButton(self.left_frame, text="Save Figure (600 dpi)", command=self.save_figure)
        self.btn_save_fig.grid(row=13, column=0, padx=20, pady=5)

    def add_param_input(self, row, label_text, default_val, attr_name):
        lbl = ctk.CTkLabel(self.left_frame, text=label_text, anchor="w")
        lbl.grid(row=row, column=0, padx=20, pady=(5, 0), sticky="w")

        entry = ctk.CTkEntry(self.left_frame)
        entry.insert(0, default_val)
        entry.grid(row=row + 100, column=0, padx=20, pady=(0, 10), sticky="ew")

        # Bezpośrednie przypisanie nazwy atrybutu
        setattr(self, attr_name, entry)

    def init_plot(self):
        # Setup clean matplotlib figure
        self.fig, self.ax1 = plt.subplots()
        self.ax2 = self.ax1.twinx()

        # Visual styling for publication
        self.fig.patch.set_facecolor('white')
        self.ax1.set_facecolor('white')

        self.setup_axes_style(self.ax1)
        self.setup_axes_style(self.ax2)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(expand=True, fill="both")

        # Add toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, self.plot_frame)
        toolbar.update()
        self.canvas.get_tk_widget().pack(expand=True, fill="both")

    def setup_axes_style(self, ax):
        ax.spines['bottom'].set_color('black')
        ax.spines['top'].set_color('black')
        ax.spines['right'].set_color('black')
        ax.spines['left'].set_color('black')
        ax.tick_params(axis='x', colors='black')
        ax.tick_params(axis='y', colors='black')
        ax.yaxis.label.set_color('black')
        ax.xaxis.label.set_color('black')
        ax.title.set_color('black')

    def load_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("Gaussian Log/Out", "*.log *.out"), ("All Files", "*.*")])
        if file_path:
            self.current_filepath = file_path
            df, error = parse_gaussian_log(file_path)
            if df is not None:
                self.df_raw = df
                messagebox.showinfo("Success", f"Loaded {len(df)} excited states.")
                self.update_view()
            else:
                messagebox.showerror("Error", error)

    def get_params(self):
        try:
            sigma = float(self.entry_fwhm.get())
            l_min = float(self.entry_min.get())
            l_max = float(self.entry_max.get())
            l_step = float(self.entry_step.get())
            thresh = float(self.entry_osc.get())

            # NOWE: Pobieranie wymiarów
            width_cm = float(self.entry_width.get())
            height_cm = float(self.entry_height.get())

            return sigma, l_min, l_max, l_step, thresh, width_cm, height_cm
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numeric values.")
            return None

    def update_view(self):
        if self.df_raw is None:
            return

        params = self.get_params()
        if not params: return
        # Rozpakowujemy teraz więcej zmiennych
        sigma, l_min, l_max, l_step, thresh, w_cm, h_cm = params

        # 1. Update Spectrum Data
        self.df_spectrum = generate_spectrum(self.df_raw, l_min, l_max, l_step, sigma)

        # 2. Update Plot
        self.ax1.clear()
        self.ax2.clear()

        # NOWE: Ustawienie rozmiaru figury (cm -> cale)
        # 1 cal = 2.54 cm
        self.fig.set_size_inches(w_cm / 2.54, h_cm / 2.54)

        # Plot Spectrum (Blue Line)
        self.ax1.plot(self.df_spectrum['Wavelength (nm)'],
                      self.df_spectrum['Molar Absorptivity'],
                      color='blue', label='Absorptivity', linewidth=1.5)

        # Plot Sticks (Red Vertical Lines)
        mask = (self.df_raw['Wavelength (nm)'] >= l_min) & (self.df_raw['Wavelength (nm)'] <= l_max)
        subset = self.df_raw[mask]

        if not subset.empty:
            self.ax2.vlines(x=subset['Wavelength (nm)'], ymin=0, ymax=subset['Oscillator Strength (f)'],
                            colors='red', linewidth=1.5, label='Osc. Strength')

        # Formatting
        self.ax1.set_xlabel("Wavelength (nm)", fontsize=12)
        self.ax1.set_ylabel(r"Molar Absorptivity ($\epsilon$) [L mol$^{-1}$ cm$^{-1}$]", fontsize=12, color='blue')

        self.ax2.set_ylabel("Oscillator Strength (f)", fontsize=12, color='red', rotation=-90, labelpad=20)
        self.ax2.yaxis.set_label_position("right")

        self.ax1.set_xlim(l_min, l_max)
        self.ax1.set_ylim(bottom=0)
        self.ax2.set_ylim(bottom=0)

        self.ax1.tick_params(axis='y', colors='blue')
        self.ax2.tick_params(axis='y', colors='red')

        title = os.path.basename(self.current_filepath) if self.current_filepath else "Spectrum"
        self.ax1.set_title(f"{title} (σ = {sigma} eV)", fontsize=14, color='black')

        # Layout update
        self.fig.tight_layout()
        self.canvas.draw()

        # 3. Update Table
        self.update_table(thresh)

    def update_table(self, threshold):
        if self.df_raw is None: return

        # Filtrowanie
        df_filtered = self.df_raw[self.df_raw['Oscillator Strength (f)'] >= threshold].copy()

        # Przygotowanie danych do wyświetlenia (zaokrąglanie)
        df_filtered['Energy (eV)'] = df_filtered['Energy (eV)'].round(4)
        df_filtered['Wavelength (nm)'] = df_filtered['Wavelength (nm)'].round(2)
        df_filtered['Oscillator Strength (f)'] = df_filtered['Oscillator Strength (f)'].round(4)
        # Zaokrąglamy S**2 (zazwyczaj 4 miejsca po przecinku wystarczą)
        df_filtered['<S**2>'] = df_filtered['<S**2>'].round(4)

        # Odblokowanie pola do edycji
        self.table_text.configure(state="normal")
        self.table_text.delete("1.0", "end")

        if df_filtered.empty:
            self.table_text.insert("end", f"No transitions found with f >= {threshold}")
        else:
            # ZMIANA: Dodano kolumnę <S**2> do nagłówka (szerokość 10 znaków)
            header = f"{'State':<8} {'E (eV)':<12} {'Lambda (nm)':<14} {'f':<12} {'<S**2>':<10} {'Transitions'}"

            self.table_text.insert("end", header + "\n")
            self.table_text.insert("end", "-" * 130 + "\n")  # Wydłużyłem nieco linię oddzielającą

            for _, row in df_filtered.iterrows():
                # ZMIANA: Dodano wartość <S**2> do wiersza
                line = (f"{row['State']:<8} "
                        f"{row['Energy (eV)']:<12} "
                        f"{row['Wavelength (nm)']:<14} "
                        f"{row['Oscillator Strength (f)']:<12} "
                        f"{row['<S**2>']:<10} "
                        f"{row['Transitions']}\n")
                self.table_text.insert("end", line)

        # Zablokowanie pola (tylko do odczytu)
        self.table_text.configure(state="disabled")

    def save_csv(self):
        if self.df_spectrum is None:
            messagebox.showwarning("Warning", "No spectrum generated yet.")
            return

        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV Files", "*.csv")])
        if path:
            self.df_spectrum.to_csv(path, index=False)
            messagebox.showinfo("Success", "Data saved successfully.")

    def save_figure(self):
        if self.df_spectrum is None:
            messagebox.showwarning("Warning", "No spectrum generated yet.")
            return

        path = filedialog.asksaveasfilename(defaultextension=".png",
                                            filetypes=[("PNG Image", "*.png"), ("PDF", "*.pdf"), ("SVG", "*.svg")])
        if path:
            self.fig.savefig(path, dpi=600, bbox_inches='tight')
            messagebox.showinfo("Success", "Figure saved with 600 DPI.")


if __name__ == "__main__":
    app = UVVisApp()
    app.mainloop()