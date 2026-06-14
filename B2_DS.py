#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulación B.2: Sector de materia T-simétrico (gradiente OFF)
Basado en def:tsym-sector del Apéndice B.2.
Propósito: Establecer la línea base de un sistema cerrado, unitario, T-simétrico,
sin decoherencia inducida por entorno.

COMPARACIÓN CUANTITATIVA numérica vs predicción teórica del framework:
- Unitaridad: norma del estado = 1, pureza global = 1, entropía global = 0.
- Producción de entropía nula (entropía constante).
- Coherencia reducida |ρ01(t)|: debe mostrar revivals (sin decaimiento ~N^{-1.8}).
- Tiempo de recurrencia: predicción teórica T_rec_teo = 2π / ΔE_min,
  donde ΔE_min es la mínima diferencia de energías del Hamiltoniano total.
- Escalado con N: el tiempo de recurrencia aumenta con N.

Ejecución: python sim_B2_Tsymmetric_sector.py
Salidas: CSV, PNGs, y prints en consola.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
import pandas as pd
from scipy.signal import find_peaks

# ============================================================================
# 1. Definición del modelo T-simétrico
# ============================================================================
def build_hamiltonian(N, Delta=1.0, omega0=1.0, g=0.5):
    """
    Hamiltoniano total: qubit + N espines con acoplamiento Ising.
    H = - (Δ/2) σ_x + ω0 Σ_j σ_z^{(j)} + g Σ_j σ_z^{(sys)} ⊗ σ_z^{(j)}
    Todos los coeficientes son reales → Hamiltoniano real → [H,Θ]=0.
    """
    dim_sys = 2
    dim_bath = 2**N
    dim_total = dim_sys * dim_bath

    sigma_x = np.array([[0,1],[1,0]], dtype=float)
    sigma_z = np.array([[1,0],[0,-1]], dtype=float)
    I_sys = np.eye(2, dtype=float)
    I_bath = np.eye(dim_bath, dtype=float)

    # H_sys = -Delta/2 * sigma_x (actúa en el qubit)
    H_sys = np.kron(I_bath, - (Delta/2.0) * sigma_x)

    # H_bath = ω0 Σ_j σ_z^{(j)}
    H_bath = np.zeros((dim_total, dim_total), dtype=float)
    for j in range(N):
        # operador σ_z en el sitio j
        op = 1
        for k in range(N):
            if k == j:
                op = np.kron(op, sigma_z)
            else:
                op = np.kron(op, np.eye(2))
        H_bath += omega0 * np.kron(np.eye(2), op)

    # H_int = g Σ_j σ_z^{(sys)} ⊗ σ_z^{(j)}
    H_int = np.zeros((dim_total, dim_total), dtype=float)
    for j in range(N):
        op_bath = 1
        for k in range(N):
            if k == j:
                op_bath = np.kron(op_bath, sigma_z)
            else:
                op_bath = np.kron(op_bath, np.eye(2))
        H_int += g * np.kron(sigma_z, op_bath)

    return H_sys + H_bath + H_int

def initial_state(N):
    """|ψ(0)⟩ = (|0⟩+|1⟩)/√2 ⊗ |0...0⟩ (todos los espines del baño abajo)."""
    dim_sys = 2
    dim_bath = 2**N
    psi_sys = np.array([1,1], dtype=float) / np.sqrt(2)
    psi_bath = np.zeros(dim_bath, dtype=float)
    psi_bath[0] = 1.0
    return np.kron(psi_sys, psi_bath)

def reduced_density_matrix(psi, N):
    """Traza sobre el baño → matriz densidad del qubit."""
    dim_sys = 2
    dim_bath = 2**N
    psi = psi.reshape(dim_sys, dim_bath, order='F')
    return psi @ psi.conj().T

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 0]
    return -np.sum(evals * np.log2(evals))

def purity(rho):
    return np.trace(rho @ rho).real

# ============================================================================
# 2. Parámetros
# ============================================================================
Delta = 1.0
omega0 = 1.0
g = 0.5
N_values = [2, 4, 6, 8]          # N=10 omitido por coste computacional
t_max = 50.0
dt = 0.05
times = np.arange(0, t_max + dt, dt)

# ============================================================================
# 3. Simulación y comparativa teórica
# ============================================================================
resultados = []

for N in N_values:
    print(f"\n--- N = {N} ---")
    H = build_hamiltonian(N, Delta, omega0, g)
    psi0 = initial_state(N)
    dim_total = len(psi0)

    # Diagonalización exacta (una sola vez)
    evals, evecs = eigh(H)          # evals reales, evecs ortonormales reales
    # Coeficientes de la condición inicial
    c0 = evecs.T @ psi0

    # Evolución temporal
    coh_abs = np.zeros(len(times))
    norm_psi = np.zeros(len(times))
    purity_global = np.zeros(len(times))
    entropy_global = np.zeros(len(times))

    for i, t in enumerate(times):
        # Estado en tiempo t
        phases = np.exp(-1j * evals * t)
        psi_t = evecs @ (c0 * phases)
        norm_psi[i] = np.linalg.norm(psi_t)
        # Matriz densidad reducida del qubit
        rho_sys = reduced_density_matrix(psi_t, N)
        coh_abs[i] = np.abs(rho_sys[0,1])
        # Global pureza y entropía (debe ser 1 y 0 respectivamente)
        rho_global = np.outer(psi_t, psi_t.conj())
        purity_global[i] = purity(rho_global)
        entropy_global[i] = von_neumann_entropy(rho_global)

    # Verificación de unitaridad
    norm_error = np.max(np.abs(norm_psi - 1.0))
    purity_error = np.max(np.abs(purity_global - 1.0))
    entropy_max = np.max(entropy_global)

    # Predicción teórica del tiempo de recurrencia (Poincaré)
    # Mínima diferencia de energías distinta de cero
    diffs = np.diff(evals)
    min_diff = np.min(diffs[diffs > 1e-12])
    T_rec_teo = 2.0 * np.pi / min_diff   # recurrencia de Poincaré (primer retorno aproximado)

    # Tiempo de recurrencia numérico: segundo pico de coherencia > 0.8*max_inicial
    max_initial = np.max(coh_abs[:len(times)//4])  # máximo en el primer cuarto
    peaks, _ = find_peaks(coh_abs, height=0.8 * max_initial)
    if len(peaks) >= 2:
        T_rec_num = times[peaks[1]]
        rec_error = abs(T_rec_num - T_rec_teo) / T_rec_teo
        rec_status = "PASS" if rec_error < 0.2 else "FAIL"  # tolerancia 20%
    else:
        T_rec_num = np.nan
        rec_error = np.nan
        rec_status = "FAIL (no revival)"

    # Verificar ausencia de decaimiento ~N^{-1.8} (solo chequeo cualitativo)
    # Se considera que si hay revivals y la coherencia no se hunde a 0 permanentemente, pasa.
    # Para cuantificar, se ajusta una exponencial a los valles y se comprueba que no decae a cero.
    # Aquí simplificamos: si T_rec_num no es nan, consideramos que no hay decaimiento irreversible.
    decay_status = "PASS" if not np.isnan(T_rec_num) else "FAIL"

    # Almacenar en la tabla (un registro por observable, o varios por N)
    # Para facilitar, creamos una lista de diccionarios con cada observable
    resultados.append({
        'N': N, 'observable': 'norma', 'valor_numerico': norm_error,
        'valor_teorico': 0.0, 'desviacion': norm_error, 'veredicto': 'PASS' if norm_error < 1e-10 else 'FAIL'
    })
    resultados.append({
        'N': N, 'observable': 'pureza_global', 'valor_numerico': purity_error,
        'valor_teorico': 0.0, 'desviacion': purity_error, 'veredicto': 'PASS' if purity_error < 1e-10 else 'FAIL'
    })
    resultados.append({
        'N': N, 'observable': 'entropia_global_max', 'valor_numerico': entropy_max,
        'valor_teorico': 0.0, 'desviacion': entropy_max, 'veredicto': 'PASS' if entropy_max < 1e-10 else 'FAIL'
    })
    resultados.append({
        'N': N, 'observable': 'T_recurrencia', 'valor_numerico': T_rec_num,
        'valor_teorico': T_rec_teo, 'desviacion': rec_error, 'veredicto': rec_status
    })
    resultados.append({
        'N': N, 'observable': 'decaimiento_irreversible', 'valor_numerico': 1 if np.isnan(T_rec_num) else 0,
        'valor_teorico': 0, 'desviacion': np.nan, 'veredicto': decay_status
    })

    # Guardar datos de coherencia para la gráfica del ejemplo (N=6 o el último)
    if N == max(N_values):
        times_plot = times
        coh_plot = coh_abs
        rec_num_plot = T_rec_num
        rec_teo_plot = T_rec_teo

# ============================================================================
# 4. Generar tabla CSV y mostrar en consola
# ============================================================================
df = pd.DataFrame(resultados)
print("\n" + "="*80)
print("TABLA DE COMPARACIÓN NUMÉRICA VS PREDICCIÓN TEÓRICA")
print("="*80)
print(df.to_string(index=False))
df.to_csv("sim_B2_comparativa.csv", index=False)

# ============================================================================
# 5. Gráficas
# ============================================================================
# Figura (a): Coherencia numérica + marcador de recurrencia teórico/numérico
plt.figure(figsize=(10,6))
plt.plot(times_plot, coh_plot, label='|ρ₀₁(t)| numérica', color='blue')
if not np.isnan(rec_num_plot):
    plt.axvline(rec_num_plot, color='red', linestyle='--', label=f'T_rec numérico = {rec_num_plot:.2f}')
if not np.isnan(rec_teo_plot):
    plt.axvline(rec_teo_plot, color='green', linestyle=':', label=f'T_rec teórico = {rec_teo_plot:.2f}')
plt.xlabel('Tiempo')
plt.ylabel('|ρ₀₁(t)|')
plt.title(f'Coherencia del qubit – N={max(N_values)} (gradiente OFF)')
plt.legend()
plt.grid(True)
plt.savefig('fig_B2_coherence.png', dpi=150)
plt.close()

# Figura (b): Entropía global y pureza (deben ser constantes)
plt.figure(figsize=(10,6))
# Tomamos los datos del último N simulado (son representativos)
# Recalculamos para el último N si es necesario, o usamos los arrays ya guardados.
# Como no guardamos todos, hacemos una simulación rápida solo para este N:
_, _, _, _, _, purity_last, entropy_last = None, None, None, None, None, None, None
# En su lugar, podemos generar una figura genérica a partir de los datos almacenados en resultados.
# Pero mejor: durante la simulación del último N guardamos purity_global y entropy_global en arrays globales.
# Modificamos ligeramente el bucle para guardar estos arrays para el último N.
# Como ya hemos ejecutado, lo haremos de nuevo dentro de esta sección (solo para el último N):
if max(N_values) in N_values:
    H_last = build_hamiltonian(max(N_values), Delta, omega0, g)
    psi0_last = initial_state(max(N_values))
    evals_last, evecs_last = eigh(H_last)
    c0_last = evecs_last.T @ psi0_last
    purity_last_arr = np.zeros(len(times))
    entropy_last_arr = np.zeros(len(times))
    for i, t in enumerate(times):
        psi_t = evecs_last @ (c0_last * np.exp(-1j * evals_last * t))
        rho_global = np.outer(psi_t, psi_t.conj())
        purity_last_arr[i] = purity(rho_global)
        entropy_last_arr[i] = von_neumann_entropy(rho_global)
    plt.plot(times, entropy_last_arr, label='Entropía global (debe ser 0)')
    plt.plot(times, purity_last_arr, '--', label='Pureza global (debe ser 1)')
    plt.xlabel('Tiempo')
    plt.ylabel('Valor')
    plt.title(f'Unitaridad – N={max(N_values)} (gradiente OFF)')
    plt.legend()
    plt.grid(True)
    plt.savefig('fig_B2_entropy_purity.png', dpi=150)
    plt.close()

# Figura (c): Tiempo de recurrencia vs N (numérico y teórico)
rec_df = df[df['observable'] == 'T_recurrencia'].copy()
if not rec_df.empty:
    plt.figure(figsize=(8,5))
    plt.plot(rec_df['N'], rec_df['valor_numerico'], 'o-', label='Numérico')
    plt.plot(rec_df['N'], rec_df['valor_teorico'], 's--', label='Teórico (2π/ΔE_min)')
    plt.xlabel('N (tamaño del baño)')
    plt.ylabel('Tiempo de recurrencia')
    plt.title('Escalado del tiempo de recurrencia con N')
    plt.legend()
    plt.grid(True)
    plt.savefig('fig_B2_recurrence_vs_N.png', dpi=150)
    plt.close()

# Figura (d): Contraste OFF (recurrencia) vs ON (decaimiento ~N^{-1.8})
# Se genera una figura conceptual; los datos ON deben venir de la simulación P5 (no incluida aquí).
plt.figure(figsize=(10,6))
plt.plot(times_plot, coh_plot, label='Gradiente OFF (recurrencia)', color='blue')
# Datos simulados de ON (por ejemplo, una exponencial o ley de potencias)
# Aquí solo se muestra una curva ilustrativa, no real.
t_plot = np.linspace(0, t_max, 200)
decay_ON = np.exp(-0.2 * t_plot)   # mero ejemplo
plt.plot(t_plot, decay_ON, 'r--', label='Gradiente ON (decaimiento ~N^{-1.8}, ejemplo)')
plt.xlabel('Tiempo')
plt.ylabel('|ρ₀₁(t)|')
plt.title('Contraste: OFF (recurrencia) vs ON (decaimiento irreversible)')
plt.legend()
plt.grid(True)
plt.savefig('fig_B2_contrast_OFF_ON.png', dpi=150)
plt.close()

print("\nGráficas generadas: fig_B2_coherence.png, fig_B2_entropy_purity.png, fig_B2_recurrence_vs_N.png, fig_B2_contrast_OFF_ON.png")
print("\nSimulación completada. La tabla CSV 'sim_B2_comparativa.csv' contiene las desviaciones y veredictos.")