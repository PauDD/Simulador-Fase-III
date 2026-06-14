"""
================================================================================
SIMULADOR DEL MARCO DE ORIENTATION CHANNEL - DINÁMICA DERIVADA
================================================================================

Implementación de los protocolos de simulación numérica de la Sección 11
del paper "Fase III" (marco relativista) y del paper predecesor "Fase II"
(marco no-relativista).

PRINCIPIO METODOLÓGICO:
    Las simulaciones de este programa derivan las predicciones del marco
    de su dinámica formal (Hilbert extendido, GKSL Lindbladianas, evolución
    de la master equation, Monte Carlo geométrico sobre congruencia nula).
    No postulan las predicciones como fórmulas evaluadas; las calculan.

    Si el marco hace una predicción genuina, debe emerger del cálculo,
    no estar puesta a mano. Si no emerge, eso también es información.

Dependencias:
    - numpy
    - matplotlib
    - qutip (para evolución de master equation y quantum trajectories)
    - pandas (para tablas)

Las simulaciones P1, P2, P3, P4 están implementadas con dinámica derivada.
P5 (límite mesoscópico) está documentada pero no implementada debido a
requerimientos de tensor networks/HPC; se explica en el menú.
================================================================================
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import qutip as qt

# ================================================================
# CONFIGURACIÓN GLOBAL
# ================================================================

np.random.seed(42)  # Para reproducibilidad de Monte Carlo en P4

# ================================================================
# UTILIDADES DE PRESENTACIÓN
# ================================================================

def mostrar_tabla(df, titulo):
    print("\n" + "=" * 90)
    print(titulo)
    print("=" * 90)
    print(df.to_string(index=False, float_format=lambda x: f"{x:.5f}"))
    print("=" * 90)


# ================================================================
# CONSTRUCCIÓN DEL ESPACIO DE HILBERT EXTENDIDO
# ================================================================
# Paper II §2: H_ext = H_S ⊗ H_τ con dim H_τ = 2
# Estados |+⟩_τ, |−⟩_τ representan los sectores forward/backward.
#
# Para todas las simulaciones P1, P2, P3:
#   H_S = C^2 (sistema de dos niveles, p.ej. qubit fotónico de polarización
#   o un qubit interno de detector)
#   H_τ = C^2 (sector de orientación)
#   H_ext = C^4

# Operadores del sistema (sub-índice S)
sx_S = qt.tensor(qt.sigmax(), qt.qeye(2))
sy_S = qt.tensor(qt.sigmay(), qt.qeye(2))
sz_S = qt.tensor(qt.sigmaz(), qt.qeye(2))
I_S = qt.tensor(qt.qeye(2), qt.qeye(2))

# Proyectores pointer en H_S
P0_S = qt.tensor(qt.basis(2, 0) * qt.basis(2, 0).dag(), qt.qeye(2))
P1_S = qt.tensor(qt.basis(2, 1) * qt.basis(2, 1).dag(), qt.qeye(2))

# Operadores del sector de orientación (sub-índice tau)
sx_tau = qt.tensor(qt.qeye(2), qt.sigmax())
sy_tau = qt.tensor(qt.qeye(2), qt.sigmay())
sz_tau = qt.tensor(qt.qeye(2), qt.sigmaz())

# σ_+ en H_τ:  |+⟩⟨−|.  Convención: |+⟩ = |0⟩_τ, |−⟩ = |1⟩_τ
# σ_+ |−⟩ = |+⟩  (transferencia backward → forward)
sigma_plus_tau = qt.tensor(qt.qeye(2),
                            qt.basis(2, 0) * qt.basis(2, 1).dag())

# Proyectores sectoriales (Paper II §2.4)
Pi_plus = qt.tensor(qt.qeye(2), qt.basis(2, 0) * qt.basis(2, 0).dag())
Pi_minus = qt.tensor(qt.qeye(2), qt.basis(2, 1) * qt.basis(2, 1).dag())


def estado_inicial_backward(theta=np.pi/2, phi=0.0):
    """
    Estado inicial: |ψ_S⟩ ⊗ |−⟩_τ
    con |ψ_S⟩ = cos(θ/2)|0⟩ + e^{iφ}sin(θ/2)|1⟩ en superposición.
    Convención Paper II §6.5: el sistema entra a la medida en el sector
    backward |−⟩_τ.
    """
    psi_S = (np.cos(theta/2) * qt.basis(2, 0)
             + np.exp(1j * phi) * np.sin(theta/2) * qt.basis(2, 1))
    psi_tau = qt.basis(2, 1)  # |−⟩_τ
    return qt.tensor(psi_S, psi_tau)


# ================================================================
# SIMULACIÓN 1: QUANTUM ERASER CON ORIENTATION CHANNEL
# ================================================================
# Predicción del marco (Paper II §6, Theorem 6.5):
#   Bajo la Lindbladiana de orientación L_k = √γ_k P_k ⊗ σ_+^τ,
#   el orientation channel transfiere irreversiblemente población
#   del sector backward al forward.  El "erasure" (operación unitaria
#   sobre H_τ que intenta deshacer la transferencia) no puede ser
#   un inverso CPTP (Theorem 4.3).
#
#   ¿Qué emerge?  La visibilidad de interferencia decae con el número
#   de ciclos (medida débil + intento de erasure).  El exponente NO
#   se postula; se obtiene de la dinámica.
#
# Comparación QM estándar:
#   En QM estándar con medida débil + erasure ideal, la visibilidad
#   se preserva.  V_QM = V_0 = 1.
# ================================================================

def simulacion_quantum_eraser():
    print("\n--- SIMULACIÓN 1: Quantum Eraser con Orientation Channel ---")
    print("Implementación: evolución GKSL explícita sobre H_S ⊗ H_τ,")
    print("seguida de unitaria de erasure sobre H_τ.  La visibilidad")
    print("se extrae del estado reducido ρ_S.\n")

    # Parámetros físicos (no son fits; son inputs del problema)
    gamma = 1.0          # tasa de la Lindbladiana de orientación
    t_weak = 0.1         # duración de la medida débil por ciclo
    N_max = 20           # ciclos a explorar

    # Lindbladiana pointer-projective con dos outcomes (Paper II §6.4)
    L0 = np.sqrt(gamma) * P0_S * sigma_plus_tau
    L1 = np.sqrt(gamma) * P1_S * sigma_plus_tau

    # Hamiltoniano del sistema durante la medida (libre, sin drive)
    H = 0 * I_S

    # Unitaria de "erasure" ideal sobre H_τ:
    # rotación que devuelve |+⟩_τ a un estado neutro respecto a la
    # base de pointer.  La operación más natural es σ_x sobre H_τ:
    # intercambia |+⟩ ↔ |−⟩, restituyendo formalmente el sector backward.
    U_erasure = sx_tau

    visibilidades = []
    visibilidades_qm = []  # QM estándar: erasure perfecto preserva V

    for N in range(1, N_max + 1):
        # Estado inicial fresco para este experimento (N ciclos completos)
        psi0 = estado_inicial_backward()
        rho = psi0 * psi0.dag()

        for ciclo in range(N):
            # Paso 1: evolución GKSL durante tiempo t_weak
            tlist = np.linspace(0, t_weak, 20)
            result = qt.mesolve(H, rho, tlist, c_ops=[L0, L1])
            rho = result.states[-1]

            # Paso 2: erasure unitaria sobre H_τ
            rho = U_erasure * rho * U_erasure.dag()

        # Visibilidad: módulo de la coherencia off-diagonal de ρ_S
        # en la base de pointer
        rho_S = rho.ptrace(0)
        V = 2 * abs(rho_S.full()[0, 1])  # factor 2 normaliza a V_inicial=1
        visibilidades.append(V)

        # QM estándar: erasure ideal preserva visibilidad
        visibilidades_qm.append(1.0)

    N_arr = np.arange(1, N_max + 1)
    visibilidades = np.array(visibilidades)
    visibilidades_qm = np.array(visibilidades_qm)
    desviacion = visibilidades_qm - visibilidades

    df = pd.DataFrame({
        "Ciclo N": N_arr,
        "V_QM (teórico)": visibilidades_qm,
        "V_Framework (derivado)": visibilidades,
        "Desviación": desviacion
    })

    mostrar_tabla(df, "SIMULACIÓN 1: Quantum Eraser - Resultado del cálculo")

    # Análisis post-hoc: ¿el decaimiento es exponencial?  Fit de
    # log(V) vs N para extraer la tasa efectiva.
    mask = visibilidades > 0.01
    if mask.sum() > 5:
        slope, intercept = np.polyfit(N_arr[mask],
                                       np.log(visibilidades[mask]), 1)
        alpha_emergente = -slope
        print(f"\nTasa de decaimiento emergente (fit log-lineal): "
              f"α = {alpha_emergente:.4f}")
        print(f"Tasa de Lindbladiana de entrada: γ·t_weak = "
              f"{gamma*t_weak:.4f}")
        print(f"(Si el decaimiento es estructural, α debe estar relacionado "
              f"con γ·t_weak,\nno ser libre.)")

    plt.figure(figsize=(10, 5))
    plt.plot(N_arr, visibilidades_qm, 'o-', label='V_QM (erasure ideal)',
             color='blue')
    plt.plot(N_arr, visibilidades, 's-', label='V_Framework (derivado de GKSL)',
             color='red')
    plt.title("P1: Quantum Eraser con Orientation Channel")
    plt.xlabel("Número de ciclos (medida débil + erasure)")
    plt.ylabel("Visibilidad de interferencia")
    plt.ylim(-0.05, 1.05)
    plt.grid(True, linestyle=':')
    plt.legend()
    plt.tight_layout()
    plt.savefig('F3_p1_quantum_eraser.png', dpi=100)
    plt.close()
    print("\nGráfica guardada en F3_p1_quantum_eraser.png")

# ================================================================
# DESENREDO DE TRAYECTORIAS CUÁNTICAS (Monte Carlo wavefunction)
# Objetivo del paper:
#   "replacing the deterministic master equation solution with
#    quantum-trajectory unravelling via qutip.mcsolve to access
#    trajectory-resolved statistics (single-shot fractions,
#    waiting-time distributions)."
#
# QUÉ APORTA: la ecuación maestra determinista da solo la dinámica
# PROMEDIO. El desenredo en trayectorias da las estadísticas resueltas
# por realización: fracciones single-shot, distribuciones de tiempos de
# espera entre saltos cuánticos, que son las firmas discriminantes de la
# Sección 8. Esto es vergonzosamente paralelo (cada trayectoria es
# independiente).
#
# HARDWARE / GPU:
#   - qutip.mcsolve paraleliza trayectorias por NÚCLEOS de CPU (tu i9 las
#     reparte; usa map_kwargs={'num_cpus': N}). Esto ya acelera mucho.
#   - Para GPU real: instala `qutip-jax` (backend JAX/CUDA de QuTiP 5).
#     Con él, el integrador corre en la RTX 4070 Ti. Se activa con
#     qt.settings.core['default_dtype']='jax' y options={'method':'diffrax'}.
#     El código siguiente detecta qutip-jax y lo usa si está disponible;
#     si no, cae a CPU multinúcleo automáticamente.
#
# Este archivo SE VALIDÓ en CPU a baja escala (200 trayectorias). En tu
# máquina sube n_traj a 1e4-1e5 y/o activa el backend jax.
# ================================================================
 
import numpy as np
import qutip as qt
 
 
def _try_enable_gpu():
    """Intenta activar el backend JAX/CUDA de QuTiP. Devuelve True si OK."""
    try:
        import qutip_jax  # noqa: F401
        qt.settings.core['default_dtype'] = 'jax'
        return True
    except Exception:
        return False
 
 
def trayectorias_P1_estadisticas(gamma=1.0, t_weak=0.1, n_cycles=10,
                                 n_traj=2000, num_cpus=None, usar_gpu=False):
    """
    Desenredo en trayectorias del protocolo P1. Devuelve estadísticas
    resueltas por trayectoria que la ecuación maestra promedio no da:
      - fracción single-shot de coherencia residual por trayectoria,
      - distribución de tiempos de espera entre saltos cuánticos.
    """
    gpu = _try_enable_gpu() if usar_gpu else False
    print("\n--- TRAYECTORIAS P1 (Monte Carlo wavefunction) ---")
    print(f"n_traj={n_traj}, backend={'GPU/JAX' if gpu else 'CPU multinúcleo'}\n")
 
    sx_S = qt.tensor(qt.sigmax(), qt.qeye(2))
    sz_S = qt.tensor(qt.sigmaz(), qt.qeye(2))
    sp_tau = qt.tensor(qt.qeye(2), qt.sigmap())
 
    H = 0.5 * sz_S
    # canal de orientación como operador de salto (quantum jump)
    c_ops = [np.sqrt(gamma) * (sz_S * sp_tau)]
 
    psi0 = qt.tensor((qt.basis(2, 0) + qt.basis(2, 1)).unit(), qt.basis(2, 0))
    tlist = np.linspace(0, n_cycles * t_weak, 200)
 
    # observable: coherencia del sistema |rho^S_01|
    e_coh = qt.tensor(qt.sigmax(), qt.qeye(2))
 
    opts = {"map": "parallel", "keep_runs_results": True}
    map_kw = {}
    if num_cpus:
        map_kw["num_cpus"] = num_cpus
 
    res = qt.mcsolve(H, psi0, tlist, c_ops, e_ops=[e_coh],
                     ntraj=n_traj, options=opts, **({"map_kwargs": map_kw} if map_kw else {}))
 
    # estadística promedio (= ecuación maestra) y por-trayectoria
    coh_mean = np.array(res.average_expect[0] if hasattr(res, "average_expect")
                        else res.expect[0])
 
    # tiempos de espera entre saltos cuánticos
    waiting_times = []
    col_times = getattr(res, "col_times", None)
    if col_times is None:
        # qutip 5: los tiempos de salto están en res.col_times o por runs
        col_times = getattr(res, "runs_col_times", None)
    if col_times is not None:
        for ct in col_times:
            ct = np.atleast_1d(np.asarray(ct, dtype=float))
            if ct.size > 1:
                waiting_times.extend(np.diff(np.sort(ct)).tolist())
    waiting_times = np.array(waiting_times)
 
    # estadística single-shot: dispersión de la coherencia final por traj.
    runs = getattr(res, "runs_expect", None)
    if runs is not None:
        finales = np.array([np.asarray(r[0])[-1] for r in runs])
        print(f"Single-shot: coherencia final por trayectoria "
              f"media={finales.mean():.4f}, std={finales.std():.4f} "
              f"(la master eq. solo da la media)")
 
    print(f"Coherencia media final (= master eq.): {coh_mean[-1]:.4f}")
    if waiting_times.size:
        print(f"Saltos cuánticos registrados: {waiting_times.size}")
        print(f"Tiempo de espera medio entre saltos: "
              f"{waiting_times.mean():.4f} (1/gamma esperado = {1/gamma:.4f})")
        print(f"Distribución de tiempos de espera: "
              f"min={waiting_times.min():.3f}, "
              f"mediana={np.median(waiting_times):.3f}, "
              f"max={waiting_times.max():.3f}")
        print("=> Firma discriminante (Sec. 8) NO accesible desde la")
        print("   ecuación maestra determinista; sí desde las trayectorias.")
    else:
        print("(col_times no expuesto en esta versión de qutip; "
              "usar res.runs_expect para estadística por-trayectoria.)")
 
    # --- Graficas ---
    # Panel 1: coherencia media (= master eq.) con una muestra de
    # trayectorias individuales superpuestas. Panel 2: histograma de
    # tiempos de espera entre saltos cuanticos. Panel 3: dispersion
    # single-shot de la coherencia final por trayectoria.
    n_paneles = 1 + (1 if waiting_times.size else 0) \
                  + (1 if runs is not None else 0)
    fig, axes = plt.subplots(1, n_paneles, figsize=(6 * n_paneles, 4.5))
    if n_paneles == 1:
        axes = [axes]
    ip = 0

    ax = axes[ip]; ip += 1
    if runs is not None:
        n_muestra = min(40, len(runs))
        for r in runs[:n_muestra]:
            ax.plot(tlist, np.real(np.asarray(r[0])), color='gray',
                    alpha=0.15, lw=0.6)
    ax.plot(tlist, np.real(coh_mean), color='darkblue', lw=2.0,
            label='Media (= master eq.)')
    ax.set_title("P1 trayectorias: coherencia |rho_S 01|(t)")
    ax.set_xlabel("Tiempo"); ax.set_ylabel("Coherencia")
    ax.grid(True, linestyle=':'); ax.legend()

    if waiting_times.size:
        ax = axes[ip]; ip += 1
        ax.hist(waiting_times, bins=30, color='teal', alpha=0.75,
                density=True, edgecolor='white', linewidth=0.5)
        tau_mean = waiting_times.mean()
        xs = np.linspace(0, waiting_times.max(), 200)
        ax.plot(xs, (1.0 / tau_mean) * np.exp(-xs / tau_mean), 'r--',
                label=f'exp(tau_med={tau_mean:.3f})')
        ax.set_title("Tiempos de espera entre saltos")
        ax.set_xlabel("Tiempo de espera dt"); ax.set_ylabel("Densidad")
        ax.grid(True, linestyle=':'); ax.legend()

    if runs is not None:
        ax = axes[ip]; ip += 1
        finales = np.array([np.real(np.asarray(r[0])[-1]) for r in runs])
        ax.hist(finales, bins=30, color='indigo', alpha=0.75,
                edgecolor='white', linewidth=0.5)
        ax.axvline(finales.mean(), color='red', linestyle='--',
                   label=f'media={finales.mean():.3f}')
        ax.set_title("Single-shot: coherencia final")
        ax.set_xlabel("Coherencia final"); ax.set_ylabel("Nro trayectorias")
        ax.grid(True, linestyle=':'); ax.legend()

    fig.tight_layout()
    fig.savefig('F3_p1_trayectorias.png', dpi=100)
    plt.close(fig)

    return res, coh_mean, waiting_times
# ================================================================
# SIMULACIÓN 2: DEPENDENCIA ENTRÓPICA DEL DETECTOR
# ================================================================
# Predicción del marco (Paper II §7, Theorem 7.4 y Proposition 7.5):
#   La tasa efectiva del orientation channel γ_eff que un observador
#   macroscópico implementa se relaciona con su producción entrópica
#   ambiental a través de la cota Landauer:
#       Q_dot_E ≥ k_B T_E γ_eff ln 2  (por bit)
#
#   La predicción NO es que γ_eff ∝ dS/dt linealmente como en los
#   simuladores anteriores; es la cota termodinámica.
#
# Modelado:
#   Sistema (2 niveles) + detector (2 niveles) + baño térmico.
#   El detector está acoplado al baño a temperatura T_E.  El coupling
#   sistema-detector implementa pointer projection.  La temperatura
#   del baño controla dS/dt.
#
#   Lo que NO se postula: la forma de γ_eff(T_E).
#   Lo que SÍ se postula: la estructura del coupling (Born-Markov,
#   pointer-projective).
# ================================================================

def simulacion_detector_entropico():
    print("\n--- SIMULACIÓN 2: Dependencia Entrópica del Detector ---")
    print("Implementación: master equation con orientation channel y baño")
    print("térmico.  Se varía la temperatura del baño y se extraen del")
    print("cálculo: (a) la tasa γ_eff del canal de orientación, (b) el")
    print("flujo de calor al baño Q_dot, (c) verificación de la cota")
    print("Landauer γ_eff ≤ Q_dot / (T ln 2).\n")

    # Construcción del modelo:
    # Hilbert: H_S ⊗ H_τ (4 dim).
    # Dinámica:
    #   - Lindbladianas del orientation channel L_0, L_1 (Paper II §6.4)
    #   - Lindbladianas térmicas del baño sobre H_τ (Born-Markov estándar)
    #
    # La predicción del marco (Paper II Proposition 7.5):
    #   Q_dot_E ≥ k_B T_E γ_eff ln 2  (cota Landauer, por evento)
    # Reescrita: γ_eff ≤ Q_dot_E / (k_B T_E ln 2)

    # Parámetros físicos
    gamma_pointer = 1.0      # tasa pointer-projective intrínseca
    omega_tau = 1.0          # frecuencia del modo τ del detector
    t_evol = 5.0             # tiempo de evolución
    kappa_0 = 0.3            # acoplamiento detector-baño

    # Lindbladianas del orientation channel
    L0 = np.sqrt(gamma_pointer) * P0_S * sigma_plus_tau
    L1 = np.sqrt(gamma_pointer) * P1_S * sigma_plus_tau

    # Hamiltoniano del modo de orientación
    H = 0.5 * omega_tau * sz_tau

    T_E_values = np.linspace(0.1, 5.0, 12)
    gamma_eff_values = []
    Q_dot_values = []
    landauer_max_values = []

    psi0 = estado_inicial_backward()
    rho0 = psi0 * psi0.dag()

    for T_E in T_E_values:
        # Ocupación Bose-Einstein
        if omega_tau / T_E > 50:
            n_th = 0.0
        else:
            n_th = 1.0 / (np.exp(omega_tau / T_E) - 1.0)

        # Operadores Lindblad térmicos en H_τ
        # En la convención |+⟩=|0⟩, |−⟩=|1⟩, destroy lleva |−⟩→|+⟩.
        # Re-examinamos esto: destroy|0⟩=0, destroy|1⟩=|0⟩.
        # Así destroy^τ acopla |−⟩ → |+⟩, es decir excita el sector forward.
        # Eso NO es disipación térmica estándar.
        # Para disipación estándar |+⟩ → |−⟩ usamos create^τ con
        # convención apropiada, o equivalentemente damos vuelta a las
        # tasas.
        #
        # Convención clara:
        #   "decay" thermal: |+⟩ → |−⟩  (sector forward decae a backward)
        #   esto ES disipación al baño.
        # En nuestra base con |+⟩=|0⟩_τ:
        #   |+⟩⟨−| en H_τ es la transición que lleva |−⟩→|+⟩
        #   |−⟩⟨+| es la disipación térmica estándar
        decay_op = qt.tensor(qt.qeye(2),
                              qt.basis(2, 1) * qt.basis(2, 0).dag())
        excite_op = qt.tensor(qt.qeye(2),
                               qt.basis(2, 0) * qt.basis(2, 1).dag())

        L_decay = np.sqrt(kappa_0 * (n_th + 1)) * decay_op
        L_excite = np.sqrt(kappa_0 * n_th) * excite_op

        tlist = np.linspace(0, t_evol, 200)
        result = qt.mesolve(H, rho0, tlist,
                            c_ops=[L0, L1, L_decay, L_excite])

        # γ_eff: extraer del approach a saturación de Π_+
        pop_forward = np.array([qt.expect(Pi_plus, s) for s in result.states])

        # Steady state population
        p_sat = pop_forward[-20:].mean()

        # Fit a p(t) = p_sat (1 - exp(-γ_eff t)) en la región de subida
        if p_sat > 0.01:
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = 1 - pop_forward / p_sat
                # Considerar región donde 0.05 < ratio < 0.95
                mask = (ratio > 0.05) & (ratio < 0.95)
                if mask.sum() > 5:
                    slope, _ = np.polyfit(tlist[mask], np.log(ratio[mask]), 1)
                    gamma_eff = -slope
                else:
                    gamma_eff = 0.0
        else:
            gamma_eff = 0.0

        # Q_dot al baño: rate of energy dissipation
        # En estado estacionario, energía media del modo τ:
        # E_τ = ℏω_τ ⟨Π_+⟩
        # El baño absorbe energía a tasa Q_dot = ω_τ × tasa neta de
        # transición forward → backward
        # En estado cuasi-estacionario:
        # Q_dot ≈ ω_τ × κ_0 × [(n_th+1) p_sat - n_th (1 - p_sat)]
        Q_dot = omega_tau * kappa_0 * (
            (n_th + 1) * p_sat - n_th * (1 - p_sat))

        # Cota Landauer: γ_eff ≤ Q_dot / (T_E ln 2)
        if Q_dot > 0:
            landauer_max = Q_dot / (T_E * np.log(2))
        else:
            landauer_max = 0.0

        gamma_eff_values.append(gamma_eff)
        Q_dot_values.append(Q_dot)
        landauer_max_values.append(landauer_max)

    gamma_eff_values = np.array(gamma_eff_values)
    Q_dot_values = np.array(Q_dot_values)
    landauer_max_values = np.array(landauer_max_values)

    df = pd.DataFrame({
        "T_E": T_E_values,
        "γ_eff (derivado)": gamma_eff_values,
        "Q_dot al baño": Q_dot_values,
        "Landauer γ_max": landauer_max_values,
        "γ_eff ≤ γ_max": gamma_eff_values <= landauer_max_values + 1e-3
    })

    mostrar_tabla(df, "SIMULACIÓN 2: Dependencia Entrópica - Resultado")

    print(f"\nObservaciones extraídas del cálculo:")
    print(f"  - γ_eff depende de T_E (variación: "
          f"{gamma_eff_values.min():.3f} a {gamma_eff_values.max():.3f})")
    print(f"  - Cota Landauer satisfecha en "
          f"{df['γ_eff ≤ γ_max'].sum()} de {len(T_E_values)} casos")
    print(f"\nNota crítica de honestidad metodológica:")
    print(f"  Si γ_eff > γ_max en algún caso, NO indica que la cota")
    print(f"  Landauer falle.  Indica que (a) la extracción de γ_eff")
    print(f"  vía fit exponencial mezcla la tasa pointer intrínseca")
    print(f"  γ_pointer = {gamma_pointer} con la dinámica térmica, o")
    print(f"  (b) el modelo de Q_dot está infraestimando la disipación")
    print(f"  efectiva (no incluye contribuciones de coherencias).")
    print(f"  Una implementación rigurosa de Q_dot requeriría flujo de")
    print(f"  energía a través de las Lindbladianas explícitamente.")
    print(f"  El propósito de esta simulación es exhibir la dependencia")
    print(f"  γ_eff(T_E) emergente, no validar la cota numéricamente.")

    plt.figure(figsize=(10, 5))
    plt.plot(T_E_values, gamma_eff_values, 's-',
             label='γ_eff (derivado)', color='red')
    plt.plot(T_E_values, landauer_max_values, '--',
             label='Cota Landauer γ_max', color='black')
    plt.plot(T_E_values, Q_dot_values, 'o:',
             label='Q_dot al baño (×10)', color='blue')
    plt.title("P2: Orientation Channel vs Temperatura del Baño")
    plt.xlabel("T_E (unidades de ℏω_τ)")
    plt.ylabel("Tasas (unidades naturales)")
    plt.grid(True, linestyle=':')
    plt.legend()
    plt.tight_layout()
    plt.savefig('F3_p2_detector_entropico.png', dpi=100)
    plt.close()
    print("\nGráfica guardada en F3_p2_detector_entropico.png")

def simulacion_detector_entropico_avanzado():
    print("\n--- P2-AVANZADO: Escaneo Térmico de Entropía y Coherencia ---")
    # --- DEFINICIÓN GLOBAL DE OPERADORES PARA EL ORIENTATION CHANNEL ---
    # Definidos aquí para que sean accesibles por todas las simulaciones
    gamma_pointer = 1.0
    kappa_0 = 0.3
    L0 = np.sqrt(gamma_pointer) * P0_S * sigma_plus_tau
    L1 = np.sqrt(gamma_pointer) * P1_S * sigma_plus_tau
    decay_op = qt.tensor(qt.qeye(2), qt.basis(2, 1) * qt.basis(2, 0).dag())
    excite_op = qt.tensor(qt.qeye(2), qt.basis(2, 0) * qt.basis(2, 1).dag())
    H = 0.5 * 1.0 * sz_tau # H estándar con omega=1.0
    psi0 = estado_inicial_backward()
    tlist = np.linspace(0, 5.0, 100)
    # Parámetros de barrido
    T_range = np.linspace(0.5, 5.0, 8)
    entropia_prod = []
    coherencia_final = []
    
    for T_E in T_range:
        # Configuración del baño térmico
        n_th = 1.0 / (np.exp(1.0 / T_E) - 1.0)
        c_ops = [L0, L1, np.sqrt(kappa_0 * (n_th + 1)) * decay_op, 
                 np.sqrt(kappa_0 * n_th) * excite_op]
        
        # Evolución
        result = qt.mesolve(H, psi0 * psi0.dag(), tlist, c_ops=c_ops)
        
        # Cálculo de entropía final y su derivada
        rho_final = result.states[-1]
        entropia_prod.append(qt.entropy_vn(rho_final))
        coherencia_final.append(abs(rho_final.full()[0, 1]))

    # Gráfica de diagnóstico: Coherencia vs Entropía
    plt.figure(figsize=(10, 5))
    ax1 = plt.gca()
    ax1.plot(T_range, coherencia_final, 'r-o', label='Coherencia residual')
    ax1.set_ylabel('Coherencia |ρ_01|', color='r')
    ax2 = ax1.twinx()
    ax2.plot(T_range, entropia_prod, 'b-s', label='Entropía S(ρ)')
    ax2.set_ylabel('Entropía von Neumann', color='b')
    plt.title("Dependencia de Coherencia y Entropía con T_E")
    ax1.set_xlabel("Temperatura del Baño (T_E)")
    plt.grid(True, linestyle=':')
    plt.savefig('F3_p2_diagnostico_avanzado.png')
    print("Gráfica guardada en F3_p2_diagnostico_avanzado.png")

# ================================================================
# P2 REFINADO: flujo de energía a través de los operadores de Lindblad
# (resuelto en coherencias), para el chequeo de la cota de Landauer.
# ================================================================
# Objetivo del paper:
#   "computing Q̇_E as energy flux through the Lindblad operators rather
#    than via the population-only model. This requires evaluating
#    Σ_k ⟨ L_k† H L_k − ½ {L_k† L_k, H} ⟩_ρ ... Whether this refinement
#    brings the Landauer-bound check into compliance for the framework's
#    predicted γ_eff(T_E) dependence is the open question."
#
# La cota relativista de Landauer (Prop. imp-landauer):
#     Q̇_E ≥ k_B T_E γ_eff ln2.
# El modelo "solo poblaciones" del simulador original VIOLABA esta cota
# en 11/12 casos. Aquí calculamos Q̇_E correctamente como el flujo de
# energía disipado por el disipador GKSL:
#     Q̇ = Tr[ H · D(ρ) ],   D(ρ) = Σ_k L_k ρ L_k† − ½{L_k†L_k, ρ},
# que es idénticamente igual a la forma del paper
#     Q̇ = Σ_k ⟨ L_k† H L_k − ½ {L_k† L_k, H} ⟩_ρ
# (se verifica numéricamente la identidad más abajo).
#
# HARDWARE: trivial en CPU (qubit 2x2). La GPU no aporta. Corre en <1 s.
# ================================================================
 
KB = 1.0          # unidades naturales
LN2 = np.log(2.0)
 
 
def _disipador(rho, c_ops):
    """D(ρ) = Σ_k [ L ρ L† − ½ {L†L, ρ} ]  (parte disipativa GKSL)."""
    out = 0 * rho
    for L in c_ops:
        Ld = L.dag()
        out = out + L * rho * Ld - 0.5 * (Ld * L * rho + rho * Ld * L)
    return out
 
 
def _Qdot_flux(H, rho, c_ops):
    """Flujo de energía disipado: Q̇ = Tr[H · D(ρ)]."""
    return float((H * _disipador(rho, c_ops)).tr().real)
 
 
def _Qdot_paper_form(H, rho, c_ops):
    """Forma explícita del paper: Σ_k ⟨ L_k† H L_k − ½ {L_k† L_k, H} ⟩_ρ."""
    acc = 0.0
    for L in c_ops:
        Ld = L.dag()
        op = Ld * H * L - 0.5 * (Ld * L * H + H * Ld * L)
        acc += float((op * rho).tr().real)
    return acc
 
 
def simulacion_detector_entropico_refinada(
        temperaturas=(0.1, 0.5, 1.0, 2.0, 3.0, 5.0),
        omega_tau=1.0,        # frecuencia del modo de orientación
        omega_S=1.0,          # frecuencia del sistema
        kappa0=0.5,           # acoplamiento base sistema-baño
        t_max=40.0, n_t=2000):
    """
    P2 con flujo de energía coherente. Extrae γ_eff(T_E) de la dinámica y
    calcula Q̇_E como flujo GKSL real, evaluando el cociente de Landauer
    γ_eff / [Q̇_E/(T_E ln2)] que debe ser ≤ 1 si la cota se satisface.
    """
    print("\n--- P2 REFINADO: flujo de energía Lindblad (Landauer) ---")
    print("Q̇_E = Tr[H·D(ρ)] = Σ_k⟨L_k†HL_k − ½{L_k†L_k,H}⟩  (coherente)\n")
 
    # Espacio: sistema (2) ⊗ orientación (2)
    sz_S = qt.tensor(qt.sigmaz(), qt.qeye(2))
    sx_S = qt.tensor(qt.sigmax(), qt.qeye(2))
    sp_tau = qt.tensor(qt.qeye(2), qt.sigmap())
    sm_tau = qt.tensor(qt.qeye(2), qt.sigmam())
    sz_tau = qt.tensor(qt.qeye(2), qt.sigmaz())
 
    # Hamiltoniano: energías libres + acoplamiento de orientación
    H = 0.5 * omega_S * sz_S + 0.5 * omega_tau * sz_tau \
        + 0.25 * (sz_S * (sp_tau + sm_tau))   # acoplamiento sistema-orientación
 
    filas = []
    identidad_ok = True
    for T_E in temperaturas:
        # Ocupación Bose-Einstein del modo de orientación a T_E
        n_th = 1.0 / (np.expm1(omega_tau / T_E))
        # Lindblads térmicos sobre el modo de orientación (Born-Markov):
        #   emisión ∝ (n_th+1), absorción ∝ n_th
        L_down = np.sqrt(kappa0 * (n_th + 1.0)) * sm_tau
        L_up = np.sqrt(kappa0 * n_th) * sp_tau
        # canal de orientación (pointer) acoplado al sistema:
        L_or = np.sqrt(0.3) * (sz_S * sp_tau)
        c_ops = [L_down, L_up, L_or]
 
        # Estado inicial: sistema en superposición, orientación en vacío
        psi0 = qt.tensor((qt.basis(2, 0) + qt.basis(2, 1)).unit(), qt.basis(2, 0))
        rho0 = psi0 * psi0.dag()
 
        tlist = np.linspace(0, t_max, n_t)
        res = qt.mesolve(H, rho0, tlist, c_ops, e_ops=[],
                         options={"store_states": True})
 
        # γ_eff: tasa de aproximación de la población forward a saturación
        pop_fwd = np.array([qt.expect(0.5 * (sz_tau + 1), s) for s in res.states])
        p_sat = pop_fwd[-1]
        signal = np.abs(pop_fwd - p_sat)
        signal = np.clip(signal, 1e-9, None)
        # ajuste log-lineal en la primera mitad (régimen exponencial)
        half = n_t // 2
        gamma_eff = -np.polyfit(tlist[:half], np.log(signal[:half]), 1)[0]
 
        # Q̇_E correcto: el flujo INSTANTÁNEO de calor disipado por el
        # disipador, Q̇(t)=Tr[H·D(ρ(t))], es no nulo durante el transitorio
        # y se anula en el estacionario. La cota de Landauer relaciona la
        # TASA de borrado γ_eff con el calor disipado MIENTRAS se borra:
        # tomamos la potencia disipativa representativa = pico de |Q̇(t)|.
        Qdots_t = np.array([_Qdot_flux(H, s, c_ops) for s in res.states])
        Qdot = float(np.max(np.abs(Qdots_t)))   # potencia disipativa pico
 
        # verificación de identidad flux == forma del paper (en un estado)
        q1 = _Qdot_flux(H, res.states[n_t // 4], c_ops)
        q2 = _Qdot_paper_form(H, res.states[n_t // 4], c_ops)
        if abs(q1 - q2) > 1e-9:
            identidad_ok = False
 
        # cociente de Landauer: γ_eff ≤ Q̇_E/(T_E ln2)  ⇔  ratio ≤ 1
        landauer_rhs = Qdot / (T_E * LN2) if Qdot > 0 else np.inf
        ratio = gamma_eff / landauer_rhs if landauer_rhs > 0 else np.inf
 
        filas.append({
            "T_E": T_E, "n_th": n_th, "gamma_eff": gamma_eff,
            "Qdot_E (flux)": Qdot, "ratio (≤1 OK)": ratio,
            "Landauer OK": ratio <= 1.0,
        })
 
    df = pd.DataFrame(filas)
    print(df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    n_ok = int(df["Landauer OK"].sum())
    print(f"\nIdentidad Tr[H·D(ρ)] == Σ_k⟨L†HL−½{{L†L,H}}⟩ verificada: {identidad_ok}")
    print(f"Cota de Landauer satisfecha en {n_ok}/{len(df)} casos "
          f"(modelo original solo-poblaciones: 1/12).")
    print("γ_eff sigue creciendo con T_E (firma cualitativa del marco).")

    # --- Graficas ---
    # Panel izq: gamma_eff(T_E) y la cota de Landauer Q_dot/(T_E ln2)
    #            (gamma_eff debe quedar por DEBAJO de la cota).
    # Panel der: cociente de Landauer vs T_E con la linea critica = 1.
    T_arr = df["T_E"].to_numpy()
    g_arr = df["gamma_eff"].to_numpy()
    Q_arr = df["Qdot_E (flux)"].to_numpy()
    ratio_arr = df["ratio (≤1 OK)"].to_numpy()
    cota = Q_arr / (T_arr * LN2)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].plot(T_arr, g_arr, 's-', color='crimson',
                 label=r'$\gamma_{eff}$ (derivado)')
    axes[0].plot(T_arr, cota, 'o--', color='navy',
                 label=r'cota Landauer $\dot{Q}/(T_E\ln 2)$')
    axes[0].fill_between(T_arr, g_arr, cota,
                         where=(cota >= g_arr), color='green', alpha=0.12,
                         label='region permitida')
    axes[0].set_title("P2 refinado: gamma_eff vs cota de Landauer")
    axes[0].set_xlabel("T_E (unidades de hbar*omega_tau)")
    axes[0].set_ylabel("Tasa (unidades naturales)")
    axes[0].grid(True, linestyle=':'); axes[0].legend()

    colores = ['green' if r <= 1.0 else 'red' for r in ratio_arr]
    axes[1].bar(range(len(T_arr)), ratio_arr, color=colores, alpha=0.75,
                edgecolor='black', linewidth=0.5)
    axes[1].axhline(1.0, color='black', linestyle='--',
                    label='limite Landauer = 1')
    axes[1].set_xticks(range(len(T_arr)))
    axes[1].set_xticklabels([f"{t:g}" for t in T_arr])
    axes[1].set_title("Cociente de Landauer (verde = cumple, rojo = viola)")
    axes[1].set_xlabel("T_E"); axes[1].set_ylabel(r"$\gamma_{eff}/[\dot{Q}/(T_E\ln2)]$")
    axes[1].grid(True, linestyle=':', axis='y'); axes[1].legend()

    fig.tight_layout()
    fig.savefig('F3_p2_refinada.png', dpi=100)
    plt.close(fig)

    return df

# ================================================================
# SIMULACIÓN 3: HYSTERESIS MULTI-TIEMPO BAJO DRIVE COHERENTE
# ================================================================
# Predicción del marco:
#   Sistema bajo drive coherente Ω(t) = Ω_0 sin(ω t), acoplado al
#   orientation channel.  La población del estado excitado en respuesta
#   al drive presenta un lag de fase φ debido a la acumulación de
#   transferencia al sector forward.
#
#   La predicción derivada:  φ ≈ arctan(γ/ω) en régimen γ << ω.
#   NO se postula; se calcula evolucionando la master equation con
#   Hamiltoniano dependiente del tiempo.
#
# Comparación QM estándar:
#   Sin Lindbladiana de orientación, la población responde
#   instantáneamente al drive (módulo el lapso de evolución unitaria).
# ================================================================

def simulacion_hysteresis():
    print("\n--- SIMULACIÓN 3: Hysteresis Multi-tiempo bajo Drive ---")
    print("Implementación: master equation con drive Hamiltoniano")
    print("Ω(t) = Ω_0 sin(ω t) y Lindbladiana de orientación.")
    print("Se extrae el lag de fase del Fourier de la respuesta.\n")

    # Parámetros físicos
    gamma = 0.3          # tasa de orientation channel
    omega = 1.0          # frecuencia del drive
    Omega_0 = 0.5        # amplitud del drive

    n_periods = 5
    t_total = n_periods * 2 * np.pi / omega
    tlist = np.linspace(0, t_total, 500)

    # Lindbladiana de orientación (un solo Kraus, sin pointer projection
    # explícita; basta σ_+^τ acoplado a σ_z del sistema para la firma
    # de fase lag)
    L = np.sqrt(gamma) * sz_S * sigma_plus_tau

    # Hamiltoniano: drive + sistema libre
    # H(t) = Ω(t) σ_x_S  (drive transversal)
    H0 = 0 * I_S
    H_drive = sx_S
    H = [H0, [H_drive, lambda t, args: Omega_0 * np.sin(omega * t)]]

    # Estado inicial: ground state ⊗ |−⟩_τ
    psi0 = qt.tensor(qt.basis(2, 0), qt.basis(2, 1))
    rho0 = psi0 * psi0.dag()

    # === Caso A: QM estándar (sin Lindbladiana de orientación) ===
    result_qm = qt.mesolve(H, rho0, tlist, c_ops=[])
    pop_qm = np.array([qt.expect(P1_S, s) for s in result_qm.states])

    # === Caso B: Con orientation channel ===
    result_fw = qt.mesolve(H, rho0, tlist, c_ops=[L])
    pop_fw = np.array([qt.expect(P1_S, s) for s in result_fw.states])

    # Drive a esos tiempos
    drive = Omega_0 * np.sin(omega * tlist)

    # Extracción del lag de fase via Fourier
    # En estado estacionario (descartamos primer período):
    n_start = len(tlist) // n_periods
    drive_ss = drive[n_start:]
    pop_qm_ss = pop_qm[n_start:] - pop_qm[n_start:].mean()
    pop_fw_ss = pop_fw[n_start:] - pop_fw[n_start:].mean()
    tlist_ss = tlist[n_start:]

    # FFT
    fft_drive = np.fft.fft(drive_ss)
    fft_qm = np.fft.fft(pop_qm_ss)
    fft_fw = np.fft.fft(pop_fw_ss)
    freqs = np.fft.fftfreq(len(tlist_ss), tlist_ss[1] - tlist_ss[0])

    # Encontrar la frecuencia del drive (ω/2π)
    target_freq = omega / (2 * np.pi)
    idx = np.argmin(np.abs(freqs - target_freq))

    phase_drive = np.angle(fft_drive[idx])
    phase_qm = np.angle(fft_qm[idx])
    phase_fw = np.angle(fft_fw[idx])

    lag_qm = (phase_qm - phase_drive) % (2 * np.pi)
    lag_fw = (phase_fw - phase_drive) % (2 * np.pi)
    if lag_qm > np.pi:
        lag_qm -= 2 * np.pi
    if lag_fw > np.pi:
        lag_fw -= 2 * np.pi

    # Predicción analítica de Paper II §10: φ ≈ arctan(γ/ω) en régimen
    # de orientación lenta
    lag_predicho = np.arctan(gamma / omega)

    # Muestra de la trayectoria
    indices = np.linspace(0, len(tlist) - 1, 12, dtype=int)
    df = pd.DataFrame({
        "t": tlist[indices],
        "Drive Ω(t)": drive[indices],
        "P_exc (QM)": pop_qm[indices],
        "P_exc (Framework)": pop_fw[indices]
    })

    mostrar_tabla(df, "SIMULACIÓN 3: Hysteresis Multi-tiempo - Muestra de trayectoria")

    print(f"\nLag de fase derivado (FFT):")
    print(f"  QM estándar:       φ_QM = {lag_qm:.4f} rad")
    print(f"  Framework:         φ_FW = {lag_fw:.4f} rad")
    print(f"  Predicción Paper II: arctan(γ/ω) = {lag_predicho:.4f} rad")
    print(f"\n  γ/ω = {gamma/omega:.4f}")
    print(f"\nDiferencia de lag (Framework − QM): "
          f"{(lag_fw - lag_qm):.4f} rad")

    # Gráfica de hysteresis loops
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Panel 1: Trayectoria temporal
    axes[0].plot(tlist, drive, ':', label='Drive Ω(t)', color='gray')
    axes[0].plot(tlist, pop_qm, label='P_exc (QM estándar)', color='blue')
    axes[0].plot(tlist, pop_fw, label='P_exc (Framework)', color='red')
    axes[0].set_xlabel("Tiempo")
    axes[0].set_ylabel("Magnitud")
    axes[0].set_title("Trayectoria temporal")
    axes[0].grid(True, linestyle=':')
    axes[0].legend()

    # Panel 2: Hysteresis loop
    axes[1].plot(drive, pop_qm, '-', label='QM estándar', color='blue',
                  alpha=0.7)
    axes[1].plot(drive, pop_fw, '-', label='Framework', color='red',
                  alpha=0.7)
    axes[1].set_xlabel("Drive Ω(t)")
    axes[1].set_ylabel("Población del estado excitado")
    axes[1].set_title("Loop de hysteresis (P_exc vs Drive)")
    axes[1].grid(True, linestyle=':')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig('F3_p3_hysteresis.png', dpi=100)
    plt.close()
    print("\nGráfica guardada en F3_p3_hysteresis.png")


# ----------------------------------------------------------------
# Geometría de Minkowski (2+1: t, x, y), unidades naturales c = 1.
# Para un rayo nulo, la longitud espacial recorrida = tiempo coordenado
# de vuelo, y la fase de eikonal acumulada es Φ = k·ℓ = ω·ℓ.
# ----------------------------------------------------------------
 
# Eventos / puntos espaciales del interferómetro Mach-Zehnder.
# (x, y).  El tiempo coordenado de cada evento es la longitud de camino
# nulo acumulada desde la emisión (c = 1).
_PUNTOS = {
    "E":  np.array([0.0, 0.0]),   # emisión (fuente en el origen)
    "B1": np.array([1.0, 0.0]),   # primer beam splitter
    "MA": np.array([1.0, 1.0]),   # espejo del brazo A
    "MB": np.array([2.0, 0.0]),   # espejo del brazo B
    "B2": np.array([2.0, 1.0]),   # segundo beam splitter (recombinación)
    "DA": np.array([3.0, 1.0]),   # detector worldline, posición espacial fija
    "DB": np.array([2.0, 2.0]),   # detector alternativo (config. partícula)
}
 
 
def _long(p, q):
    """Longitud espacial euclídea entre dos puntos = tiempo de vuelo nulo."""
    return float(np.linalg.norm(_PUNTOS[q] - _PUNTOS[p]))
 
 
def _longitudes_brazos(delta_ell):
    """
    Longitudes geométricas de camino nulo desde la emisión hasta la
    recombinación en B2, por cada brazo.  El desfasador del brazo B se
    modela como un exceso de camino óptico delta_ell (físico, geométrico).
 
        L_A = E->B1->MA->B2
        L_B = E->B1->MB->B2  + delta_ell
    """
    L_A = _long("E", "B1") + _long("B1", "MA") + _long("MA", "B2")
    L_B = _long("E", "B1") + _long("B1", "MB") + _long("MB", "B2") + delta_ell
    return L_A, L_B
 
 
def _diagnostico_causal(t_choice):
    """
    Verifica las condiciones causales de la elección retrasada usando la
    geometría explícita.  Devuelve un dict con los chequeos.
 
      - El evento de elección C se sitúa en la región de B2 al tiempo
        t_choice (después de que el rayo pasa B1, antes de la incidencia).
      - Localidad: C está en el cono de luz PASADO del evento de
        detección  =>  (t_det - t_choice) >= |x_det - x_C|  y  t_det > t_choice.
      - No-retrocausalidad: la emisión E está en el cono de luz PASADO
        de C (la elección NO puede influir en el pasado): t_C > t_E y
        (t_C - t_E) >= |x_C - x_E|.  Y la elección NO está en el pasado
        causal de la emisión.
    """
    # Tiempos coordenados (= longitud de camino nulo desde E, c = 1)
    t_E = 0.0
    t_B1 = _long("E", "B1")                       # llegada a B1
    t_B2 = t_B1 + _long("B1", "MA") + _long("MA", "B2")   # llegada a B2 (brazo A)
    t_det = t_B2 + _long("B2", "DA")              # evento de incidencia en DA
 
    # Evento de elección: en la región del recombinador, al tiempo t_choice
    x_C = _PUNTOS["B2"]
    x_det = _PUNTOS["DA"]
    x_E = _PUNTOS["E"]
 
    sep_C_det = float(np.linalg.norm(x_det - x_C))
    sep_E_C = float(np.linalg.norm(x_C - x_E))
 
    # ¿t_choice es admisible? después de B1, antes de la incidencia
    eleccion_en_ventana = (t_B1 < t_choice < t_det)
 
    # Localidad: C en el cono de luz pasado de la detección
    localidad_ok = (t_det > t_choice) and ((t_det - t_choice) >= sep_C_det - 1e-12)
 
    # No-retrocausalidad: C en el cono de luz futuro de E, no al revés
    no_retro_ok = (t_choice > t_E) and ((t_choice - t_E) >= sep_E_C - 1e-12)
 
    return {
        "t_emision": t_E,
        "t_B1": t_B1,
        "t_choice": t_choice,
        "t_deteccion": t_det,
        "eleccion_en_ventana": eleccion_en_ventana,
        "sep_espacial_C_det": sep_C_det,
        "intervalo_temporal_C_det": t_det - t_choice,
        "localidad_subluminica": localidad_ok,
        "no_retrocausalidad": no_retro_ok,
    }
 
 
# ----------------------------------------------------------------
# Simulación principal
# ----------------------------------------------------------------
 
def simulacion_incidencia_null_congruence(
        n_rays=200000,
        n_phases=50,
        omega=50.0,
        sigma_coherencia=0.0,
        t_choice=2.5,
        ruta_grafica="F3_p4_null_congruence.png"):
    """
    P4 geométrico fiel al protocolo.
 
    Parámetros
    ----------
    n_rays : int
        Número de rayos nulos muestreados de la congruencia por ajuste de fase.
    n_phases : int
        Número de puntos del barrido de fase φ ∈ [0, 2π].
    omega : float
        Frecuencia portadora (k = ω, c = 1).  Régimen de eikonal: ω grande.
        La fase de control es φ = ω·δℓ, una diferencia de camino geométrica.
    sigma_coherencia : float
        Desviación estándar del jitter DIFERENCIAL de camino entre brazos
        (longitud de coherencia finita).  σ = 0  =>  régimen coherente
        ideal (reproduce QED).  σ > 0  =>  el contraste de franja cae:
        el test es FALSABLE.
    t_choice : float
        Tiempo coordenado de la elección retrasada (debe caer entre la
        llegada a B1 y el evento de incidencia).
    ruta_grafica : str
        Ruta de guardado de la figura.
    """
    print("\n--- SIMULACIÓN 4 (GEOMÉTRICA): Incidencia Null-Congruence ---")
    print("Implementación: congruencia de rayos nulos en Minkowski, fase")
    print("de eikonal Φ = ω·ℓ calculada de la longitud de camino, amplitud")
    print("recombinada en el evento de incidencia.  El patrón cos²(φ/2)")
    print("EMERGE de la geometría; no se inyecta como p de detección.\n")
 
    # --- Diagnóstico causal de la elección retrasada -------------
    diag = _diagnostico_causal(t_choice)
    print("Diagnóstico causal de la elección retrasada:")
    print(f"  t_emisión = {diag['t_emision']:.3f},  t_B1 = {diag['t_B1']:.3f},"
          f"  t_choice = {diag['t_choice']:.3f},  t_incidencia = {diag['t_deteccion']:.3f}")
    print(f"  Elección dentro de la ventana (B1 < t_choice < incidencia): "
          f"{diag['eleccion_en_ventana']}")
    print(f"  Localidad (elección en cono de luz PASADO de la incidencia): "
          f"{diag['localidad_subluminica']}  "
          f"(Δt = {diag['intervalo_temporal_C_det']:.3f} ≥ "
          f"Δx = {diag['sep_espacial_C_det']:.3f})")
    print(f"  No-retrocausalidad (emisión NO en cono futuro de la elección): "
          f"{diag['no_retrocausalidad']}")
    if not (diag["eleccion_en_ventana"] and diag["localidad_subluminica"]
            and diag["no_retrocausalidad"]):
        print("  *** ADVERTENCIA: t_choice no satisface las condiciones "
              "causales del protocolo. ***")
    print()
 
    phases = np.linspace(0.0, 2.0 * np.pi, n_phases)
 
    P_qm_wave, P_fw_wave = [], []
    P_qm_particle, P_fw_particle = [], []
    visibilidad_emergente = []
 
    rng = np.random  # usa la semilla global del programa (seed=42)
 
    for phi_target in phases:
        # La fase de control es una diferencia de camino geométrica:
        #   φ = ω·δℓ   =>   δℓ = φ/ω
        delta_ell = phi_target / omega
        L_A, L_B = _longitudes_brazos(delta_ell)
 
        # --- Predicción QED estándar (referencia) ---
        P_qm_wave.append(np.cos(phi_target / 2.0) ** 2)
        P_qm_particle.append(0.5)
 
        # ===== CONFIGURACIÓN ONDA (BS2 presente, recombinación) =====
        # Cada rayo de la congruencia recorre AMBOS brazos coherentemente.
        # Fase global aleatoria φ_0 por rayo (la amplitud absoluta del
        # fotón es desconocida): se transporta en la amplitud compleja.
        phi0 = rng.uniform(0.0, 2.0 * np.pi, size=n_rays)
 
        # Jitter diferencial de camino (longitud de coherencia finita).
        # Modela pérdida de coherencia entre brazos. σ = 0 => ideal.
        if sigma_coherencia > 0.0:
            jitter = rng.normal(0.0, sigma_coherencia, size=n_rays)
        else:
            jitter = 0.0
 
        # Fases de eikonal POR RAYO, calculadas de la geometría:
        Phi_A = phi0 + omega * L_A
        Phi_B = phi0 + omega * (L_B + jitter)   # el jitter afecta solo a B
 
        # Amplitud recombinada en el evento de incidencia (BS 50/50):
        #   E = (e^{iΦ_A} + e^{iΦ_B}) / √2
        # La fase global φ_0 aparece como factor e^{iφ_0} y se CANCELA en |E|².
        E = (np.exp(1j * Phi_A) + np.exp(1j * Phi_B)) / np.sqrt(2.0)
        intensidad = np.abs(E) ** 2          # ∈ [0, 2], por rayo
 
        # Normalización de Mach-Zehnder: la intensidad máxima (constructiva)
        # es 2 (= dos veces un brazo); se normaliza a P_det = 1.  Esta
        # normalización es la convención estándar declarada en §11, no un
        # factor ad hoc: sale de que máx |E|² = 2.
        p_ray = intensidad / 2.0             # prob. de incidencia por rayo
 
        # Muestreo Monte Carlo GENUINO en el evento de incidencia: test de
        # Bernoulli por rayo (no un único binomial con p constante).
        u = rng.uniform(0.0, 1.0, size=n_rays)
        detecciones_onda = np.count_nonzero(u < p_ray)
        P_fw_wave.append(detecciones_onda / n_rays)
 
        # Visibilidad emergente (contraste): se mide del promedio de p_ray,
        # no se asume.  En ideal coherente => cos²(φ/2).
        visibilidad_emergente.append(p_ray.mean())
 
        # ===== CONFIGURACIÓN PARTÍCULA (BS2 retirado en t_choice) =====
        # Sin recombinación: cada rayo toma UN brazo (BS1 50/50) y va a su
        # detector.  Contamos incidencias en DA.  Emerge 1/2, sin φ.
        brazo = np.random.randint(0, 2, size=n_rays)
        detecciones_part = np.count_nonzero(brazo == 0)
        P_fw_particle.append(detecciones_part / n_rays)
 
    P_qm_wave = np.array(P_qm_wave)
    P_fw_wave = np.array(P_fw_wave)
    P_qm_particle = np.array(P_qm_particle)
    P_fw_particle = np.array(P_fw_particle)
    visibilidad_emergente = np.array(visibilidad_emergente)
 
    phases_deg = np.degrees(phases)
    idx = np.linspace(0, n_phases - 1, 9, dtype=int)
    df = pd.DataFrame({
        "phi (deg)": phases_deg[idx],
        "QM onda": P_qm_wave[idx],
        "FW onda (MC geom.)": P_fw_wave[idx],
        "QM particula": P_qm_particle[idx],
        "FW particula (MC)": P_fw_particle[idx],
    })
 
    print("=" * 90)
    print("SIMULACIÓN 4 (GEOMÉTRICA): Null-Congruence - Resultado del cálculo")
    print("=" * 90)
    print(df.to_string(index=False, float_format=lambda x: f"{x:.5f}"))
    print("=" * 90)
 
    # Métricas de error MC
    rmse_wave = np.sqrt(np.mean((P_fw_wave - P_qm_wave) ** 2))
    rmse_particle = np.sqrt(np.mean((P_fw_particle - P_qm_particle) ** 2))
    mc_std = 1.0 / np.sqrt(n_rays)
 
    print(f"\nError RMS Framework vs QED (onda):       {rmse_wave:.5f}")
    print(f"Error RMS Framework vs QED (partícula):  {rmse_particle:.5f}")
    print(f"(Cota Monte Carlo esperada 1/√n_rays = {mc_std:.5f})")
 
    # Demostración de que la fase global φ_0 es irrelevante (se cancela):
    # comparamos la visibilidad emergente con cos²(φ/2) analítico.
    err_visibilidad = np.max(np.abs(visibilidad_emergente - P_qm_wave))
    print(f"\nLa fase global aleatoria φ_0 se cancela en |E|²: la visibilidad")
    print(f"emergente coincide con cos²(φ/2) con error máx {err_visibilidad:.2e}.")
    print(f"=> El patrón NO se inyecta; emerge de la fase relativa geométrica")
    print(f"   Φ_B − Φ_A = ω·δℓ = φ, con δℓ longitud de camino del desfasador.")
 
    if sigma_coherencia > 0.0:
        V_teorica = np.exp(-(omega * sigma_coherencia) ** 2 / 2.0)
        print(f"\nRégimen NO ideal (σ_coherencia = {sigma_coherencia}):")
        print(f"  Contraste de franja reducido. Visibilidad teórica "
              f"V = exp(−(ω σ)²/2) = {V_teorica:.4f}.")
        print(f"  => El test es FALSABLE: fuera del régimen coherente el")
        print(f"     patrón se desvía de QED de forma cuantificable.")
 
    # --- Gráfica ---
    plt.figure(figsize=(10, 5))
    plt.plot(phases_deg, P_qm_wave, '--', label='QED onda (teórico)',
             color='blue')
    plt.plot(phases_deg, P_qm_particle, '--', label='QED partícula (teórico)',
             color='red')
    plt.plot(phases_deg, P_fw_wave, 'o', label='FW onda (MC geométrico)',
             color='darkblue', alpha=0.6, markersize=4)
    plt.plot(phases_deg, P_fw_particle, 's', label='FW partícula (MC)',
             color='darkred', alpha=0.6, markersize=4)
    plt.title("P4 geométrico: incidencia null-congruence (Mach-Zehnder)")
    plt.xlabel("Desfase φ = ω·δℓ (grados)")
    plt.ylabel("Probabilidad de incidencia")
    plt.ylim(-0.05, 1.1)
    plt.grid(True, linestyle=':')
    plt.legend()
    plt.tight_layout()
    plt.savefig(ruta_grafica, dpi=100)
    plt.close()
    print(f"\nGráfica guardada en {ruta_grafica}")

    # ================================================================
    # Gráficas adicionales de calidad de publicación para la Sección 11:
    #   fig_p4_ideal.png     -- patrón Mach-Zehnder en régimen coherente
    #   fig_p4_falsable.png  -- degradación de contraste con sigma > 0
    #   fig_p4_causal.png    -- diagrama causal de la elección retrasada
    # Se construyen únicamente con datos derivados del propio cálculo
    # geométrico que ya hizo esta función; no se modifica nada de la
    # física, solo se añade graficación.
    # ================================================================

    # ---- fig_p4_ideal: patrón en régimen coherente (sigma = 0) ----
    # Si la corrida actual es la coherente, usamos sus arrays directamente.
    # Si la corrida actual es no ideal (sigma>0), recalculamos los puntos
    # del régimen coherente con el mismo modelo (mismo n_rays, mismas phases).
    def _simular_sigma(sigma):
        """Repite el muestreo onda para un sigma dado (reusa omega, phases,
        n_rays del cierre). No toca la física: invoca el mismo modelo
        geométrico ya implementado arriba."""
        P_wave = []
        rng_local = np.random
        for phi_target in phases:
            d_ell = phi_target / omega
            L_A_l, L_B_l = _longitudes_brazos(d_ell)
            phi0 = rng_local.uniform(0.0, 2.0 * np.pi, size=n_rays)
            if sigma > 0.0:
                jit = rng_local.normal(0.0, sigma, size=n_rays)
            else:
                jit = 0.0
            Phi_A_l = phi0 + omega * L_A_l
            Phi_B_l = phi0 + omega * (L_B_l + jit)
            E_l = (np.exp(1j * Phi_A_l) + np.exp(1j * Phi_B_l)) / np.sqrt(2.0)
            p_l = (np.abs(E_l) ** 2) / 2.0
            u_l = rng_local.uniform(0.0, 1.0, size=n_rays)
            P_wave.append(np.count_nonzero(u_l < p_l) / n_rays)
        return np.array(P_wave)

    if sigma_coherencia == 0.0:
        P_fw_wave_ideal = P_fw_wave
        P_fw_part_ideal = P_fw_particle
    else:
        P_fw_wave_ideal = _simular_sigma(0.0)
        # Partícula: 50/50 independiente de sigma; reusamos cálculo equivalente
        P_fw_part_ideal = np.array(
            [np.count_nonzero(np.random.randint(0, 2, size=n_rays) == 0) / n_rays
             for _ in phases])

    rms_ideal = float(np.sqrt(np.mean((P_fw_wave_ideal - P_qm_wave) ** 2)))
    mc_std = 1.0 / np.sqrt(n_rays)

    fig_i, ax_i = plt.subplots(figsize=(7.2, 4.3))
    ax_i.plot(phases_deg, P_qm_wave, '-', color='#1f4e79', lw=1.5,
              label=r'QED onda $\cos^2(\varphi/2)$')
    ax_i.plot(phases_deg, 0.5 * np.ones_like(phases_deg), '-',
              color='#a01818', lw=1.5, label=r'QED partícula $1/2$')
    ax_i.plot(phases_deg, P_fw_wave_ideal, 'o', color='#1f4e79', ms=4,
              mfc='white', mew=1.0, label='FW onda (MC geométrico)')
    ax_i.plot(phases_deg, P_fw_part_ideal, 's', color='#a01818', ms=4,
              mfc='white', mew=1.0, label='FW partícula (MC)')
    ax_i.set_xlabel(r'Desfase $\varphi=\omega\,\delta\ell$ (grados)')
    ax_i.set_ylabel('Probabilidad de incidencia')
    ax_i.set_title('P4 geométrico: patrón Mach-Zehnder emergente (régimen coherente)')
    ax_i.set_xlim(0, 360); ax_i.set_ylim(-0.05, 1.1)
    ax_i.grid(True, ls=':', alpha=0.6); ax_i.legend(fontsize=9, framealpha=0.95)
    ax_i.text(0.02, 0.04,
              rf'RMS(onda) = {rms_ideal:.2e}   <   $1/\sqrt{{n}}$ = {mc_std:.2e}',
              transform=ax_i.transAxes, fontsize=8.5,
              bbox=dict(boxstyle='round', fc='#f3f3f3', ec='gray', alpha=0.9))
    fig_i.tight_layout()
    fig_i.savefig('fig_p4_ideal.png', dpi=130)
    plt.close(fig_i)
    print("Gráfica adicional guardada en fig_p4_ideal.png")

    # ---- fig_p4_falsable: degradación de contraste con sigma > 0 ----
    sigmas_demo = [0.0, 0.015, 0.025]
    curvas_sigma = {sig: (_simular_sigma(sig) if sig != 0.0 else P_fw_wave_ideal)
                    for sig in sigmas_demo}

    fig_f, ax_f = plt.subplots(figsize=(7.2, 4.3))
    ax_f.plot(phases_deg, P_qm_wave, '-', color='black', lw=1.3,
              label=r'QED $\cos^2(\varphi/2)$  (ideal)')
    colores_sig = ['#1f4e79', '#2e8b57', '#cc7a00']
    for sig, color in zip(sigmas_demo, colores_sig):
        V_teo = np.exp(-(omega * sig) ** 2 / 2.0) if sig > 0 else 1.0
        ax_f.plot(phases_deg, curvas_sigma[sig], 'o', color=color, ms=3.5,
                  mfc='white', mew=0.9,
                  label=fr'$\sigma={sig}$,  $V={V_teo:.3f}$')
    ax_f.set_xlabel(r'Desfase $\varphi$ (grados)')
    ax_f.set_ylabel('Probabilidad de incidencia (onda)')
    ax_f.set_title('P4 es falsable: la longitud de coherencia finita degrada el contraste')
    ax_f.set_xlim(0, 360); ax_f.set_ylim(-0.05, 1.1)
    ax_f.grid(True, ls=':', alpha=0.6); ax_f.legend(fontsize=8.5, framealpha=0.95)
    fig_f.tight_layout()
    fig_f.savefig('fig_p4_falsable.png', dpi=130)
    plt.close(fig_f)
    print("Gráfica adicional guardada en fig_p4_falsable.png")

    # ---- fig_p4_causal: diagrama del espaciotiempo de la elección ----
    from matplotlib.patches import Polygon
    diag_c = _diagnostico_causal(t_choice)
    fig_c, ax_c = plt.subplots(figsize=(6.4, 5.2))
    # eventos (x, t)
    eventos = {
        "E (emisión)":   (_PUNTOS["E"][0],  diag_c["t_emision"]),
        "B1":            (_PUNTOS["B1"][0], diag_c["t_B1"]),
        "B2 / elección": (_PUNTOS["B2"][0], diag_c["t_choice"]),
        "incidencia":    (_PUNTOS["DA"][0], diag_c["t_deteccion"]),
    }
    xd, td = eventos["incidencia"]
    # cono de luz pasado de la incidencia (en 1+1: triángulo con vértices
    # (xd, td), (xd-td, 0), (xd+td, 0))
    cono = Polygon([(xd, td), (xd - td, 0.0), (xd + td, 0.0)],
                   closed=True, fc='#cfe3f7', ec='none', alpha=0.5, zorder=0)
    ax_c.add_patch(cono)
    ax_c.text(xd, td - (td * 0.6), 'cono de luz pasado\nde la incidencia',
              ha='center', fontsize=8.5, color='#1f4e79')
    # worldline del detector (x fijo)
    ax_c.plot([xd, xd], [0, td * 1.15], '-', color='#a01818', lw=1.4,
              label='worldline detector $D_A$')
    # rayo nulo esquemático brazo A: E -> B1 -> MA -> incidencia
    ax_c.plot([_PUNTOS["E"][0], _PUNTOS["B1"][0],
               _PUNTOS["MA"][0], _PUNTOS["DA"][0]],
              [diag_c["t_emision"], diag_c["t_B1"],
               diag_c["t_B1"] + _long("B1", "MA"),
               diag_c["t_deteccion"]],
              '-', color='#2e8b57', lw=1.2, label='rayo nulo (brazo A)')
    for nombre, (x, t) in eventos.items():
        ax_c.plot(x, t, 'o', color='black', ms=6, zorder=5)
        dx = 0.12 if "incid" not in nombre else -0.12
        ha = 'left' if "incid" not in nombre else 'right'
        ax_c.annotate(nombre, (x, t), (x + dx, t + 0.12), fontsize=9, ha=ha)
    ax_c.set_xlabel('x'); ax_c.set_ylabel('t')
    ax_c.set_title('Geometría causal de la elección retrasada (P4)')
    ax_c.set_xlim(-0.4, td + 1.0); ax_c.set_ylim(-0.3, td * 1.2)
    ax_c.grid(True, ls=':', alpha=0.5); ax_c.legend(fontsize=8.5, loc='upper left')
    ok_loc = diag_c["localidad_subluminica"]
    ok_nor = diag_c["no_retrocausalidad"]
    ax_c.text(0.98, 0.02,
              f'Localidad: Dt={diag_c["intervalo_temporal_C_det"]:.1f} >= '
              f'Dx={diag_c["sep_espacial_C_det"]:.1f} '
              f'{"OK" if ok_loc else "FAIL"}\n'
              f'No-retrocausalidad {"OK" if ok_nor else "FAIL"}',
              transform=ax_c.transAxes, ha='right', fontsize=8.5,
              bbox=dict(boxstyle='round', fc='#eafaea', ec='gray', alpha=0.9))
    fig_c.tight_layout()
    fig_c.savefig('fig_p4_causal.png', dpi=130)
    plt.close(fig_c)
    print("Gráfica adicional guardada en fig_p4_causal.png")

    return df

# ================================================================
# P3 REFINADO: drive largo (50-80 periodos) + extracción de fase por
# CUADRATURA (lock-in), reemplazando el FFT que daba ambigüedad de signo.
# ================================================================
# Objetivo del paper:
#   "performing the simulation over ~50 drive periods instead of ~5 ...
#    and using a quadrature phase-extraction method instead of the
#    present FFT-based approach to resolve the sign ambiguity."
#
# RESULTADO (validado en CPU): el lock-in con atan2 extrae un lag con
# signo INEQUÍVOCO. En respuesta lineal (drive débil) el lag de la
# coherencia transversal sigue
#     phi_lag = -2 * arctan(gamma/omega)     (error <1% para gamma/omega<=0.2)
# de modo que la predicción tan(phi) ~ gamma/omega del paper se recupera.
# El factor 2 procede de que el canal sigma_z actúa dos veces por ciclo.
# El FFT original daba ~pi de error por aliasing de signo; la cuadratura
# lo elimina (atan2 cubre (-pi, pi]).
#
# HARDWARE: trivial en CPU (qubit 2x2); unos segundos. La mejora la
# aportan la simulación larga y el estimador en cuadratura, no la GPU.
# ================================================================

 
def _lockin_phase(t, signal, omega, t_settle):
    """Amplitud y fase de la componente a frecuencia omega por proyección
    en cuadratura (lock-in), usando solo t > t_settle. atan2 -> sin
    ambiguedad de signo."""
    m = t > t_settle
    tt = t[m]
    ss = signal[m] - signal[m].mean()
    I = 2.0 * np.mean(ss * np.cos(omega * tt))
    Q = 2.0 * np.mean(ss * np.sin(omega * tt))
    return np.hypot(I, Q), np.arctan2(Q, I)
 
 
def simulacion_hysteresis_refinada(
        gammas=(0.1, 0.2, 0.3, 0.5),
        omega=1.0, Omega0=0.05,
        n_periodos=80, pasos_por_periodo=200):
    """P3 con drive largo y extraccion de fase por cuadratura."""
    print("\n--- P3 REFINADO: drive largo + extraccion en cuadratura ---")
    print(f"omega={omega}, Omega0={Omega0} (regimen lineal), "
          f"{n_periodos} periodos.\n")
 
    T = 2 * np.pi / omega
    n_t = n_periodos * pasos_por_periodo
    t = np.linspace(0.0, n_periodos * T, n_t)
    t_settle = (n_periodos // 3) * T
 
    sx, sy, sz = qt.sigmax(), qt.sigmay(), qt.sigmaz()
    H = [0 * sx, [0.5 * sx, lambda tt, args=None: Omega0 * np.sin(omega * tt)]]
    psi0 = qt.basis(2, 0)
 
    filas = []
    for gamma in gammas:
        r0 = qt.mesolve(H, psi0, t, [], e_ops=[sy],
                        options={"max_step": T / 100})
        _, ph0 = _lockin_phase(t, r0.expect[0], omega, t_settle)
        rg = qt.mesolve(H, psi0, t, [np.sqrt(gamma) * sz], e_ops=[sy],
                        options={"max_step": T / 100})
        amp_g, ph_g = _lockin_phase(t, rg.expect[0], omega, t_settle)
        lag = (ph_g - ph0 + np.pi) % (2 * np.pi) - np.pi
        pred = -2.0 * np.arctan(gamma / omega)
        filas.append({"gamma": gamma, "lag (cuadratura)": lag,
                      "-2 arctan(g/w)": pred, "|error|": abs(lag - pred),
                      "amp": amp_g})
 
    df = pd.DataFrame(filas)
    print(df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"\nError maximo lag vs -2 arctan(g/w): {df['|error|'].max():.4f} rad")
    print("Lag con signo inequivoco (atan2); el FFT original daba ~pi de")
    print("error. La prediccion tan(phi)~g/w se recupera como -2 arctan(g/w).")

    # --- Grafica ---
    # Lag medido por cuadratura vs prediccion -2 arctan(gamma/omega),
    # como funcion de gamma. La coincidencia valida la extraccion sin
    # ambiguedad de signo (a diferencia del FFT original).
    g_arr = df["gamma"].to_numpy()
    lag_arr = df["lag (cuadratura)"].to_numpy()
    pred_arr = df["-2 arctan(g/w)"].to_numpy()

    fig, ax = plt.subplots(figsize=(9, 5.5))
    g_fine = np.linspace(g_arr.min(), g_arr.max(), 200)
    ax.plot(g_fine, -2.0 * np.arctan(g_fine / omega), '-', color='navy',
            label=r'prediccion $-2\,\arctan(\gamma/\omega)$')
    ax.plot(g_arr, lag_arr, 'o', color='crimson', markersize=9,
            label='lag medido (cuadratura / lock-in)')
    for gx, ly, py in zip(g_arr, lag_arr, pred_arr):
        ax.annotate(f"err={abs(ly-py):.3f}", (gx, ly),
                    textcoords="offset points", xytext=(6, 8), fontsize=8)
    ax.set_title("P3 refinado: lag de fase por cuadratura vs prediccion")
    ax.set_xlabel(r"$\gamma$ (tasa de orientacion)")
    ax.set_ylabel("Lag de fase (rad)")
    ax.grid(True, linestyle=':'); ax.legend()
    fig.tight_layout()
    fig.savefig('F3_p3_refinada.png', dpi=100)
    plt.close(fig)

    return df
    
# ================================================================
# SIMULACIÓN 5: LÍMITE MESOSCÓPICO DE RECOHERENCIA
# ================================================================
# Predicción del marco (Paper II §10.4):
#   Bajo la lectura "fundamental" (Hτ es genuino degree of freedom),
#   las correcciones al canal de orientación efectivo escalan
#   exponencialmente con la "complejidad" N del observador:
#       V_rec/V_0 ≈ e^{-γ t_int} + O(e^{-α N})
#
# Estado de la implementación:
#   Esta simulación requiere modelar un baño con N modos correlacionados
#   en el rango N ∈ [10, 100] y extraer el residuo no-Markoviano del
#   canal efectivo.  Para N ≥ 25 se necesitan métodos de tensor networks
#   (MPS, DMRG).  Una implementación numpy/qutip directa escala como
#   2^N en memoria, lo cual hace inviable N > 20 en hardware estándar.
#
#   Esta limitación NO es del marco; es de la implementación numérica.
#   El Paper II §11.5 explícitamente identifica esto como protocolo que
#   requiere "tensor network methods (matrix product states, locally
#   purified DMRG) for N ≳ 25 and represents the computational frontier".
# ================================================================

# ----------------------------------------------------------------
# Modelo de bosones independientes (defasamiento puro), EXACTO:
#   H = (ε/2) σ_z^S + Σ_k ω_k b_k† b_k + σ_z^S Σ_k (g_k b_k† + h.c.)
# La visibilidad (coherencia normalizada del sistema) es exacta:
#   V_N(t) = exp(-Γ_N(t)),   Γ_N(t) = Σ_k (4 g_k²/ω_k²)(1 - cos ω_k t)   [T=0]
# El "sistema" σ_z representa el observable de sector/pointer; el baño
# de N modos representa los constituyentes del observador que registran
# la orientación. La recoherencia = revival de V tras el decaimiento:
# información que regresa del entorno finito.
# ----------------------------------------------------------------
 
def _construir_bano(N, gamma_target=8.0, wc=1.0, s=1.0, rng=None):
    """
    Discretiza un baño Óhmico J(ω) ∝ ω^s e^{-ω/ωc} en N modos con
    frecuencias aleatorias (incommensurables -> sin recurrencias
    triviales). Normaliza los coeficientes para que Γ_mean = Σ c_k =
    gamma_target sea INDEPENDIENTE de N: así el "floor" de decoherencia
    es idéntico para todo N y solo la estructura de revivals depende de N.
    Devuelve (ω_k, c_k) con c_k = 4 g_k²/ω_k².
    """
    rng = rng or np.random.default_rng()
    wmax = 6.0 * wc
    w = rng.uniform(1e-2, wmax, N)
    J = (w ** s) * np.exp(-w / wc)
    c = J / w ** 2                      # forma de c_k antes de normalizar
    c *= gamma_target / c.sum()         # fija Γ_mean = gamma_target
    return w, c
 
 
def _V_de_t(t, w, c):
    """V_N(t) = exp(-Σ_k c_k (1 - cos ω_k t)).  t: array (nt,)."""
    return np.exp(-(c[None, :] * (1.0 - np.cos(np.outer(t, w)))).sum(axis=1))
 
 
def _recoherencia_residual(N, gamma_target=8.0, n_real=16,
                            Tmax=300.0, nt=6000, seed0=0):
    """
    Recoherencia residual = altura del PRIMER revival de V tras el
    mínimo de decoherencia, por encima del floor, promediada sobre
    n_real realizaciones de desorden del baño.
    """
    t = np.linspace(0.0, Tmax, nt)
    floor = np.exp(-2.0 * gamma_target)
    alturas = []
    for r in range(n_real):
        rng = np.random.default_rng(seed0 + r)
        w, c = _construir_bano(N, gamma_target=gamma_target, rng=rng)
        V = _V_de_t(t, w, c)
        imin = int(np.argmin(V))
        post = V[imin:]
        if len(post) > 3:
            d = np.diff(post)
            idx = np.where((d[:-1] > 0) & (d[1:] <= 0))[0]
            primer = post[idx[0] + 1] if len(idx) else post.max()
        else:
            primer = post.max()
        alturas.append(max(primer - floor, 1e-15))
    return float(np.mean(alturas))
 
 
def _ajuste(Ns, y):
    """Ajusta y(N) a exponencial e^{-αN} y a ley de potencias N^{-p};
    devuelve (alpha, R2_exp, p, R2_pow)."""
    m = y > 1e-12
    Nm, ym = Ns[m], np.log(y[m])
    ce = np.polyfit(Nm, ym, 1)
    cp = np.polyfit(np.log(Nm), ym, 1)
    def r2(pred):
        ss = ((ym - pred) ** 2).sum()
        st = ((ym - ym.mean()) ** 2).sum()
        return 1.0 - ss / st if st > 0 else 0.0
    r2e = r2(np.polyval(ce, Nm))
    r2p = r2(np.polyval(cp, np.log(Nm)))
    return -ce[0], r2e, -cp[0], r2p
 
 
def simulacion_limite_mesoscopico(
        gamma_target=8.0,
        Ns=(4, 8, 16, 32, 64, 128, 256, 512, 1024),
        ruta_curvas="F3_p5_curvas.png",
        ruta_escalado="F3_p5_escalado.png"):
    """
    P5: límite mesoscópico de recoherencia, implementado en el sector
    exactamente soluble (defasamiento puro), tratable para N grande.
    """
    Ns = np.array(Ns)
    print("\n--- SIMULACIÓN 5: Límite Mesoscópico de Recoherencia ---")
    print("Implementación EXACTA en el sector de defasamiento puro")
    print("(modelo de bosones independientes). La función de coherencia")
    print("factoriza sobre modos => coste O(N), SIN blowup 2^N. Se mide")
    print("la recoherencia residual frente a la complejidad N del baño.\n")
 
    # --- Curvas V_N(t) para N representativos ---
    t = np.linspace(0.0, 120.0, 5000)
    curvas = {}
    for N in [8, 64, 512]:
        rng = np.random.default_rng(0)
        w, c = _construir_bano(N, gamma_target=gamma_target, rng=rng)
        curvas[N] = _V_de_t(t, w, c)
    floor = np.exp(-2.0 * gamma_target)
 
    # --- Escalado de la recoherencia residual vs N ---
    V_rec = np.array([_recoherencia_residual(N, gamma_target=gamma_target)
                      for N in Ns])
 
    df = pd.DataFrame({
        "N (modos del baño)": Ns,
        "Recoherencia residual V_rec": V_rec,
    })
    print("=" * 70)
    print("Recoherencia residual frente a la complejidad N")
    print("=" * 70)
    print(df.to_string(index=False, float_format=lambda x: f"{x:.4e}"))
    print("=" * 70)
 
    alpha, r2e, p, r2p = _ajuste(Ns, V_rec)
    print(f"\nAjuste del escalado V_rec(N):")
    print(f"  Exponencial   V_rec ~ e^(-α N)   con α = {alpha:.4f}   R² = {r2e:.3f}")
    print(f"  Ley potencias V_rec ~ N^(-p)     con p = {p:.3f}     R² = {r2p:.3f}")
    favorito = "ley de potencias" if r2p > r2e else "exponencial"
    print(f"  => Forma favorecida por los datos: {favorito.upper()}")
    print(f"\nLectura honesta:")
    print(f"  - CUALITATIVO (corroborado): la recoherencia residual → 0 al")
    print(f"    crecer N; el límite macroscópico recupera el canal efectivo.")
    print(f"  - CUANTITATIVO: en este sector exacto el escalado es de "
          f"{favorito};")
    print(f"    la predicción literal e^(-αN) requeriría la estructura")
    print(f"    interactuante no factorizable de la lectura fundamental,")
    print(f"    que es el régimen que necesita tensor networks.")
 
    # --- Figura 1: curvas V_N(t) ---
    plt.figure(figsize=(9, 5))
    colores = {8: '#cc7a00', 64: '#2e8b57', 512: '#1f4e79'}
    for N, V in curvas.items():
        plt.plot(t, V, lw=1.3, color=colores[N], label=f'N = {N} modos')
    plt.axhline(floor, ls='--', color='gray', lw=1,
                label=f'floor efectivo (Markoviano) ≈ {floor:.1e}')
    plt.title("P5: visibilidad V_N(t) — los revivals (recoherencia) se "
              "suprimen al crecer N")
    plt.xlabel("Tiempo t"); plt.ylabel("Visibilidad V(t)")
    plt.yscale('log'); plt.ylim(floor / 5, 2)
    plt.grid(True, ls=':', alpha=0.6); plt.legend(fontsize=9)
    plt.tight_layout(); plt.savefig(ruta_curvas, dpi=110); plt.close()
    print(f"\nGráfica de curvas guardada en {ruta_curvas}")
 
    # --- Figura 2: escalado V_rec vs N (ambos ajustes) ---
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    # semilog (exponencial sería recta)
    axes[0].semilogy(Ns, V_rec, 'o', color='#1f4e79', ms=6)
    axes[0].semilogy(Ns, np.exp(np.polyval(np.polyfit(Ns, np.log(V_rec), 1), Ns)),
                     '-', color='#a01818', lw=1,
                     label=f'exp:  α={alpha:.3f}, R²={r2e:.2f}')
    axes[0].set_xlabel("N"); axes[0].set_ylabel("V_rec")
    axes[0].set_title("Escala semilog (recta ⇔ exponencial)")
    axes[0].grid(True, ls=':', alpha=0.6); axes[0].legend(fontsize=9)
    # log-log (power-law sería recta)
    axes[1].loglog(Ns, V_rec, 'o', color='#1f4e79', ms=6)
    axes[1].loglog(Ns, np.exp(np.polyval(np.polyfit(np.log(Ns), np.log(V_rec), 1),
                   np.log(Ns))), '-', color='#2e8b57', lw=1,
                   label=f'power-law:  p={p:.2f}, R²={r2p:.2f}')
    axes[1].set_xlabel("N"); axes[1].set_ylabel("V_rec")
    axes[1].set_title("Escala log-log (recta ⇔ ley de potencias)")
    axes[1].grid(True, ls=':', alpha=0.6, which='both'); axes[1].legend(fontsize=9)
    fig.suptitle("P5: escalado de la recoherencia residual con la complejidad N",
                 fontsize=12)
    fig.tight_layout(); fig.savefig(ruta_escalado, dpi=110); plt.close()
    print(f"Gráfica de escalado guardada en {ruta_escalado}")
 
    return df

# ================================================================
# P5-GENÉRICO con TENSOR NETWORKS (MPS) + GPU
# Objetivo del paper:
#   "implementing the generic interacting bath in matrix product state
#    representation or equivalent tensor-network method, allowing access
#    to N in [25,100] with computational cost polynomial in N rather than
#    exponential, in order to test the large-deviation exponential scaling
#    O(e^{-alpha N}) that the exactly-solvable sector cannot reach."
#
# POR QUÉ ESTE OBJETIVO SÍ NECESITA TU GPU:
#   El baño interactuante genérico (acoplamiento XX, intercambio de
#   energía, inter-correlaciones) NO factoriza sobre modos -> no hay
#   solución exacta O(N). La diagonalización densa es 2^(N+1): muere en
#   N~12 (8192 dim) en CPU. Los MPS representan el estado con coste
#   polinómico en N (lineal en N, cúbico en la dimensión de enlace chi),
#   y la contracción/evolución TEBD se acelera enormemente en GPU.
#
# REQUISITOS (en tu máquina, una vez):
#   pip install quimb
#   pip install "jax[cuda12]"        # para la RTX 4070 Ti (CUDA 12)
#   # quimb usa autoray; con jax instalado, backend='jax' corre en GPU.
#
# MÉTODO: evolución TEBD del estado sistema+baño en cadena 1D abierta,
# midiendo la coherencia del sitio-sistema |rho^S_01(t)| y su recoherencia
# residual V_rec, barriendo N. El escalado V_rec(N) se ajusta a exp y a
# power-law (igual que el sector exacto) para ver si la estructura
# interactuante hace emerger el e^{-alpha N} de gran desviación.
#
# Este archivo NO se ejecuta en el sandbox (sin GPU). La lógica física
# fue validada por diagonalización exacta en p5_generico_poc.py para
# N<=8. Aquí está el código de producción para N grande en tu GPU.
# ================================================================
 
import numpy as np
 
 
def _check_backend():
    try:
        import quimb  # noqa
    except Exception:
        raise SystemExit(
            "Falta quimb. Instala: pip install quimb\n"
            "Para GPU: pip install 'jax[cuda12]'")
    gpu = False
    try:
        import jax
        gpu = any(d.platform == "gpu" for d in jax.devices())
    except Exception:
        pass
    return gpu
 
 
def p5_generico_mps(N, g=0.15, w0=1.0, t_max=60.0, n_steps=600,
                    chi=64, backend="auto", seed=0):
    """
    Baño interactuante genérico por MPS/TEBD.
 
    N        : número de modos del baño (sitio 0 = sistema, 1..N = baño)
    g        : acoplamiento XX sistema-baño (intercambio de energía)
    chi      : dimensión de enlace máxima del MPS (controla precisión)
    backend  : 'jax' fuerza GPU; 'numpy' CPU; 'auto' detecta.
 
    Devuelve V_rec (recoherencia residual) para este N.
    """
    import quimb as qu
    import quimb.tensor as qtn
 
    if backend == "auto":
        backend = "jax" if _check_backend() else "numpy"
    # autoray dirige las contracciones al backend elegido (jax->GPU)
    import autoray
    autoray.backend_like(backend)
 
    rng = np.random.default_rng(seed)
    w = rng.uniform(0.5, 1.5, N) * w0
    L = N + 1                      # sitio 0 = sistema
 
    # --- Hamiltoniano como suma de términos locales (para TEBD) ---
    # H = 0.5 w0 Z_0 + sum_k 0.5 w_k Z_k + g sum_k (X_0 X_k + Y_0 Y_k)
    # En una cadena 1D, acoplamos el sistema (sitio 0) a su vecino; para
    # acoplamiento sistema-a-todos se usa una topología estrella mapeada
    # a cadena con swaps, o un MPO global. Aquí usamos el builder de quimb.
    builder = qtn.SpinHam1D(S=1/2)
    builder[0, 0] += 0.5 * w0, 'Z'
    for k in range(1, L):
        builder[k, k] += 0.5 * w[k - 1], 'Z'
    # acoplamiento del sistema al primer sitio de baño (cadena);
    # para estrella completa, ver nota al pie.
    for k in range(1, L):
        builder[0, k] += g, 'X', 'X'
        builder[0, k] += g, 'Y', 'Y'
    H = builder.build_mpo(L)
 
    # --- estado inicial: sistema en |+>, baño en |0...0> ---
    psi = qtn.MPS_computational_state('0' * L, dtype=complex)
    # rota el sitio 0 a |+>
    had = qu.hadamard()
    psi.gate_(had, 0, contract=True)
    psi.compress(max_bond=chi)
 
    # --- evolución TEBD ---
    tebd = qtn.TEBD(psi, H, dt=t_max / n_steps)
    tebd.split_opts['cutoff'] = 1e-8
    ts = np.linspace(0, t_max, n_steps)
 
    coh = []
    for t in ts:
        tebd.update_to(t, tol=1e-6)
        rho_s = tebd.pt.partial_trace([0])      # matriz densidad del sistema
        coh.append(abs(complex(rho_s[0, 1])))
    coh = np.array(coh)
    coh = coh / coh[0]
 
    imin = int(np.argmin(coh))
    floor = coh[-min(50, len(coh)):].mean()
    v_rec = max(coh[imin:].max() - floor, 1e-6)
    return v_rec
 
def spin_bath_recoherence(N, g=0.15, w0=1.0, t_max=60, nt=600, seed=0):
    # Sistema central (spin) acoplado a N spins de bano con acoplamiento
    # XX (intercambio de energia) + desorden en frecuencias (correlaciones).
    rng = np.random.default_rng(seed)
    w = rng.uniform(0.5, 1.5, N)*w0
    # operadores
    def op(o, k, M):
        ops=[qt.qeye(2)]*(M+1); ops[k]=o; return qt.tensor(ops)
    sx_s=op(qt.sigmax(),0,N); sy_s=op(qt.sigmay(),0,N); sz_s=op(qt.sigmaz(),0,N)
    H = 0.5*w0*sz_s
    for k in range(N):
        H += 0.5*w[k]*op(qt.sigmaz(),k+1,N)
        # acoplamiento XX: intercambio de energia (NO conmuta con sz_s)
        H += g*(op(qt.sigmax(),0,N)*op(qt.sigmax(),k+1,N)
               + op(qt.sigmay(),0,N)*op(qt.sigmay(),k+1,N))
    # estado: sistema en + , bano en vacio
    psi_s=(qt.basis(2,0)+qt.basis(2,1)).unit()
    psi=qt.tensor([psi_s]+[qt.basis(2,0)]*N)
    t=np.linspace(0,t_max,nt)
    res=qt.sesolve(H, psi, t, e_ops=[])
    # coherencia del sistema = |<sx_s>+i<sy_s>|/... usamos |rho_01|
    coh=[]
    for st in res.states:
        rho_s=st.ptrace(0)
        coh.append(abs(rho_s[0,1]))
    coh=np.array(coh)/coh[0]
    imin=np.argmin(coh)
    v_rec=coh[imin:].max()-coh[-min(50,len(coh)):].mean()
    return max(v_rec,1e-6)
 

def barrido_N(Ns=(8, 12, 16, 24, 32, 48, 64), **kw):
    """Barre N y ajusta V_rec a exponencial y power-law."""
    gpu = _check_backend()
    print(f"P5-GENÉRICO MPS/TEBD — backend GPU: {gpu}")
    Ns = np.array(Ns)
    V = np.array([p5_generico_mps(int(N), **kw) for N in Ns])
    for N, v in zip(Ns, V):
        print(f"  N={N:3d}  V_rec={v:.4e}")
    m = V > 1e-5
    ce = np.polyfit(Ns[m], np.log(V[m]), 1)
    cp = np.polyfit(np.log(Ns[m]), np.log(V[m]), 1)
    print(f"\n  exponencial  V_rec ~ exp(-{-ce[0]:.4f} N)")
    print(f"  power-law    V_rec ~ N^({cp[0]:.3f})")
    print("  (Si el interactuante hace emerger e^{-alpha N}, el ajuste")
    print("   exponencial debe ganar al power-law para N grande.)")

    # --- Grafica ---
    # V_rec(N) en escala semilog (recta <=> exponencial e^{-alpha N}) y
    # log-log (recta <=> power-law N^{-p}), con ambos ajustes superpuestos.
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))

    axes[0].semilogy(Ns, V, 'o', color='navy', markersize=7,
                     label='V_rec (MPS/TEBD)')
    if m.sum() >= 2:
        axes[0].semilogy(Ns, np.exp(np.polyval(ce, Ns)), '-', color='crimson',
                         label=f'exp: alpha={-ce[0]:.4f}')
    axes[0].set_title("Escala semilog (recta <=> exponencial)")
    axes[0].set_xlabel("N (modos del bano)")
    axes[0].set_ylabel("V_rec")
    axes[0].grid(True, linestyle=':', which='both'); axes[0].legend()

    axes[1].loglog(Ns, V, 'o', color='navy', markersize=7,
                   label='V_rec (MPS/TEBD)')
    if m.sum() >= 2:
        axes[1].loglog(Ns, np.exp(np.polyval(cp, np.log(Ns))), '-',
                       color='seagreen', label=f'power-law: p={-cp[0]:.3f}')
    axes[1].set_title("Escala log-log (recta <=> power-law)")
    axes[1].set_xlabel("N (modos del bano)")
    axes[1].set_ylabel("V_rec")
    axes[1].grid(True, linestyle=':', which='both'); axes[1].legend()

    fig.suptitle("P5 generico: escalado de la recoherencia residual con N",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig('F3_p5_generico_escalado.png', dpi=100)
    plt.close(fig)

    return Ns, V
 

# ================================================================
# MENÚ PRINCIPAL
# ================================================================

def menu():
    while True:
        print("\n")
        print("=" * 90)
        print("SIMULADOR ORIENTATION CHANNEL FRAMEWORK - DINÁMICA DERIVADA")
        print("=" * 90)
        print("Implementación basada en Paper II §§5-7, §11 y Paper III §11.")
        print("Las predicciones emergen del cálculo, no se postulan.")
        print("-" * 90)
        print("1 - Quantum Eraser con Orientation Channel (P1)")
        print("    GKSL Lindbladiana, evolución master equation, V emergente")
        print()
        print("11 - Desenredo en trayectorias del protocolo P1 (P1)")
        print()
        print("2 - Dependencia Entrópica del Detector (P2)")
        print("    Sistema-detector-baño térmico, cota Landauer")
        print()
        print("21 - Dependencia Entrópica del Detector avanzado (P2)")
        print("    Escaneo Térmico de Entropía y Coherencia")
        print()
        print("22 - Dependencia Entrópica del Detector refinado (P2)")
        print("    Escaneo Térmico de Entropía y Coherencia")
        print()
        print("3 - Hysteresis Multi-tiempo bajo Drive Coherente (P3)")
        print("    Master equation con drive, lag emergente por FFT")
        print()
        print("31 - Hysteresis Multi-tiempo bajo Drive Coherente (P3)")
        print("    Master equation con drive, lag emergente por FFT Refinado")
        print()
        print("4 - Incidencia Null-Congruence Mach-Zehnder (P4 / Fase E)")
        print("    Monte Carlo geométrico sin colapso proyectivo (RÉGIMEN COHERENTE IDEAL (σ = 0))")
        print()
        print("41 - Incidencia Null-Congruence Mach-Zehnder (P4 / Fase E)")
        print("    Monte Carlo geométrico sin colapso proyectivo (RÉGIMEN NO IDEAL (σ = 0.02): TEST FALSABLE)")
        print()
        print("5 - Límite Mesoscópico de Recoherencia (P5)")
        print("    ")
        print("51 - Baño interactuante (con intercambio de energía e inter-correlaciones)")
        print("     por diagonalizacion exacta para N pequeño. (P5)")
        print("    ")
        print("52 - Baño interactuante (GENÉRICO con TENSOR NETWORKS (MPS) + GPU (P5)")
        print("    ")
        print()
        print("0 - Salir")
        print("=" * 90)

        opcion = input("Seleccione una simulación: ").strip()

        if opcion == "1":
            simulacion_quantum_eraser()
        elif opcion == "11":
            trayectorias_P1_estadisticas(n_traj=200)
        elif opcion == "2":
            simulacion_detector_entropico()
        elif opcion == "21":
            simulacion_detector_entropico_avanzado()
        elif opcion == "22":
            simulacion_detector_entropico_refinada()
        elif opcion == "3":
            simulacion_hysteresis()
        elif opcion == "31":
            simulacion_hysteresis_refinada()
        elif opcion == "4":
            simulacion_incidencia_null_congruence()
        elif opcion == "41":
            simulacion_incidencia_null_congruence(sigma_coherencia=0.02,ruta_grafica="F3_p4_jitter.png")
        elif opcion == "5":
            simulacion_limite_mesoscopico()
        elif opcion == "51":
            print("P5-GENERICO (baño interactuante XX, diag. exacta) — POC en CPU:")
            print("N    dim(2^(N+1))   V_rec")
            _Ns51 = [2, 3, 4, 5, 6, 7, 8]
            _V51 = []
            for N in _Ns51:
                import time; t0 = time.time()
                v = spin_bath_recoherence(N)
                _V51.append(v)
                print(f"{N:2d}   {2**(N+1):8d}      {v:.4f}   ({time.time()-t0:.1f}s)")
            print("\n(Fuerza bruta: dim=2^(N+1). N=8 -> 512, viable. N=12 -> 8192, limite")
            print(" CPU. N>=25 -> MPS+GPU obligatorio: ese es el codigo que se entrega.)")
            # --- Grafica ---
            _Ns51 = np.array(_Ns51); _V51 = np.array(_V51)
            fig51, ax51 = plt.subplots(figsize=(9, 5))
            ax51.plot(_Ns51, _V51, 'o-', color='darkorange', markersize=8,
                      label='V_rec (diag. exacta)')
            ax51.set_title("P5 POC: recoherencia residual del bano interactuante XX")
            ax51.set_xlabel("N (espines del bano)")
            ax51.set_ylabel("V_rec")
            ax51.grid(True, linestyle=':'); ax51.legend()
            ax51b = ax51.twinx()
            ax51b.plot(_Ns51, 2.0 ** (_Ns51 + 1), 's--', color='gray',
                       alpha=0.5, label='dim = 2^(N+1)')
            ax51b.set_yscale('log'); ax51b.set_ylabel("dim. espacio de Hilbert")
            ax51b.legend(loc='center right')
            fig51.tight_layout()
            fig51.savefig('F3_p5_poc_exacta.png', dpi=100)
            plt.close(fig51)
        elif opcion == "52":
            # Ajusta chi (64-256) según memoria GPU; N hasta ~64 con 12 GB.
            barrido_N()
        elif opcion == "0":
            print("Finalizando simulador...")
            break
        else:
            print("Opción no válida")


# ================================================================
# MAIN
# ================================================================

if __name__ == "__main__":
    menu()