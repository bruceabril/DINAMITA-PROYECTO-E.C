"""
Proyecto "Dinamita" - Modelo Térmico Dinámico de un Procesador
Asignatura: Ecuaciones Diferenciales


Ecuación diferencial modelada (Ley de Enfriamiento de Newton):
    dT/dt = -k * (T - Ta) + Q(t)

Este script entrega los 3 resultados requeridos:
    1) GRÁFICAS -> temperatura simulada vs. datos reales
    2) CÁLCULO DE K -> constante de disipación térmica estimada
    3) MÉTRICA DE ERROR -> RMSE y R^2 entre el modelo y los datos

Cómo correrlo:
    python modelo_termico.py
    (usa datos sintéticos de ejemplo si no se le pasa un CSV real)

    python modelo_termico.py mis_datos.csv
    (se usan los log de HWInfo/Open Hardware Monitor propios, con columnas
     'tiempo_s' y 'temperatura_c')
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit


# 1) SOLUCIÓN ANALÍTICA (para carga/potencia Q constante)
def solucion_analitica(t, temp_inicial, temp_ambiente, calor_generado, k):
    """
    T(t) = (Ta + Q0/k) + (T0 - Ta - Q0/k) * e^(-k t)

    temp_inicial     -> T0, temperatura inicial del CPU (°C)
    temp_ambiente    -> Ta, temperatura del entorno (°C)
    calor_generado   -> Q0, generación de calor por carga de trabajo
    k                -> constante de disipación térmica
    """
    temp_equilibrio = temp_ambiente + calor_generado / k
    return temp_equilibrio + (temp_inicial - temp_ambiente - calor_generado / k) * np.exp(-k * t)


 
# 2) SOLUCIÓN NUMÉRICA (RK4/RK45) - válida también con carga variable

def simular_rk4(tiempos_evaluar, temp_inicial, temp_ambiente, k, calor_generado):
    """Resuelve la EDO numéricamente usando Runge-Kutta (RK45 de scipy)."""

    def ecuacion_diferencial(t, temperatura):
        return -k * (temperatura - temp_ambiente) + calor_generado

    solucion = solve_ivp(
        ecuacion_diferencial,
        t_span=(tiempos_evaluar[0], tiempos_evaluar[-1]),
        y0=[temp_inicial],
        method="RK45",
        t_eval=tiempos_evaluar,
        max_step=0.5,
    )
    return solucion.y[0]


# 3) CÁLCULO DE K - ajuste por mínimos cuadrados a los datos reales
def calcular_constante_k(tiempo_datos, temperatura_datos, temp_inicial, temp_ambiente):
    """
    Ajusta k y Q0 (calor generado) a los datos reales medidos,
    usando el método de mínimos cuadrados (curve_fit).
    """

    def modelo(t, calor_generado, k):
        return solucion_analitica(t, temp_inicial, temp_ambiente, calor_generado, k)

    valores_iniciales = [5.0, 0.05]  # punto de arranque para el optimizador
    parametros_ajustados, _ = curve_fit(
        modelo, tiempo_datos, temperatura_datos, p0=valores_iniciales, maxfev=5000
    )
    calor_ajustado, k_ajustado = parametros_ajustados
    return calor_ajustado, k_ajustado


# 4) MÉTRICAS DE ERROR - RMSE y R^2
def calcular_metricas_error(temperatura_real, temperatura_predicha):
    """Calcula RMSE (error cuadrático medio), MAE (error absoluto medio) y R^2."""
    error = temperatura_real - temperatura_predicha
    rmse = np.sqrt(np.mean(error ** 2))
    mae = np.mean(np.abs(error))

    suma_residuos = np.sum(error ** 2)
    suma_total = np.sum((temperatura_real - np.mean(temperatura_real)) ** 2)
    r2 = 1 - suma_residuos / suma_total

    return rmse, mae, r2



# 5) CARGA DE DATOS REALES (opcional, si tienes un CSV de HWInfo/OHM)

def cargar_datos_csv(ruta_archivo):
    """
    Lee un CSV con columnas 'tiempo_s' y 'temperatura_c'
    (exportado, por ejemplo, de HWInfo o Open Hardware Monitor).
    """
    datos = pd.read_csv(ruta_archivo)
    return datos["tiempo_s"].values, datos["temperatura_c"].values



# GENERAR DATOS DE PRUEBA (mientras no tengan el CSV real del sensor)
def generar_datos_sinteticos(temp_inicial, temp_ambiente):
    """Simula datos 'medidos' con ruido, para probar el flujo completo."""
    k_real, calor_real = 0.08, 6.0
    tiempo = np.linspace(0, 120, 60)
    temperatura = solucion_analitica(tiempo, temp_inicial, temp_ambiente, calor_real, k_real)
    temperatura += np.random.normal(0, 0.4, size=tiempo.shape)  # ruido de sensor
    return tiempo, temperatura


# GRAFICAR RESULTADOS
def graficar_resultados(tiempo_real, temperatura_real, tiempo_simulado, temperatura_simulada, k_estimado):
    plt.figure(figsize=(9, 5))
    plt.plot(tiempo_real, temperatura_real, "o", markersize=3, alpha=0.6, label="Datos (medidos)")
    plt.plot(
        tiempo_simulado, temperatura_simulada, "-", linewidth=2,
        label=f"Modelo RK4 (k={k_estimado:.4f})"
    )
    plt.xlabel("Tiempo (s)")
    plt.ylabel("Temperatura (°C)")
    plt.title("Modelo térmico del CPU: simulación vs. datos")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("comparacion_modelo_real.png", dpi=150)
    print("\nGráfica guardada como 'comparacion_modelo_real.png'")
    plt.show()


# PROGRAMA PRINCIPAL
def principal():
    temp_ambiente = 25.0   # temperatura ambiente supuesta (°C)
    temp_inicial = 30.0    # temperatura inicial del CPU (°C)

    # --- Obtener datos: reales si se pasa un CSV, si no, sintéticos ---
    if len(sys.argv) > 1:
        tiempo_datos, temperatura_datos = cargar_datos_csv(sys.argv[1])
        print(f"Datos cargados desde: {sys.argv[1]}")
    else:
        print("No se dio un CSV -> usando datos sintéticos de ejemplo.")
        tiempo_datos, temperatura_datos = generar_datos_sinteticos(temp_inicial, temp_ambiente)

    # --- (2) CÁLCULO DE K ---
    calor_ajustado, k_ajustado = calcular_constante_k(
        tiempo_datos, temperatura_datos, temp_inicial, temp_ambiente
    )
    print("\n=== CÁLCULO DE LA CONSTANTE K ===")
    print(f"k estimado  = {k_ajustado:.5f} 1/s")
    print(f"Q0 estimado = {calor_ajustado:.5f} °C/s")

    # --- Generar la curva simulada con el k ya calculado (RK4) ---
    temperatura_simulada = simular_rk4(
        tiempo_datos, temp_inicial, temp_ambiente, k_ajustado, calor_ajustado
    )

    # --- (3) MÉTRICA DE ERROR ---
    rmse, mae, r2 = calcular_metricas_error(temperatura_datos, temperatura_simulada)
    print("\n=== MÉTRICA DE ERROR (modelo vs. datos) ===")
    print(f"RMSE = {rmse:.4f} °C")
    print(f"MAE  = {mae:.4f} °C")
    print(f"R^2  = {r2:.4f}")

    # --- (1) GRÁFICAS ---
    graficar_resultados(tiempo_datos, temperatura_datos, tiempo_datos, temperatura_simulada, k_ajustado)


if __name__ == "__main__":
    principal()