"""
Proyecto "Dinamita" - Modelo Térmico Dinámico de un Procesador
Asignatura: Ecuaciones Diferenciales
"""

import json
import os
import re
import sys
import time
import urllib.request
import warnings

try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp
    from scipy.optimize import curve_fit
except ImportError as e:
    print("\nERROR: faltan librerías necesarias.")
    print("Instálalas con:")
    print("    python -m pip install numpy pandas scipy matplotlib")
    print(f"\nDetalle: {e}")
    sys.exit(1)


# --------------------------------------------
# 1) SOLUCIÓN ANALÍTICA
# --------------------------------------------
def solucion_analitica(t, temp_inicial, temp_ambiente, calor_generado, k):
    """
    T(t) = (Ta + Q0/k) + (T0 - Ta - Q0/k) * e^(-k t)
    """
    temp_equilibrio = temp_ambiente + calor_generado / k
    return temp_equilibrio + (temp_inicial - temp_ambiente - calor_generado / k) * np.exp(-k * t)


# --------------------------------------------
# 2) SOLUCIÓN NUMÉRICA (RK45)
# --------------------------------------------
def simular_rk4(tiempos_evaluar, temp_inicial, temp_ambiente, k, calor_generado):
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

    if not solucion.success:
        raise RuntimeError(f"Error al integrar RK45: {solucion.message}")

    return solucion.y[0]


# --------------------------------------------
# 3) CÁLCULO DE K
# --------------------------------------------
def calcular_constante_k(tiempo_datos, temperatura_datos, temp_inicial, temp_ambiente):
    def modelo(t, calor_generado, k):
        return solucion_analitica(t, temp_inicial, temp_ambiente, calor_generado, k)

    valores_iniciales = [5.0, 0.05]
    try:
        parametros_ajustados, _ = curve_fit(
            modelo, tiempo_datos, temperatura_datos, p0=valores_iniciales, maxfev=5000
        )
    except Exception as e:
        raise RuntimeError(f"No se pudo ajustar el modelo: {e}")

    calor_ajustado, k_ajustado = parametros_ajustados
    return calor_ajustado, k_ajustado


# --------------------------------------------
# 4) MÉTRICAS DE ERROR
# --------------------------------------------
def calcular_metricas_error(temperatura_real, temperatura_predicha):
    error = temperatura_real - temperatura_predicha
    rmse = np.sqrt(np.mean(error ** 2))
    mae = np.mean(np.abs(error))

    suma_residuos = np.sum(error ** 2)
    suma_total = np.sum((temperatura_real - np.mean(temperatura_real)) ** 2)

    if suma_total == 0:
        r2 = 1.0
    else:
        r2 = 1 - suma_residuos / suma_total

    return rmse, mae, r2


# --------------------------------------------
# 5) CARGA DE CSV
# --------------------------------------------
def cargar_datos_csv(ruta_archivo):
    if not os.path.exists(ruta_archivo):
        raise FileNotFoundError(f"No existe el archivo: {ruta_archivo}")

    try:
        datos = pd.read_csv(ruta_archivo)
    except Exception as e:
        raise ValueError(f"No se pudo leer el CSV: {e}")

    columnas_esperadas = ["tiempo_s", "temperatura_c"]

    faltantes = [c for c in columnas_esperadas if c not in datos.columns]
    if faltantes:
        raise ValueError(
            f"El CSV no tiene las columnas requeridas: {columnas_esperadas}. "
            f"Columnas encontradas: {list(datos.columns)}"
        )

    tiempo = pd.to_numeric(datos["tiempo_s"], errors="raise").to_numpy()
    temperatura = pd.to_numeric(datos["temperatura_c"], errors="raise").to_numpy()

    if len(tiempo) < 2:
        raise ValueError("El CSV debe contener al menos 2 puntos.")

    return tiempo, temperatura


# --------------------------------------------
# 6) DATOS SINTÉTICOS
# --------------------------------------------
def generar_datos_sinteticos(temp_inicial, temp_ambiente):
    k_real, calor_real = 0.08, 6.0
    tiempo = np.linspace(0, 120, 60)
    temperatura = solucion_analitica(tiempo, temp_inicial, temp_ambiente, calor_real, k_real)
    temperatura += np.random.normal(0, 0.4, size=tiempo.shape)
    return tiempo, temperatura


# --------------------------------------------
# 7) GRAFICAR RESULTADOS
# --------------------------------------------
def graficar_resultados(tiempo_real, temperatura_real, tiempo_simulado, temperatura_simulada, k_estimado):
    plt.figure(figsize=(9, 5))
    plt.plot(tiempo_real, temperatura_real, "o", markersize=3, alpha=0.6, label="Datos (medidos)")
    plt.plot(
        tiempo_simulado,
        temperatura_simulada,
        "-",
        linewidth=2,
        label=f"Modelo RK4 (k={k_estimado:.4f})",
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
    plt.close()


def leer_temperatura_ohm(url):
    with urllib.request.urlopen(url, timeout=3) as respuesta:
        datos = json.loads(respuesta.read().decode("utf-8"))

    valores = []

    def buscar_sensor(nodo):
        if isinstance(nodo, dict):
            nombre = str(nodo.get("Text", "")).strip().casefold()
            if nombre == "cpu package" and nodo.get("Value") is not None:
                valores.append(nodo["Value"])
            for hijo in nodo.values():
                buscar_sensor(hijo)
        elif isinstance(nodo, list):
            for hijo in nodo:
                buscar_sensor(hijo)

    buscar_sensor(datos)
    if not valores:
        raise ValueError("No se encontró el sensor 'CPU Package' en los datos de Open Hardware Monitor.")

    coincidencia = re.search(r"[-+]?(?:\d+(?:[.,]\d*)?|[.,]\d+)", str(valores[0]))
    if coincidencia is None:
        raise ValueError(f"El valor del sensor 'CPU Package' no es una temperatura válida: {valores[0]}")

    return float(coincidencia.group().replace(",", "."))


def calibrar_parametros_en_vivo(tiempos, temperaturas, temp_ambiente):
    if len(tiempos) < 10 or tiempos[-1] - tiempos[0] < 8:
        return None
    if np.ptp(temperaturas) < 0.3:
        return None

    tiempo_inicio = tiempos[0]
    temp_inicial = temperaturas[0]

    def modelo(t, calor_generado, k):
        return solucion_analitica(
            t - tiempo_inicio, temp_inicial, temp_ambiente, calor_generado, k
        )

    calor_inicial = max((np.mean(temperaturas) - temp_ambiente) * 0.05, 0.1)
    parametros, _ = curve_fit(
        modelo,
        tiempos,
        temperaturas,
        p0=(calor_inicial, 0.05),
        bounds=((0, 1e-6), (np.inf, 2.0)),
        maxfev=10000,
    )
    return parametros


def graficar_temperatura_en_vivo(url, temp_ambiente):
    from matplotlib.animation import FuncAnimation

    temperatura_inicial = leer_temperatura_ohm(url)
    inicio = time.monotonic()
    tiempos = [0.0]
    temperaturas = [temperatura_inicial]
    ultimo_error = [None]

    figura, eje = plt.subplots(figsize=(9, 5))
    linea_datos, = eje.plot(
        tiempos, temperaturas, "-", marker="o", markersize=3, label="CPU Package medido"
    )
    linea_modelo, = eje.plot([], [], "--", linewidth=2, label="Modelo analítico")
    estado = eje.set_title("Esperando lecturas para calibrar el modelo...")
    eje.set_xlabel("Tiempo (s)")
    eje.set_ylabel("Temperatura (°C)")
    eje.grid(True, alpha=0.3)
    eje.legend()

    def actualizar(_):
        tiempo_actual = time.monotonic() - inicio
        try:
            temperatura = leer_temperatura_ohm(url)
        except (OSError, ValueError) as error:
            mensaje = f"No se pudo leer Open Hardware Monitor: {error}"
            estado.set_text(mensaje)
            if ultimo_error[0] != mensaje:
                print(mensaje)
                ultimo_error[0] = mensaje
        else:
            tiempos.append(tiempo_actual)
            temperaturas.append(temperatura)
            limite_tiempo = tiempo_actual - 60
            while tiempos[0] < limite_tiempo:
                tiempos.pop(0)
                temperaturas.pop(0)
            linea_datos.set_data(tiempos, temperaturas)
            eje.set_xlim(max(0, tiempo_actual - 60), max(10, tiempo_actual))
            print(f"\rCPU Package: {temperatura:.1f} °C", end="", flush=True)
            ultimo_error[0] = None

            try:
                parametros = calibrar_parametros_en_vivo(
                    tiempos, temperaturas, temp_ambiente
                )
            except (RuntimeError, ValueError, FloatingPointError) as error:
                linea_modelo.set_data([], [])
                estado.set_text(f"No se pudo calibrar el modelo: {error}")
            else:
                if parametros is None:
                    linea_modelo.set_data([], [])
                    estado.set_text(
                        f"CPU Package: {temperatura:.1f} °C | "
                        "esperando suficiente variación para calibrar"
                    )
                else:
                    calor_ajustado, k_ajustado = parametros
                    tiempos_modelo = np.linspace(tiempos[0], tiempos[-1], 100)
                    temperaturas_modelo = solucion_analitica(
                        tiempos_modelo - tiempos[0],
                        temperaturas[0],
                        temp_ambiente,
                        calor_ajustado,
                        k_ajustado,
                    )
                    linea_modelo.set_data(tiempos_modelo, temperaturas_modelo)
                    estado.set_text(
                        f"CPU Package: {temperatura:.1f} °C | "
                        f"k={k_ajustado:.5f} 1/s, Q={calor_ajustado:.3f} °C/s"
                    )

            eje.relim()
            eje.autoscale_view(scalex=False)

        return linea_datos, linea_modelo, estado

    animacion = FuncAnimation(figura, actualizar, interval=1000, blit=False, cache_frame_data=False)
    print(
        f"Leyendo CPU Package desde {url}; temperatura ambiente asumida: "
        f"{temp_ambiente:.1f} °C. Cierra la ventana para terminar."
    )
    plt.tight_layout()
    plt.show()
    plt.close(figura)


# --------------------------------------------
# 8) PROGRAMA PRINCIPAL
# --------------------------------------------
def principal():
    temp_ambiente = 25.0
    temp_inicial = 30.0

    if len(sys.argv) > 1 and sys.argv[1] == "--live":
        url_ohm = sys.argv[2] if len(sys.argv) > 2 else "http://localhost:8085/data.json"
        try:
            graficar_temperatura_en_vivo(url_ohm, temp_ambiente)
        except (OSError, ValueError, json.JSONDecodeError) as e:
            print(f"\nERROR al conectar con Open Hardware Monitor: {e}")
            print("Activa el servidor web de Open Hardware Monitor y verifica la dirección del archivo data.json.")
            sys.exit(1)
        return

    if len(sys.argv) > 1:
        ruta_csv = sys.argv[1]
        print(f"Intentando cargar datos desde: {ruta_csv}")
        try:
            tiempo_datos, temperatura_datos = cargar_datos_csv(ruta_csv)
            print("Datos cargados correctamente.")
        except Exception as e:
            print(f"\nERROR al leer el CSV: {e}")
            print("Asegúrate de que el archivo exista y tenga columnas: tiempo_s, temperatura_c")
            sys.exit(1)
    else:
        print("No se dio un CSV -> usando datos sintéticos de ejemplo.")
        tiempo_datos, temperatura_datos = generar_datos_sinteticos(temp_inicial, temp_ambiente)

    try:
        calor_ajustado, k_ajustado = calcular_constante_k(
            tiempo_datos, temperatura_datos, temp_inicial, temp_ambiente
        )
    except Exception as e:
        print(f"\nERROR al ajustar el modelo: {e}")
        sys.exit(1)

    print("\n=== CÁLCULO DE LA CONSTANTE K ===")
    print(f"k estimado  = {k_ajustado:.5f} 1/s")
    print(f"Q0 estimado = {calor_ajustado:.5f} °C/s")

    temperatura_simulada = simular_rk4(
        tiempo_datos, temp_inicial, temp_ambiente, k_ajustado, calor_ajustado
    )

    rmse, mae, r2 = calcular_metricas_error(temperatura_datos, temperatura_simulada)
    print("\n=== MÉTRICA DE ERROR (modelo vs. datos) ===")
    print(f"RMSE = {rmse:.4f} °C")
    print(f"MAE  = {mae:.4f} °C")
    print(f"R^2  = {r2:.4f}")

    graficar_resultados(tiempo_datos, temperatura_datos, tiempo_datos, temperatura_simulada, k_ajustado)


if __name__ == "__main__":
    principal()