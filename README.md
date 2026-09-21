# DINAMITA-PROYECTO-E.C
# DINAMITA 🧨
Proyecto 4.° semestre (Ing. Sistemas - Unilibre): modelado de temperatura en procesadores con la Ley de Enfriamiento de Newton y validación empírica [source: 2].

## ⚡ Flujo
1. **Captura:** CSV de Open Hardware Monitor (CPU Temp + Carga) [source: 2].
2. **Ingesta/Limpieza:** Pandas.
3. **Modelo:** EDO (`dT/dt = -k(T - Ta) + Q(t)`) resuelta con `scipy.integrate.solve_ivp` [source: 2].
4. **Validación:** Estimación de `k` por mínimos cuadrados (`curve_fit`) + RMSE.
5. **Plot:** Matplotlib [source: 2].

## 🛠️ Stack
Python | NumPy | SciPy | Pandas | Matplotlib | Open Hardware Monitor [source: 2]

## 👥 Equipo
* Joseph David Gómez Argote [source: 2]
* Bruce Adrián Abril Ramírez [source: 2]
