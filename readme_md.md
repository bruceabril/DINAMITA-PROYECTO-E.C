# Dinamita: Modelado Dinámico y Simulación Térmica de un Procesador Mediante Ecuaciones Diferenciales Ordinarias

Proyecto académico desarrollado para la asignatura de **Ecuaciones Diferenciales** (Ingeniería en Sistemas - Universidad Libre). Este repositorio contiene la implementación matemática, analítica y de software para modelar, simular y monitorear el comportamiento térmico de un microprocesador frente a variaciones de carga de trabajo.

##  Modelo Matemático

El comportamiento térmico del procesador se fundamenta en la **Ley de Enfriamiento de Newton** y el balance de energía del sistema, representándose mediante una Ecuación Diferencial Ordinaria (EDO) lineal de primer orden:

$$
\frac{dT}{dt} = -k(T - T_a) + Q(t)
$$

Donde:
* $T$: Temperatura instantánea del procesador.
* $T_a$: Temperatura ambiente.
* $k$: Constante de disipación térmica.
* $Q(t)$: Potencia térmica generada por el procesador.

Bajo la suposición de una potencia constante $Q$, la solución analítica implementada en el sistema es:

$$
T(t) = T_a + \frac{Q}{k} + \left(T_0 - T_a - \frac{Q}{k}\right)e^{-kt}
$$

##  Instrucciones de Ejecución

Si alguien descarga el proyecto, puede ejecutarlo así:

1. Descarga el repositorio como ZIP desde GitHub (`Code` $\rightarrow$ `Download ZIP`) y descomprímelo.
2. Instala Python. Luego abre PowerShell dentro de la carpeta descomprimida y ejecuta:

```powershell
python -m pip install numpy pandas scipy matplotlib requests
```

3. Ejecuta el programa normal:

```powershell
python ".\Proyecto.py.py"
```

Eso usa datos sintéticos y genera `comparacion_modelo_real.png`. El código está en `Proyecto.py.py`.

##  Monitoreo en Vivo (Tiempo Real)

Para la gráfica en vivo, también debe instalar y abrir **Open Hardware Monitor**, activar su servidor web y luego ejecutar:

```powershell
python ".\Proyecto.py.py" --live
```

Por defecto busca los datos en `http://localhost:8085/data.json`. Si el servidor está en otra dirección, se puede indicar después de `--live`. La lectura en vivo requiere que Open Hardware Monitor exponga el sensor **CPU Package**; si no está disponible, esa parte no funcionará en ese equipo.

##  Autores

* **Joseph David Gomez Argote**
* **Bruce Adrián Abril Ramírez**
* **Luis Hernando Rico Batanero**

**Docente Asesor:** Jorge Enrique Taboada Álvarez  
**Institución:** Universidad Libre