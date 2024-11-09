# Precesión del Perihelio de Mercurio con Métodos de Runge-Kutta
## resumen
## introducción
### Métodos de Runge-Kutta

Los métodos de Runge-Kutta son una familia de algoritmos numéricos iterativos utilizados para aproximar soluciones de ecuaciones diferenciales ordinarias (EDOs). Estos métodos son ampliamente utilizados en física computacional y otras disciplinas científicas.

Los métodos de Runge-Kutta funcionan calculando incrementos intermedios entre cada paso de integración, lo que permite una mejor aproximación de la solución real. Para este estudio, utilizaremos el método de Runge-Kutta de segundo orden (RK2), también conocido como el método del punto medio. Este método ofrece un buen equilibrio entre precisión y eficiencia computacional para nuestro problema.

La fórmula general para el método RK2 es:

y_{n+1} = y_n + k_2

Donde:
- y_n es el valor actual de la función
- h es el tamaño del paso
- k_1 = h * f(t_n, y_n)
- k_2 = h * f(t_n + h/2, y_n + k_1/2)

En el caso de la órbita de Mercurio, el método RK2 nos permite resolver numéricamente las ecuaciones de movimiento derivadas de la fuerza gravitacional.

La aplicación del método RK2 a este problema implica descomponer las ecuaciones de movimiento en un sistema de ecuaciones diferenciales de primer orden y luego aplicar el algoritmo iterativamente para avanzar la posición y velocidad de Mercurio en el tiempo. Aunque el RK2 es menos preciso que métodos de orden superior como el RK4, es suficiente para nuestro estudio y nos permite un cálculo más eficiente, lo cual es beneficioso para simulaciones de largo plazo como la que necesitamos para observar la precesión del perihelio de Mercurio.

### La órbita de Mercurio

Para modelar con precisión el movimiento de Mercurio, es necesario considerar no solo la fuerza gravitatoria del Sol, sino también los efectos de otros planetas.

La ecuación para la fuerza gravitatoria sobre Mercurio, que incluye el efecto de otros planetas como una fuerza externa, se puede expresar de la siguiente manera:

F = G * M_s * M_m * (1/r^2) * (1 + α / r^2)

Donde:
- F es la fuerza gravitatoria
- G es la constante de gravitación universal
- M_s es la masa del Sol
- M_m es la masa de Mercurio
- r es la distancia entre el centro del Sol y el centro de Mercurio
- α es un parámetro que representa las correcciones debidas a los efectos de otros planetas.

## antecedentes
## metodología
### Aplicando los métodos de Runge-Kutta a la órbita de Mercurio

Para aplicar el método de Runge-Kutta de segundo orden (RK2) a la órbita de Mercurio, necesitamos primero expresar las ecuaciones de movimiento en un sistema de ecuaciones diferenciales de primer orden.

Partiendo de la fuerza gravitatoria que actúa sobre Mercurio:

F = G * M_s * M_m * (1/r^2) * (1 + α / r^2)

Podemos expresar las ecuaciones de movimiento en coordenadas polares (r, θ) de la siguiente manera:

1. d²r/dt² - r * (dθ/dt)² = -G * M_s * (1/r²) * (1 + α / r²)
2. (1/r) * d/dt(r² * dθ/dt) = 0

Para aplicar el método RK2, necesitamos convertir estas ecuaciones de segundo orden en un sistema de ecuaciones de primer orden. Definimos las siguientes variables:

- u = r
- v = dr/dt
- w = dθ/dt

Ahora, podemos reescribir nuestras ecuaciones como un sistema de primer orden:

1. du/dt = v
2. dv/dt = u * w² - G * M_s * (1/u²) * (1 + α / u²)
3. dw/dt = -2 * v * w / u

Este sistema de ecuaciones es el que utilizaremos con el método RK2. Para cada paso de tiempo, aplicaremos el método a estas tres ecuaciones simultáneamente.

El algoritmo RK2 para este sistema se puede expresar de la siguiente manera:

1. Calcular los valores intermedios:
   k1_u = δt * v_n
   k1_v = δt * (u_n * w_n² - G * M_s * (1/u_n²) * (1 + α / u_n²))
   k1_w = δt * (-2 * v_n * w_n / u_n)

2. Calcular los valores finales:
   k2_u = δt * (v_n + k1_v/2)
   k2_v = δt * ((u_n + k1_u/2) * (w_n + k1_w/2)² - G * M_s * (1/(u_n + k1_u/2)²) * (1 + α / (u_n + k1_u/2)²))
   k2_w = δt * (-2 * (v_n + k1_v/2) * (w_n + k1_w/2) / (u_n + k1_u/2))

3. Actualizar los valores para el siguiente paso:
   u_{n+1} = u_n + k2_u
   v_{n+1} = v_n + k2_v
   w_{n+1} = w_n + k2_w

Donde δt es el tamaño del paso de tiempo.

Este algoritmo se aplicará iterativamente para calcular la trayectoria de Mercurio a lo largo del tiempo, permitiéndonos observar y medir la precesión del perihelio.

### Implementando una función de iteración

Para implementar el método RK2 que hemos derivado, podemos crear una función en Python que realice las iteraciones. A continuación, se presenta un pseudocódigo que ilustra cómo podría ser esta función:

```python
import numpy as np

def simular_orbita_mercurio(
        x_inicial
        , y_inicial
        , vx_inicial
        , vy_inicial
        , tiempo_total
        , delta_t
        , alpha = 0.0008
        ):
    # Inicializar variables
    x = x_inicial
    y = y_inicial
    vx = vx_inicial
    vy = vy_inicial
    t = 0  # Tiempo inicial

    # Constantes
    G = 6.67430e-11  # Constante gravitacional (m^3 kg^-1 s^-2)
    M_s = 1.989e30   # Masa del Sol (kg)

    # Listas para almacenar resultados
    tiempos = [t]
    posiciones = [(x, y)]
    velocidades = [(vx, vy)]

    # Bucle principal de la simulación
    while t < tiempo_total:
        # Calcular variables radiales
        r = np.sqrt(x**2 + y**2)
        v = np.sqrt(vx**2 + vy**2)
        
        # Calcular los valores intermedios (k1)
        k1_x = delta_t * vx
        k1_y = delta_t * vy
        k1_vx = delta_t * (-G * M_s * x / r**3) * (1 + alpha / r**2)
        k1_vy = delta_t * (-G * M_s * y / r**3) * (1 + alpha / r**2)

        # Calcular los valores finales (k2)
        x_mid = x + k1_x/2
        y_mid = y + k1_y/2
        vx_mid = vx + k1_vx/2
        vy_mid = vy + k1_vy/2
        r_mid = np.sqrt(x_mid**2 + y_mid**2)

        k2_x = delta_t * vx_mid
        k2_y = delta_t * vy_mid
        k2_vx = delta_t * (-G * M_s * x_mid / r_mid**3) * (1 + alpha / r_mid**2)
        k2_vy = delta_t * (-G * M_s * y_mid / r_mid**3) * (1 + alpha / r_mid**2)

        # Actualizar los valores para el siguiente paso
        x += k2_x
        y += k2_y
        vx += k2_vx
        vy += k2_vy
        t += delta_t

        # Guardar los resultados en cada paso
        tiempos.append(t)
        posiciones.append((x, y))
        velocidades.append((vx, vy))

    return tiempos, posiciones, velocidades
```

### Encontrando la precesión para un α dado

Para determinar la precesión del perihelio de Mercurio a partir de los resultados de nuestra simulación, comenzamos por identificar los puntos de perihelio en la órbita. Estos ocurren cuando Mercurio está más cerca del Sol, lo que corresponde a mínimos locales en la distancia radial r = √(x² + y²). En la práctica, buscamos puntos donde r(t) sea menor que sus valores adyacentes.

Una vez identificados los perihelios, calculamos el ángulo entre dos puntos de perihelio consecutivos utilizando el producto escalar. Si (x₁, y₁) y (x₂, y₂) son dos perihelios consecutivos, el ángulo θ entre ellos se obtiene de cos(θ) = (x₁x₂ + y₁y₂) / (r₁r₂), donde r₁ y r₂ son las magnitudes de los vectores correspondientes.

La precesión por órbita, Δθ, es la diferencia entre este ángulo y 2π radianes (una vuelta completa). Para expresar la precesión en unidades más convencionales, la convertimos a segundos de arco por siglo. Esto implica multiplicar Δθ por factores de conversión apropiados y por el número de órbitas en un siglo.

El valor de α en nuestra simulación afecta directamente esta precesión. Ajustando α y comparando la precesión resultante con el valor observado de aproximadamente 43 segundos de arco por siglo, podemos calibrar nuestro modelo para que coincida con las observaciones reales del perihelio de Mercurio. Este proceso nos permite validar nuestra simulación y explorar cómo diferentes valores de α afectan la órbita de Mercurio.


## resultados
## discusión
## conclusiones
Fuimos capaces de calcular el perihelio de Mercurio para 
## bibliografía 