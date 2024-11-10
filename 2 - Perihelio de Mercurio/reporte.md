# Precesión del Perihelio de Mercurio con Métodos de Runge-Kutta
## Resumen

Este estudio investiga la precesión del perihelio de Mercurio utilizando métodos numéricos, específicamente el método de Runge-Kutta de segundo orden (RK2). Se desarrolla un modelo computacional que incorpora tanto la fuerza gravitacional del Sol como los efectos perturbadores de otros planetas, representados por un parámetro α. El trabajo comienza con una introducción a los métodos de Runge-Kutta y una descripción detallada de la órbita de Mercurio, incluyendo las ecuaciones que gobiernan su movimiento.

La metodología empleada implica la aplicación del método RK2 a un sistema de ecuaciones diferenciales de primer orden derivadas de las ecuaciones de movimiento de Mercurio. Se implementa una función en Python, `simular_orbita_mercurio`, que simula la órbita del planeta con diversos parámetros iniciales. Esta simulación permite observar y medir la precesión del perihelio bajo diferentes condiciones.

El estudio analiza cómo el parámetro α, que representa las correcciones relativistas y los efectos de otros cuerpos celestes, influye en la precesión del perihelio. Se desarrolla una función adicional, `calcular_precesion`, para cuantificar esta precesión basándose en los resultados de la simulación.

Este enfoque computacional proporciona una herramienta valiosa para estudiar y visualizar el fenómeno de la precesión del perihelio de Mercurio, ofreciendo insights sobre los factores que influyen en este importante efecto astronómico y relativista.

## Introducción

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

## Metodología
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

Implementamos la función `simular_orbita_mercurio` en Python, la cual simula la órbita de Mercurio utilizando el método de Runge-Kutta de segundo orden (RK2). Esta función toma los siguientes parámetros:

- `x_inicial` (float, opcional): Posición inicial en el eje x. Si no se proporciona, se calcula como (1+e)*a, donde e es la excentricidad y a es el semieje mayor de la órbita.
- `y_inicial` (float, opcional): Posición inicial en el eje y. Por defecto es 0.
- `vx_inicial` (float, opcional): Velocidad inicial en el eje x. Por defecto es 0.
- `vy_inicial` (float, opcional): Velocidad inicial en el eje y. Si no se proporciona, se calcula utilizando la ecuación de la velocidad orbital.
- `iteraciones` (int, opcional): Número de iteraciones para la simulación. Por defecto es 20000.
- `delta_t` (float, opcional): Tamaño del paso de tiempo para la simulación. Por defecto es 0.0001.
- `alpha` (float, opcional): Parámetro que representa la corrección relativista. Por defecto es 0.0008.

La función inicializa las variables necesarias, incluyendo las constantes como GM_s (producto de la constante gravitacional y la masa del Sol), el semieje mayor (a), y la excentricidad (e) de la órbita de Mercurio.

Luego, realiza la simulación utilizando el método RK2, calculando en cada iteración:
1. Las variables radiales (r y v).
2. Los valores intermedios (k1) para posición y velocidad.
3. Los valores finales (k2) para posición y velocidad.
4. Actualiza los valores de posición, velocidad y tiempo para el siguiente paso.

La función devuelve tres listas:
1. `tiempos`: Una lista con los tiempos de cada paso de la simulación.
2. `posiciones`: Una lista de tuplas (x, y) representando las posiciones en cada paso.
3. `velocidades`: Una lista de tuplas (vx, vy) representando las velocidades en cada paso.

Esta implementación nos permite simular la órbita de Mercurio con diferentes parámetros iniciales y observar cómo estos afectan la precesión del perihelio.
### Encontrando la precesión para un α dado

Para analizar la precesión del perihelio de Mercurio, graficamos la relación entre el parámetro alfa y el ángulo de la posición con respecto al eje x en el perihelio. Esto nos permite visualizar cómo el parámetro alfa afecta a la precesión de la órbita.

Luego implementamos una función `calcular_precesion` para determinar la precesión del perihelio basados en los resultados de una simulación. Esta función realiza los siguientes pasos:

1. Identifica los puntos de perihelio en la órbita simulada. Estos son los puntos donde la distancia al Sol es mínima en cada revolución. Utilizamos dos criterios para identificar las mínimos:
- Son aquellos puntos que son menores que sus dos vecinos
- Son aquellos puntos donde la razón de cambio de la distancia radial se desvanece.
\frac{dr}{dt} = \frac{xv^x + yv^y}{r} = 0

2. Para cada perihelio se calcula el ángulo de la posición con respecto al eje x. La razón de cambio de este ángulo con respecto al tiempo es la precesión buscada.

### Extrapolando el valor de la Precesión para un Valor de \alpha Difícil de simular

Una vez que tenemos la función calcular_precesion podemos utilizarla para simular la órbita de Mercurio para diferentes valores de \alpha, y luego encontrar la pendiente de la recta que relaciona el parámetro alfa con la precesión, para luego extrapolar el valor el valor de la precesión cuando \alpha = 1.1E-8, el cual es difícil de simular dado que tomaría mucho tiempo observar la precesión.

## resultados
![Órbita simulada de Mercurio](./figures/orbita.png)

La figura anterior muestra la órbita simulada de Mercurio utilizando nuestro método de Runge-Kutta de segundo orden. Como se puede observar, la órbita presenta una forma elíptica característica, con el Sol ubicado en uno de los focos de la elipse.

![Ángulo vs. Tiempo para la órbita de Mercurio](./figures/angulo_vs_tiempo.png)

En la figura anterior podemos observar la relación lineal entre el ángulo y el tiempo, lo que nos permite reconocer la precesión del ángulo como la razón de cambio del ángulo con respecto al tiempo, es decir la pendiente de la recta trazada por el ángulo.

Finalmente, la tabla1 muestra la precesión calculada para los valores de alfa requeridos (0.0008, 0.001, 0.002, 0.004) basados en los cuales podemos encontrar la pendiente de 10848.00 para la recta que relaciona alfa con la precesión. Por lo que para \alpha = 1.1E-8 observado tendríamos una precesión de 42.96 segundos de arco por siglo.

## Discusión
Notemos que este valor de 42.96 segundos por siglo es lo que se obtiene luego de convertir la precesión devuelta por nuestro código en grados por año. Los detalles se pueden ver en el notebook, donde la conversión es explícita.

El cálculo de la pendiente que relaciona la precesión con el parámetro alfa se hizo tomando en cuenta los dos primeros puntos de la lista, dado que es una relación lineal, utilizar métodos más avanzados no nos mostraba ninguna mejora en el valor de la precesión, peor aún, dado que al aumentar alfa, la órbita se aleja cada vez más de los valores observados, los puntos con \alpha = 0.002, 0.004 tienden a exagerar más la precesión, como se puede observar al comparar la línea de tendencia lineal en la figura3 con los valores de la simulación.


## Conclusiones
Fuimos capaces de implementar la función simular_orbita_mercurio y verificar el correcto funcionamiento de la simulación. Utilizando esta función pudimos encontrar la relación entre el parámetro alfa de correción relativista y la precesión de la órbita. Verificando así nuestra implementación del algoritmo RK2.

## bibliografía 