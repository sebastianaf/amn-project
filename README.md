# Proyecto Análisis y Métodos Numéricos

## Versión 1.0.0
Codigo base
- **Entrega 1:** No culminada, no es necesario revisarla aún
- **Entrega 2:** En ejecución, sólo se está implementando Richardson, Jocabi se hará luego de Richardson

## Versión 1.0.1
Se cambio la velocidad inicial, la cual estaba generando errores (V0 = 25 -> V0 = 0.3).

## Versión 1.0.2
Se ha disminuido el omega, para que la convergencia se amas rápida.
Las fronteras fueron puestas dentro del bucle que calcula las iteraciones de Richardson

## Versión 1.0.3
Los omega se dividiero en dos, hay un omega para cada matriz. Esto se hace para acelerar la convergencia de y, ya que suma valores extremadamente pequeños.

## Versión 1.0.4
Es esta implementando la viga 1. Está terminada la pared izquierda y la superior.

## Versión 1.1.0
Implementación con viga 1 completa

## Versión 1.1.1
Las vigas tienen mejor comportamiento

## Versión 1.2.0
Se corrigieron errores en la forntera superior. La matriz de magnitudes fué corregida, ya que pitágoras no aplica en este caso.
Para solucionar esto, se sumaron los tamaños de cada vector en u y w, y se saco su promedio.

## Versión 1.2.1
Se ha implementado la viga dos
Se solucionaron algunos errores
- Está en capacidad de ejecutar más iteraciones (antes sólo hacia 30)
- Ahora es mas rápido en el cálculo

Aún hay algunos errores. El más relevante es que no funciona la pared izquierda de la viga 2

## Versión 1.2.2
Ya están en funcionamiento todas las vigas (Aún hay errores)
Se añadieron archivos en blanco para las demás entregas
Hay un prototipo de la adaptación del código del libro guía

## Versión 1.2.3
Se corrigió la adaptación del código del libro

## Versión 1.2.4
Código final Richardosn y jacobi

## Versión 1.2.5
Corrección en condición de rontera de la entrega 2.

## Versión 1.2.6
Actualizacipon código libro. Prototipo de NewtonRaphson