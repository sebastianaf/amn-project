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