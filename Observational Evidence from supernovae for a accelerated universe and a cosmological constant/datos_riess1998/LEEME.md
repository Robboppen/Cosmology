# Datos de Riess et al. (1998), AJ 116, 1009

Tablas extraídas de la versión publicada del paper (AJ). Formato CSV; los signos negativos están corregidos y los campos vacíos indican "..." en el original.

| Archivo | Tabla | Contenido |
|---|---|---|
| tabla1_espectroscopia.csv | 1 | Espectros de identificación: telescopio, rango espectral (nm), z y su origen, SN de comparación de la Fig. 1 (edad en días respecto del máximo en B) |
| tabla2_fotometria_observada.csv | 2 | Fotometría observada en los filtros B45, V45, B35, V35 (81 épocas). Columna *_filtro_alt: R, I o V si se observó en otra banda |
| tabla3_estrellas_estandar_locales.csv | 3 | Magnitudes de las estrellas estándar locales de cada campo (Fig. 2). Columna *_nota: V o R si la magnitud está en esa banda |
| tabla4_curvas_luz_reposo.csv | 4 | Curvas de luz en B y V **en el sistema en reposo**, con las correcciones K aplicadas (K_B, K_V). Base de la Fig. 3 |
| tabla5_MLCS_lejanas.csv | 5 | Ajuste MLCS de las lejanas: m_max, Delta, A_B, mu0 |
| tabla6_dm15_lejanas.csv | 6 | Ajuste de plantillas Δm15(B) de las lejanas |
| tabla7_snapshot_lejanas.csv | 7 | Distancias "snapshot" de las lejanas con pocos datos |
| tabla9_snapshot_cercanas.csv | 9 | Distancias snapshot de las cercanas |
| tabla10_cercanas_MLCS_dm15.csv | 10 | 27 SNe cercanas: log cz, MLCS y Δm15, tipo de galaxia (E/L) |

Notas
- JD está dado como JD − 2 450 000.
- z de las cercanas: z = 10^(log_cz) / 299792.458.
- La fecha del máximo (t_max), necesaria para el eje "edad" de la Fig. 3, no está tabulada: sale de ajustar las curvas de la Tabla 4. Edad en reposo = (JD − JD_max)/(1 + z).
- Inconsistencia del paper: la Tabla 1 da z = 0.24 (1995ao) y 0.30 (1995ap); la Tabla 7 da 0.30 y 0.23.
- Para el χ²: σ_v = 200 km/s, y 2500 km/s adicionales en cuadratura cuando el z proviene de los rasgos de la SN (columna origen_z de la Tabla 1).
