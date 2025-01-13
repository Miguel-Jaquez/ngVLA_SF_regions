El folder "utils" contiene herramientas necesarias para correr tanto los modelos de transferencia radiativa con RADMC-3D, como las observaciones sinteticas del ngLVA con CASA. 

Breve descripcion de los archivos del folder utils:

- Hydrogen_recom_lines_in_use.csv: Este archivo contiene las propiedades basicas de las lineas de recombinacion que caen dentro de  la banda de 93 GHz del ngVLA. Las columnas del archivo son:
	a) Line : nombre de la linea de recombinacion.
	b) nsup: nivel de energia antes de la transicion
	c) ninf: nivel de energia despues de la transicion
	d) Frequency_splatalogue : Frequencia del foton resultante de la transicion del nivel "nsup" a "ninf". Esta frequencia fue obtenida del splatalogue.
	e) Frequency_radmc : Frequencia del foton resultante de la transicion del nivel "nsup" a "ninf". Esta frequencia es con la que trabaja internamente RADMC-3D.
	f) Desfase[GHz]: Diferencia entre la frequencia de RADMC y el splatalogue.

- ngVLA_sensitive_calculator.py: Script de Viviana Rosero que calcula las tablas de performance del ngVLA. Este script usa los archivos "receiver_data.pkl" y "subarray_data.pkl". 

- plot_helpers.py y myplt_style_andizq.py: Script de Andres Izquierdo para graficar. 

- simmos: folder donde vienen los archivos ".cfg" con las configuraciones de antenas del ngVLA. Si se actualizan los archivos d eCASA podemos prescindir de este folder.

- utils_run_radmc.py: En este script se define la funcion que corre RADMC-3D y genera las imagenes fits resultantes de la transferencia radiativa que son (usando como ejemplo la linea de H38):
	a) img_rl_jet_H38.fits: modelo output de RADMC.
	b) img_rl_jet_line_H38.fits : solo la emision de linea del archivo a)
	c) img_rl_jet_cont_H38.fits : solo emision de continuo del archivo a)
	d) img_rl_jet_line_convl_H38.fits : modelo "idealizado" de observacion sintetica del modelo b). Se usa una PSF con un beam gaussiano. Se agrega el ruido usando el script de Rosero.
	e) img_rl_jet_cont_convl_H38.fits: modelo "idealizado" de observacion sintetica del modelo c). Se usa una PSF con un beam gaussiano.   Se agrega el ruido usando el script de Rosero.

 
