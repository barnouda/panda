# PANDA (Average Density from Nodes Along Path, in reverse order)

Computes the contribution of grid node densities to the average density along muon paths.

## References

When using this code, please cite:
* Barnoud A., Cayol V., Niess V., Cârloganu C., Lelièvre P., Labazuy P., Le Ménédeu E., Bayesian joint muographic and gravimetric inversion applied to volcanoes, Geophysical Journal International, Volume 218, Issue 3, September 2019, Pages 2179–2194, https://doi.org/10.1093/gji/ggz300

## PANDA VS PANDA2

panda (ancien code)  
- contribution of densities along a single line of sight per bin  
  
panda2 (nouveau code)  
- analytical integration instead of numerical integration along lines of sight  
- multiple sub lines of sight per bin (ie per data)  
- output files with averaged densities for a given model and rock depths  
- TO DO weighting par bin en fonction de la distance au bin principal  
- threading openmp  
- [TO DO possibility of densities defined by nodes or constant by blocs]  
- TO DO possibility of skiping 'none' line in parameter file (count nb of lines first)  
- WARNING: code en dur: pas de subdivision pour les bins en azimuth et elevation (eg 0.2 deg)

## Inputs

pdd				| prefix for output files
------------------------------- | ----------------------------------------------------------  
composite_DEM_PDD_50cm.bin      | topography file name  
pdd_dxyz_25_rho_1000.par        | grid node parameters ascii file  
pdd_dxyz_25_rho_1000.bin	| grid nodes densities binary file  
geometry_muog_seff15_rate70.txt | muographic data geometry (x0,y0,z0,azim,elev,dazim,delev)  
muog_bin_weights.txt            | weighted linear combinations between bins  
  
NB: for panda, data geometry = (x0, y0, z0, azim, elev)  

## Outputs

- ascii file with averaged rock thickness and density per bin  
- bin file with sensitivity matrix with contribution of all nodes to the averaged densities  

## Warning

Maximum of 9 telescopes => not true anymore because geometry file modifies  
=> now, the telescope position is repeated for each line of sight, so that the number of telescopes is transparent to the code (this is why the telescope caracteristics are repeated as many times as data in the data file...)  

## Algorithm (panda)  

rho moyen = somme ds * rhos / somme ds  
rhos = somme csi * rhoi  

rho moyen = somme ds * (somme csi * rhoi) / somme ds

nl rho moyen - index = il - rl rho
nn density nodes - index = in - rn rho
ns integration samples along path

rl = somme ds(is) * rs(is) / somme ds(is)
rs(is) = somme cs(in) * rn(in)
=> rl = somme ds(is) * (somme cs(in) * rn(in)) / somme ds(is)

M(il,in) contribution of nodes to rho along path
M(il,in) = somme_s ds(is) * (somme_n cs(in)) / somme_s ds(is)

* boucle sur les nl paths - on rempli uniquement M(il,?)
  * boucle sur les ns samples along path
    * boucle sur les nn nodes of grid (dans trilinint) - on rempli uniquement M(il,in)
      C(is,in) => somme * somme ds => M(il,in)

## Alogirthm panda2  

- read input parameters and files  
  - readparam  
  - readgridparam  
  - readdensitygrid  
  - numberofmuogdata  
  - readmuogdata2  
  - readtopoparam  
  - readtopo  
- determine all blocs delimited by 8 nodes (blocs of non uniform density...)  
  - construct and index blocs  
    - loop x = xmin:dx:xmax-dx  
      - loop y  
        - loop z  
          - bloc limits  
          - indices of blocs summits    
- determine maximum length for lines of sight (whatever the telescope is)  
  - distance between the telescopes and the 8 cornes of the full mesh  
  - keep largest distance = smax  
- initialize output variables  
- loop over data (ie bins)  
  - determine optimal angular step theta for lines of sight  
    - period sampling = period nyquist / 2 = dxyz / 2  
    - max length = smax  
    - tan(thetatmp/2) = dxyz / 2 / smax  
    - thetatmp = 2 * arctan(dxyz/2/smax)  
    - nLoS per bin = ceil((azim step = elev step) / thetatmp)  
    - theta optimal = (azim step = elev step) / nLoS  
    - define sub lines of sight inside bin  
  - double loop over (azim, elev) of lines of sight  
    - from degrees to radians  
    - direction (a,b,c) defining the line of sight  
    - intersections with the 6 faces of the full mesh  
       - if 0 intersection => ignore (LoS outside model)  
       - if 1 intersection => ignore (LoS tangantial to model)  
       - if 2 intersections => in and out points for the full mesh  
    - identify topography crossings  
      - coordinates of the line of sight, step along line of sight = topography resolution step  
      - bilinint: topography at the location of line of sight  
      - find points entering and exiting topography  
      - interpolate between two consecutive points inside and outside topo  
    - loop over segments below topography  
      - loop over blocs defined at the begining  
        - if segment does not cross the bloc => ignore bloc, go to next bloc 
        - interscube: intersection pts with bloc (if no intersection => ignore bloc)  
        - compute length traveled inside bloc  
        - trilinintcoeffs2: compute the contribution of the 8 nodes to the INTEGRATED density along the line of sight inside this bloc  
    - update output variables (rock depth, sensitivity matrix)  
      - update values  
      - weight by cos(elevation) 
      - divied by traveled length (from integrated density to averaged density)  
  - divied by sum of cos(elevation)  
- write output files  
  - sensitivity matrix (binary file)  
  - rock depth and averaged densities (ascii file)  

## Example (./example)

The example dataset is a synthetic volcano (conic topography) with laterally varying densities and two muon telescopes.

The synthetic topography, model, data and input files for modelling are generated with the example_toy_cone.m matlab code and its dependencies. You may ignore the *gravi* and *mag* files for use with another code.

To run the code with the example dataset:
* build code with the makefile
* run the code in the example folder, ex: ../bin/panda input_muog_modelling.in

