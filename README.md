# py-microwave
A microwave toolbox written in Python3

The repository includes three modules:

* **Microwave Engineering Tools** mwave.py
* **Smith Chart** Creation in Python smith.py
* **Filter Synthesis** in Python filter.py

### mwave.py

A module including a lot of methods for:

 * creating and cascading two-port networks
 * convert two-port networks between S, Z, Y, ABCD representation
 * plot frequency response
 * import/export touchstone and mdif
 * calculating microstrip line geometries, vias, coupled lines 
 * convert and handle complex values
 * design amplifiers with stability, gain calculations, gain circles, stability circles, noise behaviour

### smith.py

  create and manipulate nice looking smith charts 
  
  ![Smith Chart](https://github.com/sfpeik/py-microwave/blob/main/examples/schematic.svg "Smith")
  
  ![Smith Chart](https://github.com/sfpeik/py-microwave/blob/main/examples/smithchart.svg "Smith")
  
  ## filter.py
  
  Tools for filter synthesis
  
   * g-parameters
   * elliptic filter prototypes
   * Coupling matrix creation
   * Coupling matrix analysis, response plots
   * Coupling value visualisaton
   * Coupling Matrix extraction


   
  
  
