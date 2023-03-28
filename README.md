# QuPath_ORF1P

* **Developed for:** Tom
* **Team:** Proschiantz
* **Date:** March 2022
* **Software:** QuPath

### Images description

2D images of mouse brain sections taken with the Axioscan

4 channels: 
  1. *DAPI:* nuclei
  2. *EGFP:* microglia (?)
  3. *Cy3:* ORF1p cells
  4. *Cy5:* NeuN cells

### Plugin description

Version 1:
* Detect nuclei, microglia and ORF1p and NeuN cells with Stardist
* Compute their colocalization

Version 2:
* Ignore EGFP channel
* Detect nuclei with Stardist and ORF1p and NeuN cells with Cellpose
* Compute their colocalization
* Compute background noise automatically in each annotation by taking into account every pixel that is not belonging to a cell (after Stardist/Cellpose detection)
* Handle hierarchy by skipping the analysis of parents regions

### Dependencies

* **Stardist** QuPath extension +  *dsb2018_heavy_augment.pb* model
* **Cellpose** QuPath extension + *Cellpose** conda environment + *cyto2* model

### Version history

Version 1 (*getnuclei.groovy*) released on September 1, 2022.

Version 2 (*getnuclei_improved.groovy*) released on March 28, 2022.
