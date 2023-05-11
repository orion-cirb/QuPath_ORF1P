import qupath.ext.biop.cellpose.Cellpose2D
import qupath.lib.gui.dialogs.Dialogs
import qupath.lib.scripting.QP

// Customize how the Cellpose detection should be applied
cellpose = Cellpose2D.builder('cyto2')
            .preprocess(ImageOps.Filters.median(1))   // List of preprocessing ImageOps to run on the images before exporting them
            .normalizePercentiles(1, 99)              // Percentile normalization
            .pixelSize(0.5)                           // Resolution for detection
            .channels('Cy5')                          // Select detection channel(s)
//          .maskThreshold(-0.2)                      // Threshold for the mask detection, defaults to 0.0
//          .flowThreshold(0.5)                       // Threshold for the flows, defaults to 0.4
            .diameter(30)                             // Median object diameter. Set to 0.0 for the `bact_omni` model or for automatic computation
            .setOverlap(60)                           // Overlap between tiles (in pixels) that the QuPath Cellpose Extension will extract. Defaults to 2x the diameter or 60 px if the diameter is set to 0
            .simplify(1.5)
            .constrainToParent(false)
            .measureShape()                           // Add shape measurements
            .measureIntensity()                       // Add intensity measurements (in all compartments)
//          .doLog()
            .build()
	
// Define which objects will be used as the 'parents' for detection
// Use QP.getAnnotationObjects() if you want to use all annotations, rather than selected objects
def pathObjects = QP.getSelectedObjects()

// Run detection for the selected objects
def imageData = QP.getCurrentImageData()
if (pathObjects.isEmpty()) {
    QP.getLogger().error("No parent objects are selected!")
    return
}
cellpose.detectObjects(imageData, pathObjects)
println('Done!')