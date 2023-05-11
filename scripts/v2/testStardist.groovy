import qupath.ext.stardist.StarDist2D
import qupath.lib.gui.dialogs.Dialogs
import qupath.lib.scripting.QP

def pathProject = buildFilePath(PROJECT_BASE_DIR)
def modelPath = buildFilePath(pathProject, 'models', 'dsb2018_heavy_augment.pb')

// Customize how the StarDist detection should be applied
def stardist = StarDist2D
    .builder(modelPath)
    .channels('DAPI')   
    .preprocess(ImageOps.Filters.median(2))
    .normalizePercentiles(1, 99)       // Percentile normalization
    .pixelSize(0.5)                 // Select detection channel
    .threshold(0.75)          // Prediction threshold
    .simplify(0)
    .constrainToParent(false)
    .measureShape()                    // Add shape measurements
    .measureIntensity()                // Add intensity measurements
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
stardist.detectObjects(imageData, pathObjects)
stardist.close() // This can help clean up & regain memory
println('Done!')