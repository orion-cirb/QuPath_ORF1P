import org.apache.commons.io.FilenameUtils
import static qupath.lib.scripting.QP.*

// Get current image name
def imageDir = new File(project.getImageList()[0].getURIs()[0]).getParent()
def imageName = getCurrentImageData().getServer().getMetadata().getName()
def imgNameWithOutExt = FilenameUtils.removeExtension(imageName)

// Clear all detections
clearDetections()

// Get all annotations
def annotations = getAnnotationObjects()

// Find annotations files for current image
def p = ~/${imgNameWithOutExt}.*\.annot/
def resultsDir = new File(buildFilePath(imageDir+'/Results'))
resultsDir.eachFileMatch(p) {file ->
    new File(file.path).withObjectInputStream {
        def detections = it.readObject()
        def parentName = detections.getAt(0).getName()
        println 'Adding detections in ' + parentName + ' region'
        def parent = annotations.find{it.getName() == parentName}
        parent.addChildObjects(detections.getAt(0).getChildObjects())
    }
}
fireHierarchyUpdate()
println 'Enjoy!'

