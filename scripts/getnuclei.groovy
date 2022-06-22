import org.apache.commons.io.FilenameUtils
import qupath.lib.objects.classes.PathClassFactory
import static qupath.lib.gui.scripting.QPEx.*
import qupath.ext.stardist.StarDist2D
import qupath.lib.gui.dialogs.Dialogs
import qupath.lib.objects.*
import qupath.lib.geom.Point2

// init project
setImageType('Fluorescence')
def project = getProject()
def pathProject = buildFilePath(PROJECT_BASE_DIR)
def pathModel = buildFilePath(pathProject,'models','dsb2018_heavy_augment.pb')
if (pathModel == null)
    print('No model found')
def imageDir = new File(project.getImageList()[0].getUris()[0]).getParent()

// create results file and write headers
def resultsDir = buildFilePath(imageDir, 'Results')
if (!fileExists(resultsDir)) {
    mkdirs(resultsDir)
}
def resultsFile = new File(buildFilePath(resultsDir, 'Results.csv'))
resultsFile.createNewFile()
def resHeaders = 'Image Name\tAnnotation Name\tArea\tnb DAPI\tDAPI Mean intensity\tnb EGFP\tEGFP Mean intensity\tnb Cy3\tCy3 Mean intensity\tnb Cy5\tCy5 Mean intensity\tnb Cy5Cy3\n'
resultsFile.write(resHeaders)

// Classpath definition
def dapiCellsClass = PathClassFactory.getPathClass('DAPI', makeRGB(0,0,255))
def egfpCellsClass = PathClassFactory.getPathClass('EGFP', makeRGB(0,255,255))
def cy5CellsClass = PathClassFactory.getPathClass('Cy5', makeRGB(255,0,0))
def cy3CellsClass = PathClassFactory.getPathClass('Cy3',  makeRGB(0,255,0))
def cy5cy3CellsClass = PathClassFactory.getPathClass('Cy5-Cy3', makeRGB(255,0,255))

// cells size filter
def toDelete = getDetectionObjects().findAll {
    measurement(it, 'Area µm^2') < 10
}

def stardistDAPI = StarDist2D.builder(pathModel)
          .threshold(0.40)         // Prediction threshold
          .normalizePercentiles(60, 99)     // Percentile normalization
          .pixelSize(1)             // Resolution for detection
          .channels('DAPI') 
          .doLog() 
          .measureShape()                   // Add shape measurements
          .measureIntensity() 
          .build()

def stardistEgfp = StarDist2D.builder(pathModel)
          .threshold(0.40)         // Prediction threshold
          .normalizePercentiles(60, 99)     // Percentile normalization
          .pixelSize(1)             // Resolution for detection
          .channels('EGFP')
          .doLog()
          .measureShape()                   // Add shape measurements
          .measureIntensity()
          .build()


def stardistCy5 = StarDist2D.builder(pathModel)
      .threshold(0.4)             // Prediction threshold
      .normalizePercentiles(60, 99)         // Percentile normalization
      .pixelSize(1)                 // Resolution for detection
      .channels('Cy5') 
      .doLog()
      .measureShape()                       // Add shape measurements
      .measureIntensity()
      .build()

def stardistCy3 = StarDist2D.builder(pathModel)
      .threshold(0.2)              // Prediction threshold
      .normalizePercentiles(60, 99)         // Percentile normalization
      .pixelSize(1)                 // Resolution for detection
      .channels('Cy3') 
      .doLog()
      .measureShape()                       // Add shape measurements
      .measureIntensity() 
      .build()
 
 def coloc(cell1, cell2, cal) {
     def cellColoc = []
     for (c1 in cell1) {
         c1X = c1.getROI().getCentroidX()
         c1Y = c1.getROI().getCentroidY()
         p1 = new Point2(c1X, c1Y)
         for (c2 in cell2) {
             c2X = c2.getROI().getCentroidX()
             c2Y = c2.getROI().getCentroidY()
             p2 = new Point2(c2X, c2Y)
             dist = p1.distance(p2)  * cal.getAveragedPixelSizeMicrons()
             if (dist < 10) {
                 cell = new PathDetectionObject(c1.getROI())
                 cellColoc << cell
             }
         }
     }
     fireHierarchyUpdate()
     return(cellColoc)
 }

// Save annotations
def saveAnnotations(imgName) {
    def path = buildFilePath(imgName + '.annot')
    def annotations = getAnnotationObjects().findAll {it.hasChildren()}
    new File(path).withObjectOutputStream {
        it.writeObject(annotations)
    }
    println('Annotations saved...')
}

// Get objects intensity
def getObjectsIntensity(cells, channel) {
    def measure = channel + ': Mean'
    def means = 0
    def nbCells = cells.size()
    if (nbCells) {
        for (cell in cells) {
            means += cell.getMeasurementList().getMeasurementValue(measure)
        }
        means = means / cells.size()
    }
    return means
}

// loop over images in project

for (entry in project.getImageList()) {
    def imageData = entry.readImageData()
    def server = imageData.getServer()
    def cal = server.getPixelCalibration()
    def pixelWidth = cal.getPixelWidth().doubleValue()
    def pixelUnit = cal.getPixelWidthUnit()
    def imgName = entry.getImageName()
    def imgNameWithOutExt = FilenameUtils.removeExtension(imgName)

    setBatchProjectAndImage(project, imageData)

    // find annotations
    def hierarchy = imageData.getHierarchy()
    def annotations = getAnnotationObjects()
    if (annotations.isEmpty()) {
        Dialogs.showErrorMessage("Problem", "Please create ROIs to analyze")
        return
    }
    def index = 0

    for (an in annotations) {
        index++
        def regionName = "Region_" + index
        if (an.getName() == null)
            an.setName(regionName)
        //def regPathClass = PathClassFactory.getPathClass(regionName, makeRGB(0, 255, 255))
        //an.setPathClass(regPathClass)
        println '-- Analysing ' + an.getName() + ' of image' + imgName + ' --'

        // get annotation area
        def regionArea = an.getROI().getScaledArea(pixelWidth, pixelWidth)
        println an.getName() + ' area = ' + regionArea + ' ' + pixelUnit

        // Do DAPI
        println 'Finding nucleus in ' + an.getName()
        stardistDAPI.detectObjects(imageData, an, true)
        removeObjects(toDelete, true)
        def dapiCells = getDetectionObjects().each{it.setPathClass(dapiCellsClass) &&
            hierarchy.addPathObjectBelowParent(an, it, true)}
        def dapiCellsNb = dapiCells.size()
        println 'DAPI cells in region = ' + dapiCellsNb

        // Do EGFP
        println 'Finding EGFP cells in ' + an.getName()
        stardistEgfp.detectObjects(imageData, an, true)
        removeObjects(toDelete, true)
        def egfpCells = getDetectionObjects().each{it.setPathClass(egfpCellsClass) &&
            hierarchy.addPathObjectBelowParent(an, it, true)}
        def egfpCellsNb = egfpCells.size()
        println 'EGFP cells in region = ' + egfpCellsNb

        // Do Cy3
        println 'Finding cy3 cells in ' + an.getName()
        stardistCy3.detectObjects(imageData, an, true)
        removeObjects(toDelete, true)
        def cy3Cells = getDetectionObjects().each{it.setPathClass(cy3CellsClass) &&
            hierarchy.addPathObjectBelowParent(an, it, true)}
        def cy3CellsNb = cy3Cells.size()
        println 'Cy3 cells in region = ' + cy3CellsNb

        // Do cy5
        println 'Finding cy5 in ' + an.getName()
        stardistCy5.detectObjects(imageData, an, true)
        removeObjects(toDelete, true)
        def cy5Cells = getDetectionObjects().each{it.setPathClass(cy5CellsClass) &&
            hierarchy.addPathObjectBelowParent(an, it, true)}
        def cy5CellsNb = cy5Cells.size()
        println 'Cy5 cells in region = ' + cy5CellsNb

        // Find Cy3 cells colocalized with cy5
        def cy3Cy5CellsNb = 0
        if (cy3CellsNb == 0 || cy5CellsNb == 0) {
            def cy5cy3 = coloc(cy3Cells, cy5Cells, cal)
            cy3Cy5CellsNb = cy5cy3.size()
        }
        print(cy3Cy5CellsNb + ' cy3 cells colocalized with cy5 cells')

        // Find cells means intensities
        def dapiMeanInt = getObjectsIntensity(dapiCells, 'DAPI')
        print('DAPI mean intensity = ' + dapiMeanInt)
        def egfpMeanInt = getObjectsIntensity(egfpCells, 'EGFP')
        print('EGFP mean intensity = ' + egfpMeanInt)
        def cy3MeanInt = getObjectsIntensity(cy3Cells, 'Cy3')
        print('Cy3 mean intensity = ' + cy3MeanInt)
        def cy5MeanInt = getObjectsIntensity(cy5Cells, 'Cy5')
        print('Cy5 mean intensity = ' + cy5MeanInt)

        // Results
        def results = imgNameWithOutExt + '\t' + an.getName() + '\t' + regionArea + '\t' + dapiCellsNb + '\t' + dapiMeanInt + '\t' + egfpCellsNb + '\t' + egfpMeanInt + '\t' + cy3CellsNb + '\t' + cy3MeanInt + '\t' + cy5CellsNb + '\t' + cy5MeanInt + '\t' + cy3Cy5CellsNb + '\n'
        resultsFile << results

        // add detections and save detections
        addObjects(dapiCells)
        addObjects(egfpCells)
        addObjects(cy3Cells)
        addObjects(cy5Cells)

        resolveHierarchy()

        // Save annotations in Shapes format
        saveAnnotations(buildFilePath(resultsDir, imgNameWithOutExt+"_"+an.getName()))
        clearDetections()
    }

}