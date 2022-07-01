import org.apache.commons.io.FilenameUtils
import qupath.lib.objects.classes.PathClassFactory

import static qupath.lib.gui.scripting.QPEx.*
import qupath.ext.stardist.StarDist2D
import qupath.lib.objects.*
import qupath.lib.geom.Point2
import qupath.lib.gui.dialogs.Dialogs

// Init project
setImageType('Fluorescence')
def project = getProject()
def pathProject = buildFilePath(PROJECT_BASE_DIR)
def pathModel = buildFilePath(pathProject,'models','dsb2018_heavy_augment.pb')
if (pathModel == null) {
    Dialogs.showErrorMessage("Problem", 'No StarDist model found')
    return
}
def imageDir = new File(project.getImageList()[0].getUris()[0]).getParent()

// Create results file and write headers
def resultsDir = buildFilePath(imageDir, '/../Results')
if (!fileExists(resultsDir)) mkdirs(resultsDir)
def resultsFile = new File(buildFilePath(resultsDir, 'Results.csv'))
resultsFile.createNewFile()
def resHeaders = 'Image name\tAnnotation name\tArea (um2)\tNb DAPI\tNb EGFP\tEGFP mean intensity' +
        '\tNb EGFP-DAPI\tEGFP-DAPI mean intensity\tNb Cy3\tCy3 mean intensity\tNb Cy3-DAPI\tCy3-DAPI mean intensity\tNb Cy5\tCy5 mean intensity' +
        '\tNb Cy5-DAPI\tCy5-DAPI mean intensity\tNb Cy3-Cy5-DAPI\n'
resultsFile.write(resHeaders)

// Define ClassPaths
def dapiCellsClass = PathClassFactory.getPathClass('DAPI', makeRGB(0,0,255))
def egfpCellsClass = PathClassFactory.getPathClass('EGFP', makeRGB(0,255,0))
def cy5CellsClass = PathClassFactory.getPathClass('Cy5', makeRGB(255,0,0))
def cy3CellsClass = PathClassFactory.getPathClass('Cy3',  makeRGB(255,165,0))

// Build StarDist model
def buildStarDistModel(pathModel, threshold, channel, cellClass) {
    return StarDist2D.builder(pathModel)
            .threshold(threshold)              // Prediction threshold
            .normalizePercentiles(1, 99)       // Percentile normalization
            .pixelSize(0.5)           // Resolution for detection
            .channels(channel)
            .measureShape()                  // Add shape measurements
            .measureIntensity()              // Add intensity measurements
            .classify(cellClass)
            .build()
}

// Detect cells in a specific annotation and channel
def detectCells(imageData, an, channel, hierarchy, pathModel, probThreshold, cellsClass, pixelWidth, minCellSize, bgInt, minIntensityPercentage) {
    println '--- Finding ' + channel + ' cells ---'
    println 'Background median intensity in ' + channel + ' channel = ' + bgInt
    def stardist = buildStarDistModel(pathModel, probThreshold, channel, cellsClass)
    stardist.detectObjects(imageData, an, true)
    def cells = getDetectionObjects().findAll{it.getPathClass() == cellsClass
            && it.getROI().getScaledArea(pixelWidth, pixelWidth) > minCellSize
            && it.getMeasurementList().getMeasurementValue(channel +': Median') > (minIntensityPercentage*bgInt)}
    cells.each{hierarchy.addPathObjectBelowParent(an, it, true)}
    println 'Nb ' + channel + ' cells = ' + cells.size() + ' (' + (getDetectionObjects().findAll{it.getPathClass() == cellsClass}.size() - cells.size()) + ' filtered out)'
    return cells
}

// Get colocalized cells among two cell populations
 def coloc(cell1, cell2, cal) {
     def cellColoc = []
     for (c1 in cell1) {
         def roiC1 = c1.getROI()
         def p1 = new Point2(roiC1.getCentroidX(), roiC1.getCentroidY())
         for (c2 in cell2) {
             def roiC2 = c2.getROI()
             def p2 = new Point2(roiC2.getCentroidX(), roiC2.getCentroidY())
             def dist = p1.distance(p2)  * cal.getAveragedPixelSizeMicrons()
             if (dist < 5) {
                 cellColoc << c1
             }
         }
     }
     fireHierarchyUpdate()
     return(cellColoc)
 }

// Get mean intensity of a population of objects
def getObjectsIntensity(cells, channel, bgInt) {
    def measure = channel + ': Mean'
    def means = 0
    def nbCells = cells.size()
    if (nbCells) {
        for (cell in cells) {
            means += cell.getMeasurementList().getMeasurementValue(measure)
        }
        means /= cells.size()
        means -= bgInt
    }
    return means
}

// Save annotations
def saveAnnotations(imgName) {
    def path = buildFilePath(imgName + '.annot')
    def annotations = getAnnotationObjects()
    new File(path).withObjectOutputStream {
        it.writeObject(annotations)
    }
    println('Annotations saved')
}

// Loop over images in project
for (entry in project.getImageList()) {
    def imageData = entry.readImageData()
    def server = imageData.getServer()
    def cal = server.getPixelCalibration()
    def pixelWidth = cal.getPixelWidth().doubleValue()
    def imgName = entry.getImageName()
    def imgNameWithOutExt = FilenameUtils.removeExtension(imgName)
    setBatchProjectAndImage(project, imageData)
    println ''
    println ''
    println '--------- ANALYZING IMAGE ' + imgName + ' ---------'

    // Find annotations
    def hierarchy = imageData.getHierarchy()
    def allAnnotations = getAnnotationObjects()
    def  annotations = allAnnotations.findAll{! it.getName().contains('_bg')}
    if (annotations.isEmpty()) {
        Dialogs.showErrorMessage("Problem", "Please create ROIs to analyze in image " + imgName)
        continue
    }
    def  backgrounds = allAnnotations.findAll{it.getName().contains('_bg')}
    if (backgrounds.isEmpty() || backgrounds.size() < annotations.size()) {
        Dialogs.showErrorMessage("Problem", "Please provide background ROI for each ROI to analyze in image " + imgName)
        continue
    }
    selectAnnotations()
    runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons": 2.0,  "region": "ROI",  "tileSizeMicrons": 25.0,  "channel1": true,  ' +
            '"channel2": true,  "channel3": true,  "channel4": true,  "doMean": false,  "doStdDev": false,  "doMinMax": false,  "doMedian": true,  "doHaralick": false}')

    def index = 0
    for (an in annotations) {
        index++
        if (an.getName() == null)
            an.setName("Region_" + index)
        println ''
        println '------ Analyzing ROI ' + an.getName() + ' of image ' + imgName + ' ------'

        // Get background corresponding to current annotation
        def bg = backgrounds.find{it -> it.getName().contains(an.getName())}
        def bgDapiInt = bg.getMeasurementList().getMeasurementValue('ROI: 2.00 µm per pixel: DAPI: Median')
        def bgEgfpInt = bg.getMeasurementList().getMeasurementValue('ROI: 2.00 µm per pixel: EGFP: Median')
        def bgCy3Int = bg.getMeasurementList().getMeasurementValue('ROI: 2.00 µm per pixel: Cy3: Median')
        def bgCy5Int = bg.getMeasurementList().getMeasurementValue('ROI: 2.00 µm per pixel: Cy5: Median')

        // Detect cells in each channel
        def dapiCells = detectCells(imageData, an, 'DAPI', hierarchy, pathModel, 0.6, dapiCellsClass,
                pixelWidth, 20, bgDapiInt, 1.2)
        def egfpCells = detectCells(imageData, an, 'EGFP', hierarchy, pathModel, 0.4, egfpCellsClass,
                pixelWidth, 20, bgEgfpInt, 1.4)
        def cy3Cells = detectCells(imageData, an, 'Cy3', hierarchy, pathModel, 0.4, cy3CellsClass,
                pixelWidth, 40, bgCy3Int, 1.4)
        def cy5Cells = detectCells(imageData, an, 'Cy5', hierarchy, pathModel, 0.6, cy5CellsClass,
                pixelWidth, 40, bgCy5Int, 1.2)

        println '--- Colocalization ---'
        // Find EGFP cells colocalized with DAPI nuclei
        def egfpDapiCells = coloc(egfpCells, dapiCells, cal)
        print(egfpDapiCells.size() + '/' + egfpCells.size() + ' EGFP cells colocalized with DAPI nuclei')

        // Find Cy3 cells colocalized with DAPI nuclei
        def cy3DapiCells = coloc(cy3Cells, dapiCells, cal)
        print(cy3DapiCells.size() + '/' + cy3Cells.size() + ' Cy3 cells colocalized with DAPI nuclei')

        // Find Cy5 cells colocalized with DAPI nuclei
        def cy5DapiCells = coloc(cy5Cells, dapiCells, cal)
        print(cy5DapiCells.size() + '/' + cy5Cells.size() + ' Cy5 cells colocalized with DAPI nuclei')

        // Find Cy3-DAPI cells colocalized with Cy5-DAPI cells
        def cy3Cy5Cells = coloc(cy3DapiCells, cy5DapiCells, cal)
        print(cy3Cy5Cells.size() + '/' + cy3Cells.size() + ' Cy3-DAPI cells colocalized with Cy5-DAPI cells')

        // Save results
        def results = imgNameWithOutExt + '\t' + an.getName() + '\t' + an.getROI().getScaledArea(pixelWidth, pixelWidth) + '\t' + dapiCells.size() +
                '\t' + egfpCells.size() + '\t'+ getObjectsIntensity(egfpCells, 'EGFP', bgEgfpInt) + '\t' + egfpDapiCells.size() + '\t' + getObjectsIntensity(egfpDapiCells, 'EGFP', bgEgfpInt) +
                '\t' + cy3Cells.size() + '\t' + getObjectsIntensity(cy3Cells, 'Cy3', bgCy3Int) + '\t' + cy3DapiCells.size() + '\t' + getObjectsIntensity(cy3DapiCells, 'Cy3', bgCy3Int) +
                '\t' + cy5Cells.size() + '\t' + getObjectsIntensity(cy5Cells, 'Cy5', bgCy5Int) + '\t' + cy5DapiCells.size() + '\t' + getObjectsIntensity(cy5DapiCells, 'Cy5', bgCy5Int) +
                '\t' + cy3Cy5Cells.size() + '\n'
        resultsFile << results

        // Save detections
        clearDetections()
        addObject(bg)
        addObjects(dapiCells)
        addObjects(egfpDapiCells)
        addObjects(cy3DapiCells)
        addObjects(cy5DapiCells)
        resolveHierarchy()
        saveAnnotations(buildFilePath(resultsDir, imgNameWithOutExt+"_"+an.getName()))
        println ''
    }
}