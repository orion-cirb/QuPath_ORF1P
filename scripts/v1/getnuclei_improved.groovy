//------------------------------ PARAMETERS TO TUNE ------------------------------------------------------------------//
// DAPI cells
minCellSizeDapi = 20

// Cy3 cells
minCellSizeCy3 = 40

// Cy5 cells
minCellSizeCy5 = 40

// All cells
maxCellSize = 600

//--------------------------------------------------------------------------------------------------------------------//


// Imports
import java.util.stream.Collectors

import org.apache.commons.io.FilenameUtils
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.RoiTools
import static qupath.lib.gui.scripting.QPEx.*
import qupath.ext.stardist.StarDist2D
import qupath.lib.objects.*
import qupath.lib.gui.dialogs.Dialogs
import qupath.opencv.ops.ImageOps

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
def resultsDir = buildFilePath(imageDir, '/Results')
if (!fileExists(resultsDir)) mkdirs(resultsDir)
def resultsFile = new File(buildFilePath(resultsDir, 'Results.csv'))
resultsFile.createNewFile()
def resHeaders = 'Image name\tAnnotation name\tArea (um2)\tNb DAPI\tCy3 bg mean intensity\tCy3 bg intensity std\tNb Cy3\tCy3 mean intensity' +
        '\tNb Cy3-DAPI\tCy3-DAPI mean intensity\tCy5 bg mean intensity\tCy5 bg intensity std\tNb Cy5\tCy5 mean intensity' +
        '\tNb Cy5-DAPI\tCy5-DAPI mean intensity\tNb Cy3-Cy5-DAPI\tCy3-Cy5-DAPI mean intensity in Cy3 channel' +
        '\tNb Cy5-Cy3-DAPI\tCy5-Cy3-DAPI mean intensity in Cy5 channel\n'
resultsFile.write(resHeaders)

// Define ClassPaths
def dapiCellsClass = PathClassFactory.getPathClass('DAPI', makeRGB(0,0,255))
def cy5CellsClass = PathClassFactory.getPathClass('Cy5', makeRGB(255,0,0))
def cy3CellsClass = PathClassFactory.getPathClass('Cy3',  makeRGB(255,165,0))

// Build StarDist model
def buildStarDistModel(pathModel, threshold, channel, cellClass) {
    return StarDist2D.builder(pathModel)
            .preprocess(ImageOps.Filters.median(2))
            .threshold(threshold)              // Prediction threshold
            .normalizePercentiles(1, 99)       // Percentile normalization
            .pixelSize(0.5)           // Resolution for detection
            .channels(channel)
            .constrainToParent(false)
            .measureShape()                  // Add shape measurements
            .measureIntensity()              // Add intensity measurements
            .classify(cellClass)
            .build()
}

// Detect cells in a specific annotation and channel
def detectCells(imageData, an, channel, pathModel, probThreshold, cellsClass) {
    println '--- Finding ' + channel + ' cells ---'
    def stardist = buildStarDistModel(pathModel, probThreshold, channel, cellsClass)
    stardist.detectObjects(imageData, an, true)
    def cells = getDetectionObjects().findAll{it.getPathClass() == cellsClass
            && an.getROI().contains(it.getROI().getCentroidX(), it.getROI().getCentroidY())}
    println 'Nb ' + channel + ' cells detected = ' + cells.size()
    stardist.close()
    return cells
}

def getBackground(an, cells) {
    def geometryConverter = new GeometryTools.GeometryConverter.Builder().build()
    def union = geometryConverter.factory.buildGeometry(cells.stream().map(p -> p.getROI().getGeometry()).collect(Collectors.toList())).buffer(0)
    def geometry = an.getROI().getGeometry().difference(union)

    // Create the new ROI
    def bg = PathObjects.createAnnotationObject(GeometryTools.geometryToROI(geometry, an.getROI().getImagePlane()))
    an.addPathObject(bg)
    fireHierarchyUpdate()

    return bg
}

def filterCells(cells, channel, pixelWidth, minCellSize, maxCellSize, bgMean, bgStd) {
    println 'Background mean intensity in ' + channel + ' channel = ' + bgMean + " (std = " + bgStd + ")"
    def filteredCells = cells.findAll{it.getROI().getScaledArea(pixelWidth, pixelWidth) > minCellSize
            && it.getROI().getScaledArea(pixelWidth, pixelWidth) < maxCellSize
            && it.getMeasurementList().getMeasurementValue(channel+': Mean') > (bgMean + bgStd)}
    println 'Nb ' + channel + ' cells remaining = ' + filteredCells.size() + ' (' + (cells.size() - filteredCells.size()) + ' filtered out by shape and intensity)'
    return filteredCells
}

def getIntensityMeasure(bg, channel, measure) {
    return bg.getMeasurementList().getMeasurementValue('ROI: 2.00 µm per pixel: ' + channel + ': ' + measure)
}

// Get colocalized cells among two cell populations
def coloc(cell1, cell2, colocParam) {
    def tool = new RoiTools()
    def cellColoc = []
    if (cell1.size() != 0 && cell2.size() != 0) {
        for (c1 in cell1) {
            def roiC1 = c1.getROI()
            for (c2 in cell2) {
                def roiC2 = c2.getROI()
                if (tool.areaContains(roiC1, roiC2.getCentroidX(), roiC2.getCentroidY())) {
                    if (colocParam)
                        c1.getMeasurementList().addMeasurement("Cy3-Cy5-DAPI colocalization", 1)
                    cellColoc << PathObjectTools.transformObject(c1, null, true)
                    break
                }
            }
        }
    }
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
    def annotations = getAnnotationObjects()
    if (annotations.isEmpty()) {
        Dialogs.showErrorMessage("Problem", "Please create ROIs to analyze in image " + imgName)
        continue
    }

    def index = 0
    for (an in annotations) {
        index++
        if (an.getName() == null)
            an.setName("Region_" + index)
        println ''
        println '------ Analyzing ROI ' + an.getName() + ' of image ' + imgName + ' ------'

        // Detect cells in each channel
        clearAllObjects()
        addObject(an)
        def dapiCells = detectCells(imageData, an, 'DAPI', pathModel, 0.6, dapiCellsClass)
        def cy3Cells = detectCells(imageData, an, 'Cy3', pathModel, 0.6, cy3CellsClass)
        def cy5Cells = detectCells(imageData, an, 'Cy5', pathModel, 0.6, cy5CellsClass)

        // Get background corresponding to current annotation for each channel
        def dapiBg = getBackground(an, dapiCells)
        def cy3Bg = getBackground(an, cy3Cells)
        def cy5Bg = getBackground(an, cy5Cells)

        deselectAll()
        selectObjects([dapiBg, cy3Bg, cy5Bg])
        runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons": 2.0,  "region": "ROI",  "tileSizeMicrons": 25.0,  "channel1": true,  ' +
                  '"channel2": true,  "channel3": true,  "channel4": true,  "doMean": true,  "doStdDev": true,  "doMinMax": false,  "doMedian": false,  "doHaralick": false}')

        // Filter cells by shape and intensity
        dapiCells = filterCells(dapiCells, 'DAPI', pixelWidth, minCellSizeDapi, maxCellSize,
                getIntensityMeasure(dapiBg, 'DAPI', 'Mean'), getIntensityMeasure(dapiBg, 'DAPI', 'Std.dev.'))
        cy3Cells = filterCells(cy3Cells, 'Cy3', pixelWidth, minCellSizeCy3, maxCellSize,
                getIntensityMeasure(dapiBg, 'Cy3', 'Mean'), getIntensityMeasure(dapiBg, 'Cy3', 'Std.dev.'))
        cy5Cells = filterCells(cy5Cells, 'Cy5', pixelWidth, minCellSizeCy5, maxCellSize,
                getIntensityMeasure(dapiBg, 'Cy5', 'Mean'), getIntensityMeasure(dapiBg, 'Cy5', 'Std.dev.'))

        println '--- Colocalization ---'
        // Find Cy3 cells colocalized with DAPI nuclei
        def cy3DapiCells = coloc(cy3Cells, dapiCells, false)
        print(cy3DapiCells.size() + '/' + cy3Cells.size() + ' Cy3 cells colocalized with DAPI nuclei')

        // Find Cy5 cells colocalized with DAPI nuclei
        def cy5DapiCells = coloc(cy5Cells, dapiCells, false)
        print(cy5DapiCells.size() + '/' + cy5Cells.size() + ' Cy5 cells colocalized with DAPI nuclei')

        /*// Find Cy3-DAPI cells colocalized with Cy5-DAPI cells
        def cy3Cy5Cells = coloc(cy3DapiCells, cy5DapiCells, true)
        print(cy3Cy5Cells.size() + '/' + cy3DapiCells.size() + ' Cy3-DAPI cells colocalized with Cy5-DAPI cells')

        // Find Cy5-DAPI cells colocalized with Cy3-DAPI cells
        def cy5Cy3Cells = coloc(cy5DapiCells, cy3DapiCells, true)
        print(cy5Cy3Cells.size() + '/' + cy5DapiCells.size() + ' Cy5-DAPI cells colocalized with Cy3-DAPI cells')

        // Save results
        def results = imgNameWithOutExt + '\t' + an.getName() + '\t' + an.getROI().getScaledArea(pixelWidth, pixelWidth) + '\t' + dapiCells.size() +
                '\t' + bgCy3Int + '\t' + cy3Cells.size() + '\t' + getObjectsIntensity(cy3Cells, 'Cy3', bgCy3Int) + '\t' + cy3DapiCells.size() + '\t' + getObjectsIntensity(cy3DapiCells, 'Cy3', bgCy3Int) +
                '\t' + bgCy5Int + '\t' + cy5Cells.size() + '\t' + getObjectsIntensity(cy5Cells, 'Cy5', bgCy5Int) + '\t' + cy5DapiCells.size() + '\t' + getObjectsIntensity(cy5DapiCells, 'Cy5', bgCy5Int) +
                '\t' + cy3Cy5Cells.size() + '\t' + getObjectsIntensity(cy3Cy5Cells, 'Cy3', bgCy3Int) + '\t' + cy5Cy3Cells.size() + '\t' + getObjectsIntensity(cy5Cy3Cells, 'Cy5', bgCy5Int) + '\n'
        resultsFile << results*/

        an.clearPathObjects()
        an.addPathObjects(dapiCells)
        an.addPathObjects(cy3DapiCells)
        an.addPathObjects(cy5DapiCells)
        fireHierarchyUpdate()

        // Save detections
        clearAllObjects()
        addObject(an)
        saveAnnotations(buildFilePath(resultsDir, imgNameWithOutExt+"_"+an.getName()))
        println ''
    }
    return
}
