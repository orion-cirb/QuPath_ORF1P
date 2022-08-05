//------------------------------ PARAMETERS TO TUNE ------------------------------------------------------------------//
// DAPI cells
minCellSizeDapi = 20
minIntensityDapi = 1.2 // DAPI cell mean intensity > minIntensityDapi * background mean intensity in DAPI channel

// EGFP cells
minCellSizeEgfp = 20
minIntensityEgfp = 1.4

// Cy3 cells
minCellSizeCy3 = 40
minIntensityCy3 = 1.4

// Cy5 cells
minCellSizeCy5 = 40
minIntensityCy5 = 1.2

// Max cells size
maxCellSize = 600

//--------------------------------------------------------------------------------------------------------------------//


// Imports
import org.apache.commons.io.FilenameUtils
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.roi.RoiTools

import static qupath.lib.gui.scripting.QPEx.*
import qupath.ext.stardist.StarDist2D
import qupath.lib.objects.*
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
def resultsDir = buildFilePath(imageDir, '/Results')
if (!fileExists(resultsDir)) mkdirs(resultsDir)
def resultsFile = new File(buildFilePath(resultsDir, 'Results.csv'))
resultsFile.createNewFile()
def resHeaders = 'Image name\tAnnotation name\tArea (um2)\tNb DAPI\tEGFP bg median intensity\tNb EGFP\tEGFP mean intensity' +
        '\tNb EGFP-DAPI\tEGFP-DAPI mean intensity\tCy3 bg median intensity\tNb Cy3\tCy3 mean intensity\tNb Cy3-DAPI\tCy3-DAPI mean intensity' +
        '\tCy5 bg median intensity\tNb Cy5\tCy5 mean intensity\tNb Cy5-DAPI\tCy5-DAPI mean intensity\tNb Cy3-Cy5-DAPI' +
        '\tCy3-Cy5-DAPI mean intensity in Cy3 channel\tNb Cy5-Cy3-DAPI\tCy5-Cy3-DAPI mean intensity in Cy5 channel' +
        '\tNb Cy3-EGFP-DAPI\t\'Nb Cy5-EGFP-DAPI\n'
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
            .constrainToParent(false)
            .measureShape()                  // Add shape measurements
            .measureIntensity()              // Add intensity measurements
            .classify(cellClass)
            .build()
}

// Detect cells in a specific annotation and channel
def detectCells(imageData, an, channel, pathModel, probThreshold, cellsClass, pixelWidth, minCellSize, bgInt, minIntensityPercentage) {
    println '--- Finding ' + channel + ' cells ---'
    println 'Background median intensity in ' + channel + ' channel = ' + bgInt
    def stardist = buildStarDistModel(pathModel, probThreshold, channel, cellsClass)
    stardist.detectObjects(imageData, an, true)
    def cells = getDetectionObjects().findAll{it.getPathClass() == cellsClass
            && it.getROI().getScaledArea(pixelWidth, pixelWidth) > minCellSize
            && it.getROI().getScaledArea(pixelWidth, pixelWidth) < maxCellSize
            && it.getMeasurementList().getMeasurementValue(channel +': Median') > (minIntensityPercentage*bgInt)
            && an.getROI().contains(it.getROI().getCentroidX(), it.getROI().getCentroidY())}
    println 'Nb ' + channel + ' cells = ' + cells.size() + ' (' + (getDetectionObjects().findAll{it.getPathClass() == cellsClass}.size() - cells.size()) + ' filtered out)'
    return cells
}

// Get colocalized cells among two cell populations
def coloc(cell1, cell2, colocParam) {
    def tool = new RoiTools()
    def cellColoc = []
    if (cell1.size() !=0 && cell2 !=0) {
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
        clearAllObjects()
        addObject(an)
        def dapiCells = detectCells(imageData, an, 'DAPI', pathModel, 0.6, dapiCellsClass,
                pixelWidth, minCellSizeDapi, bgDapiInt, minIntensityDapi)
        def egfpCells = detectCells(imageData, an, 'EGFP', pathModel, 0.4, egfpCellsClass,
                pixelWidth, minCellSizeEgfp, bgEgfpInt, minIntensityEgfp)
        def cy3Cells = detectCells(imageData, an, 'Cy3', pathModel, 0.4, cy3CellsClass,
                pixelWidth, minCellSizeCy3, bgCy3Int, minIntensityCy3)
        def cy5Cells = detectCells(imageData, an, 'Cy5', pathModel, 0.6, cy5CellsClass,
                pixelWidth, minCellSizeCy5, bgCy5Int, minIntensityCy5)

        println '--- Colocalization ---'
        // Find EGFP cells colocalized with DAPI nuclei
        def egfpDapiCells = coloc(egfpCells, dapiCells, false)
        print(egfpDapiCells.size() + '/' + egfpCells.size() + ' EGFP cells colocalized with DAPI nuclei')

        // Find Cy3 cells colocalized with DAPI nuclei
        def cy3DapiCells = coloc(cy3Cells, dapiCells, false)
        print(cy3DapiCells.size() + '/' + cy3Cells.size() + ' Cy3 cells colocalized with DAPI nuclei')

        // Find Cy5 cells colocalized with DAPI nuclei
        def cy5DapiCells = coloc(cy5Cells, dapiCells, false)
        print(cy5DapiCells.size() + '/' + cy5Cells.size() + ' Cy5 cells colocalized with DAPI nuclei')

        // Find Cy3-DAPI cells colocalized with Cy5-DAPI cells
        def cy3Cy5Cells = coloc(cy3DapiCells, cy5DapiCells, true)
        print(cy3Cy5Cells.size() + '/' + cy3DapiCells.size() + ' Cy3-DAPI cells colocalized with Cy5-DAPI cells')

        // Find Cy5-DAPI cells colocalized with Cy3-DAPI cells
        def cy5Cy3Cells = coloc(cy5DapiCells, cy3DapiCells, true)
        print(cy5Cy3Cells.size() + '/' + cy5DapiCells.size() + ' Cy5-DAPI cells colocalized with Cy3-DAPI cells')
        // ????
        //cy5Cy3Cells.get(0).getMeasurementList().getMeasurementValue("Coloc")
        //cy5Cy3Cells.get(1).getMeasurementList().getMeasurementValue("Coloc")
        //cy5Cy3Cells.get(2).getMeasurementList().getMeasurementValue("Coloc")

        // Find EGFP-DAPI cells colocalized with Cy3-DAPI cells
        def egfpCy3Cells = coloc(egfpDapiCells, cy3DapiCells, false)
        print(egfpCy3Cells.size() + '/' + egfpDapiCells.size() + ' EGFP-DAPI cells colocalized with Cy3-DAPI cells')

        // Find EGFP-DAPI cells colocalized with Cy3-DAPI cells
        def egfpCy5Cells = coloc(egfpDapiCells, cy5DapiCells, false)
        print(egfpCy5Cells.size() + '/' + egfpDapiCells.size() + ' EGFP-DAPI cells colocalized with Cy5-DAPI cells')

        // Save results
        def results = imgNameWithOutExt + '\t' + an.getName() + '\t' + an.getROI().getScaledArea(pixelWidth, pixelWidth) + '\t' + dapiCells.size() +
                '\t' + bgEgfpInt + '\t' + egfpCells.size() + '\t'+ getObjectsIntensity(egfpCells, 'EGFP', bgEgfpInt) + '\t' + egfpDapiCells.size() + '\t' + getObjectsIntensity(egfpDapiCells, 'EGFP', bgEgfpInt) +
                '\t' + bgCy3Int + '\t' + cy3Cells.size() + '\t' + getObjectsIntensity(cy3Cells, 'Cy3', bgCy3Int) + '\t' + cy3DapiCells.size() + '\t' + getObjectsIntensity(cy3DapiCells, 'Cy3', bgCy3Int) +
                '\t' + bgCy5Int + '\t' + cy5Cells.size() + '\t' + getObjectsIntensity(cy5Cells, 'Cy5', bgCy5Int) + '\t' + cy5DapiCells.size() + '\t' + getObjectsIntensity(cy5DapiCells, 'Cy5', bgCy5Int) +
                '\t' + cy3Cy5Cells.size() + '\t' + getObjectsIntensity(cy3Cy5Cells, 'Cy3', bgCy3Int) + '\t' + cy5Cy3Cells.size() + '\t' + getObjectsIntensity(cy5Cy3Cells, 'Cy5', bgCy5Int) +
                '\t' + egfpCy3Cells.size() + '\t' + egfpCy5Cells.size()  + '\n'
        resultsFile << results

        an.clearPathObjects()
        an.addPathObjects(dapiCells)
        an.addPathObjects(egfpDapiCells)
        an.addPathObjects(cy3DapiCells)
        an.addPathObjects(cy5DapiCells)
        fireHierarchyUpdate()

        // Save detections
        clearAllObjects()
        addObject(bg)
        addObject(an)
        saveAnnotations(buildFilePath(resultsDir, imgNameWithOutExt+"_"+an.getName()))
        println ''
    }
}
