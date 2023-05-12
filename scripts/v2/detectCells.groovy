//------------------------------------------- PARAMETERS TO TUNE -----------------------------------------------------//
// NUCLEI DETECTION WITH STARDIST
def probabilityThreshold = 0.75 // between 0 and 1: if you decrease it, you will detect more cells

// SIZE FILTERING
// DAPI nuclei
def dapiMinArea = 20 // in µm2
def dapiMaxArea = 200
// NeuN cells
def cy5MinArea = 40
def cy5MaxArea = 400
// ORF1p cells
def cy3MinArea = 40
def cy3MaxArea = 400

// INTENSITY FILTERING
// Cell mean int > (C1 * bg mean int + C2 * bg int SD)
// DAPI nuclei
def dapiC1 = 0
def dapiC2 = 0
// NeuN cells
def cy5C1 = 0
def cy5C2 = 0
// ORF1p cells
def cy3C1 = 0
def cy3C2 = 0


//----------------------------------------------- PIPELINE -----------------------------------------------------------//
// Imports
import org.apache.commons.io.FilenameUtils
import qupath.lib.objects.classes.PathClass
import qupath.lib.gui.dialogs.Dialogs
import java.util.stream.Collectors
import qupath.lib.roi.GeometryTools
import qupath.lib.roi.RoiTools
import qupath.ext.stardist.StarDist2D
import qupath.ext.biop.cellpose.Cellpose2D
import qupath.lib.objects.*
import qupath.opencv.ops.ImageOps

// Avoid TopologyException errors: see https://forum.image.sc/t/stardist-error-message-topologyexception/67708/6
org.locationtech.jts.geom.GeometryOverlay.isOverlayNG = true
if (!org.locationtech.jts.geom.GeometryOverlay.isOverlayNG) {
    Dialogs.showErrorMessage('Problem', 'JTS option avoiding TopologyException errors not set')
    return
}

// Init project
setImageType('Fluorescence')
def project = getProject()
def pathProject = buildFilePath(PROJECT_BASE_DIR)
def stardistPathModel = buildFilePath(pathProject, 'models', 'dsb2018_heavy_augment.pb')
if (stardistPathModel == null) {
    Dialogs.showErrorMessage('Problem', 'StarDist model not found')
    return
}
def imageDir = new File(project.getImageList()[0].getURIs()[0]).getParent()


// Create results files and write headers
def resultsDir = buildFilePath(imageDir, '/Results')
if (!fileExists(resultsDir)) mkdirs(resultsDir)

def globalResultsFile = new File(buildFilePath(resultsDir, 'globalResults.csv'))
globalResultsFile.createNewFile()
def globalHeaders = 'Image name\tRegion name\tRegion area\tDAPI bg int mean\tDAPI bg int sd' +
        '\tCy5 bg int mean\tCy5 bg int sd\tCy3 bg int mean\tCy3 bg int sd' +
        '\tNb DAPI+\tNb DAPI+ Cy5+\tNb DAPI+ Cy3+\tNb DAPI+ Cy5- Cy3-\tNb DAPI+ Cy5+ Cy3-\tNb DAPI+ Cy5- Cy3+\tNb DAPI+ Cy5+ Cy3+\n'
globalResultsFile.write(globalHeaders)

def cellsResultsFile = new File(buildFilePath(resultsDir, 'cellsResults.csv'))
cellsResultsFile.createNewFile()
def cellsHeaders = 'Image name\tRegion name\tNuc area\tNuc circularity\t' +
        'Nuc DAPI int mean\tNuc DAPI int sd\tNuc Cy5 int mean\tNuc Cy5 int sd\tNuc Cy3 int mean\tNuc Cy3 int sd\t' +
        'is Cy5?\tCy5 area\tCy5 int mean\tCy5 int sd\tCy5 Cy3 int mean\tCy5 Cy3 int sd\t' +
        'is Cy3?\tCy3 area\tCy3 int mean\tCy3 int sd\tCy3 Cy5 int mean\tCy3 Cy5 int sd\n'
cellsResultsFile.write(cellsHeaders)


// Define ClassPaths
def dapiCellsClass = PathClass.fromString('DAPI', makeRGB(0, 0, 255))
def cy5CellsClass = PathClass.fromString('Cy5', makeRGB(255,0,0))
def cy3CellsClass = PathClass.fromString('Cy3',  makeRGB(255,165,0))

// Loop over images in project
for (entry in project.getImageList()) {
    def imageData = entry.readImageData()
    def cal = imageData.getServer().getPixelCalibration()
    def pixelWidth = cal.getPixelWidth().doubleValue()
    def imgName = entry.getImageName()
    def imgNameWithOutExt = FilenameUtils.removeExtension(imgName)
    setBatchProjectAndImage(project, imageData)
    println ''
    println ''
    println '------ ANALYZING IMAGE ' + imgName + ' ------'

    // Find annotations
    def annotations = getAnnotationObjects()
    if (annotations.isEmpty()) {
        Dialogs.showErrorMessage("Problem", "No region to analyze in image " + imgName)
        continue
    }

    def anNb = annotations.size()
    def anId = 1
    for (an in annotations) {
        def regionName = an.getName().replaceAll("/", "-").replaceAll(",", "-").replaceAll(" ", "")
        if(an.hasChildObjects()) {
            println '--- Skipping region ' + regionName + ' ('+anId+'/'+anNb+') ---'
        } else {
            println ''
            println '--- Analyzing region ' + regionName + ' ('+anId+'/'+anNb+') ---'

            // Get region parents
            def parentsName = regionName
            def parent = an.getParent()
            while (parent != null && !parent.isRootObject()) {
                parentsName = parent.getName() + '_' + parentsName
                parent = parent.getParent()
            }
            println 'Position in hierarchy: ' + parentsName

            // Detect nuclei and cells
            clearAllObjects()
            addObject(an)
            def dapiCells = detectNuclei(imageData, an, stardistPathModel, 'DAPI', probabilityThreshold, dapiCellsClass)
            def cy5Cells = detectCells(imageData, an, 'cyto2', 'Cy5', 30, 60, cy5CellsClass)
            def cy3Cells = detectCells(imageData, an, 'cyto2', 'Cy3', 30, 60, cy3CellsClass)

            // Compute background for each channel: difference(current region, union(detected cells in channel))
            def dapiBg = getBackground(an, dapiCells)
            def cy5Bg = getBackground(an, cy5Cells)
            def cy3Bg = getBackground(an, cy3Cells)

            deselectAll()
            selectObjects([dapiBg, cy5Bg, cy3Bg])
            runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons": 2.0,  "region": "ROI",  "tileSizeMicrons": 25.0,  "channel1": true,  ' +
                    '"channel2": false,  "channel3": true,  "channel4": true,  "doMean": true,  "doStdDev": true,  "doMinMax": false,  "doMedian": false,  "doHaralick": false}')

            // Get background intensity (mean and SD) for each channel
            def dapiBgInt = getBackgroundIntensity(dapiBg, 'DAPI')
            def dapiBgMean = dapiBgInt[0], dapiBgSD = dapiBgInt[1]
            def cy5BgInt = getBackgroundIntensity(cy5Bg, 'Cy5')
            def cy5BgMean = cy5BgInt[0], cy5BgSD = cy5BgInt[1]
            def cy3BgInt = getBackgroundIntensity(cy3Bg, 'Cy3')
            def cy3BgMean = cy3BgInt[0], cy3BgSD = cy3BgInt[1]

            // Filter cells by size and intensity
            dapiCells = filterCells(dapiCells, 'DAPI', pixelWidth, dapiMinArea, dapiMaxArea, dapiBgMean, dapiBgSD, dapiC1, dapiC2)
            cy5Cells = filterCells(cy5Cells, 'Cy5', pixelWidth, cy5MinArea, cy5MaxArea, cy5BgMean, cy5BgSD, cy5C1, cy5C2)
            cy3Cells = filterCells(cy3Cells, 'Cy3', pixelWidth, cy3MinArea, cy3MaxArea, cy3BgMean, cy3BgSD, cy3C1, cy3C2)

            // Perform colocalization between nuclei and cells
            def cells = colocalization(dapiCells, cy5Cells, cy3Cells)

            // Save results
            println '- Saving results -'
            for (cell in cells) {
                def nucleusParams = cell.getNucleusParams(dapiBgMean, cy5BgMean, cy3BgMean)
                def cy5Params = cell.getCy5Params(cy5BgMean, cy3BgMean)
                def cy3Params = cell.getCy3Params(cy3BgMean, cy5BgMean)

                def cellResults = imgNameWithOutExt + '\t' + parentsName + '\t' + nucleusParams[0] + '\t' + nucleusParams[1] +
                        '\t' + nucleusParams[2] + '\t' + nucleusParams[3] + '\t' + nucleusParams[4] + '\t' + nucleusParams[5] + '\t' + nucleusParams[6] + '\t' + nucleusParams[7] +
                        '\t' + cy5Params[0] + '\t' + cy5Params[1] + '\t' + cy5Params[2] + '\t' + cy5Params[3] + '\t' + cy5Params[4] + '\t' + cy5Params[5] +
                        '\t' + cy3Params[0] + '\t' + cy3Params[1] + '\t' + cy3Params[2] + '\t' + cy3Params[3] + '\t' + cy3Params[4] + '\t' + cy3Params[5] + '\n'
                cellsResultsFile << cellResults
            }

            def cellsCount = countCells(cells)
            def globalResults = imgNameWithOutExt + '\t' + parentsName + '\t' + an.getROI().getScaledArea(pixelWidth, pixelWidth).round(3) +
                    '\t' + dapiBgMean.round(3) + '\t' + dapiBgSD.round(3) + '\t' + cy5BgMean.round(3) + '\t' + cy5BgSD.round(3) + '\t' +
                    cy3BgMean.round(3) + '\t' + cy3BgSD.round(3) + '\t' + dapiCells.size() + '\t' + cellsCount[0] + '\t' + cellsCount[1] +
                    '\t' + cellsCount[2] + '\t' + cellsCount[3] + '\t' + cellsCount[4] + '\t' + cellsCount[5] + '\n'
            globalResultsFile << globalResults

            // Save detections
            an.clearChildObjects()
            an.addChildObjects(dapiCells)
            an.addChildObjects(cy5Cells.findAll{it.retrieveMetadataValue("Colocalized") == "true"})
            an.addChildObjects(cy3Cells.findAll{it.retrieveMetadataValue("Colocalized") == "true"})
            fireHierarchyUpdate()

            clearAllObjects()
            addObject(an)
            saveAnnotations(buildFilePath(resultsDir, imgNameWithOutExt+"_"+regionName))
            println ''
        }
        anId++
    }
}
println '------ ANALYSIS DONE ------'

//------------------------------------------------- UTILS ------------------------------------------------------------//

// Build StarDist model
def buildStarDistModel(pathModel, channel, probThreshold, cellClass) {
    return StarDist2D.builder(pathModel)
            .preprocess(ImageOps.Filters.median(1))     // List of preprocessing ImageOps to run on the images before exporting them
            .normalizePercentiles(1, 99)       // Percentile normalization
            .pixelSize(0.5)            // Resolution for detection
            .channels(channel)                 // Select detection channel
            .threshold(probThreshold)          // Prediction threshold
            .simplify(0)
            .constrainToParent(false)
            .measureShape()                    // Add shape measurements
            .measureIntensity()                // Add intensity measurements
            .classify(cellClass)
            .build()
}

// Build Cellpose model
def buildCellposeModel(pathModel, channel, diameter, overlap, cellClass) {
    return Cellpose2D.builder(pathModel)
            .preprocess(ImageOps.Filters.median(2))     // List of preprocessing ImageOps to run on the images before exporting them
            .normalizePercentiles(1, 99)                // Percentile normalization
            .pixelSize(0.5)                             // Resolution for detection
            .channels(channel)                          // Select detection channel(s)
//          .maskThreshold(-0.2)                        // Threshold for the mask detection, defaults to 0.0
//          .flowThreshold(0.5)                         // Threshold for the flows, defaults to 0.4
            .diameter(diameter)                         // Median object diameter. Set to 0.0 for the `bact_omni` model or for automatic computation
            .setOverlap(overlap)                             // Overlap between tiles (in pixels) that the QuPath Cellpose Extension will extract. Defaults to 2x the diameter or 60 px if the diameter is set to 0
            .simplify(1.5)
            .constrainToParent(false)
            .measureShape()                             // Add shape measurements
            .measureIntensity()                         // Add intensity measurements (in all compartments)
            .classify(cellClass)
//          .doLog()
            .build()
}


// Detect nuclei in a specific annotation and channel with Stardist
def detectNuclei(imageData, an, pathModel, channel, probThreshold, cellClass) {
    println '- Finding ' + channel + ' nuclei with Stardist -'
    def stardist = buildStarDistModel(pathModel, channel, probThreshold, cellClass)
    stardist.detectObjects(imageData, an, true)
    def cells = getDetectionObjects().findAll {
        it.getPathClass() == cellClass && an.getROI().contains(it.getROI().getCentroidX(), it.getROI().getCentroidY())
    }
    println 'Nb ' + channel + ' nuclei detected = ' + cells.size()
    stardist.close()
    return cells
}


// Detect cells in a specific annotation and channel with Cellpose
def detectCells(imageData, an, pathModel, channel, diameter, overlap, cellClass) {
    println '- Finding ' + channel + ' cells with Cellpose -'
    def cellpose = buildCellposeModel(pathModel, channel, diameter, overlap, cellClass)
    deselectAll()
    selectObjects(an)

    try {
        cellpose.detectObjects(imageData, getSelectedObjects())

        def cells = getDetectionObjects().findAll {
            it.getPathClass() == cellClass && an.getROI().contains(it.getROI().getCentroidX(), it.getROI().getCentroidY())
        }
        println 'Nb ' + channel + ' cells detected = ' + cells.size()
        return cells
    } catch(IOException) {
        Dialogs.showWarningNotification('Warning', 'Tiling error catched, Cellpose relaunched.')
        detectCells(imageData, an, pathModel, channel, diameter, overlap-1, cellClass)
    }
}


// Compute background in a specific annotation and channel: difference(annotation, union(detected cells in channel))
def getBackground(an, cells) {
    def geometryConverter = new GeometryTools.GeometryConverter.Builder().build()
    def union = geometryConverter.factory.buildGeometry(cells.stream().map(p -> p.getROI().getGeometry()).collect(Collectors.toList())).buffer(0)
    def geometry = an.getROI().getGeometry().difference(union)

    // Create the new ROI
    def bg = PathObjects.createAnnotationObject(GeometryTools.geometryToROI(geometry, an.getROI().getImagePlane()))
    an.addChildObject(bg)
    fireHierarchyUpdate()

    return bg
}

// Get background intensity (mean and SD) in a specific channel
def getBackgroundIntensity(bg, channel) {
    def mean = bg.getMeasurementList().getMeasurementValue('ROI: 2.00 µm per pixel: ' + channel + ': Mean')
    def SD = bg.getMeasurementList().getMeasurementValue('ROI: 2.00 µm per pixel: ' + channel + ': Std.dev.')
    println 'Background mean intensity in ' + channel + ' channel = ' + mean + ' (SD = ' + SD + ')'
    return [mean, SD]
}

// Filter cells by size (minCellSize < cell size < maxCellSize) and intensity (mean cell intensity > (C1 * bgMean + C2 * bgSD))
def filterCells(cells, channel, pixelWidth, minCellSize, maxCellSize, bgMean, bgSD, C1, C2) {
    def filteredCells = cells.findAll{
                it.getROI().getScaledArea(pixelWidth, pixelWidth) > minCellSize
                && it.getROI().getScaledArea(pixelWidth, pixelWidth) < maxCellSize
                && it.getMeasurementList().getMeasurementValue(channel + ': Mean') > (C1 * bgMean + C2 * bgSD)
    }
    println 'Nb ' + channel + ' cells remaining = ' + filteredCells.size() + ' (' + (cells.size() - filteredCells.size()) + ' filtered out by size and intensity)'
    return filteredCells
}

// Colocalize nuclei with Cy5 and Cy3 cells
def colocalization(nuclei, cellsCy5, cellsCy3) {
    println '- Colocalizing nuclei with Cy5 and Cy3 cells -'

    def tool = new RoiTools()
    def cells = []
    for(nucleus in nuclei) {
        def cell = new Cell(nucleus)
        def roiNuc = nucleus.getROI()

        for (cy5 in cellsCy5) {
            if (tool.areaContains(cy5.getROI(), roiNuc.getCentroidX(), roiNuc.getCentroidY())) {
                cy5.storeMetadataValue("Colocalized", "true")
                cell.setCy5(cy5)
                break
            }
        }

        for (cy3 in cellsCy3) {
            if (tool.areaContains(cy3.getROI(), roiNuc.getCentroidX(), roiNuc.getCentroidY())) {
                cy3.storeMetadataValue("Colocalized", "true")
                cell.setCy3(cy3)
                break
            }
        }

        cells << cell
    }
    return cells
}

// Class handling colocalization between nuclei and Cy5 and Cy3 cells
class Cell {
    def nucleus = null
    def isCy5 = false
    def cy5 = null
    def isCy3 = false
    def cy3 = null

    Cell(nucleus) {
        this.nucleus = nucleus
    }

    def getNucleusParams(dapiBg, cy5Bg, cy3Bg) {
        def area = this.nucleus.getMeasurementList().getMeasurementValue('Area µm^2').round(3)
        def circ = this.nucleus.getMeasurementList().getMeasurementValue('Circularity').round(3)
        def dapiIntMean = (this.nucleus.getMeasurementList().getMeasurementValue('DAPI: Mean') - dapiBg).round(3)
        def dapiIntSD = this.nucleus.getMeasurementList().getMeasurementValue('DAPI: Std.Dev.').round(3)
        def cy5IntMean = (this.nucleus.getMeasurementList().getMeasurementValue('Cy5: Mean') - cy5Bg).round(3)
        def cy5IntSD = this.nucleus.getMeasurementList().getMeasurementValue('Cy5: Std.Dev.').round(3)
        def cy3IntMean = (this.nucleus.getMeasurementList().getMeasurementValue('Cy3: Mean') - cy3Bg).round(3)
        def cy3IntSD = this.nucleus.getMeasurementList().getMeasurementValue('Cy3: Std.Dev.').round(3)
        return [area, circ, dapiIntMean, dapiIntSD, cy5IntMean, cy5IntSD, cy3IntMean, cy3IntSD]
    }

    def setCy5(cy5) {
        this.isCy5 = true
        this.cy5 = cy5
    }

    def getCy5Params(cy5Bg, cy3Bg) {
        if (this.isCy5) {
            def area = this.cy5.getMeasurementList().getMeasurementValue('Area µm^2').round(3)
            def cy5IntMean = (this.cy5.getMeasurementList().getMeasurementValue('Cy5: Mean') - cy5Bg).round(3)
            def cy5IntSD = this.cy5.getMeasurementList().getMeasurementValue('Cy5: Std.Dev.').round(3)
            def cy5cy3IntMean = (this.cy5.getMeasurementList().getMeasurementValue('Cy3: Mean') - cy3Bg).round(3)
            def cy5cy3IntSD = this.cy5.getMeasurementList().getMeasurementValue('Cy3: Std.Dev.').round(3)
            return [this.isCy5, area, cy5IntMean, cy5IntSD, cy5cy3IntMean, cy5cy3IntSD]
        } else {
            return [this.isCy5, null, null, null, null, null]
        }
    }

    def setCy3(cy3) {
        this.isCy3 = true
        this.cy3 = cy3
    }

    def getCy3Params(cy3Bg, cy5Bg) {
        if (this.isCy3) {
            def area = this.cy3.getMeasurementList().getMeasurementValue('Area µm^2').round(3)
            def cy3IntMean = (this.cy3.getMeasurementList().getMeasurementValue('Cy3: Mean') - cy3Bg).round(3)
            def cy3IntSD = this.cy3.getMeasurementList().getMeasurementValue('Cy3: Std.Dev.').round(3)
            def cy3Cy5IntMean = (this.cy3.getMeasurementList().getMeasurementValue('Cy5: Mean') - cy5Bg).round(3)
            def cy3cy5IntSD = this.cy3.getMeasurementList().getMeasurementValue('Cy5: Std.Dev.').round(3)
            return [this.isCy3, area, cy3IntMean, cy3IntSD, cy3Cy5IntMean, cy3cy5IntSD]
        } else {
            return [this.isCy3, null, null, null, null, null]
        }
    }
}

// Count nuclei being Cy5-negative/positive and Cy3-negative/positive
def countCells(cells) {
    def cy5 = 0, cy3 = 0, noneNone = 0, cy5None = 0, noneCy3 = 0, cy5Cy3 = 0
    for(cell in cells) {
        if(cell.isCy5) cy5++
        if(cell.isCy3) cy3++
        if(!cell.isCy3 && !cell.isCy5) noneNone++
        if(cell.isCy5 && !cell.isCy3) cy5None++
        if(!cell.isCy5 && cell.isCy3) noneCy3++
        if(cell.isCy5 && cell.isCy3) cy5Cy3++
    }
    return [cy5, cy3, noneNone, cy5None, noneCy3, cy5Cy3]
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
