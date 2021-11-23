#########################################################################################################################################################
####################################################### PUSH "RUN" TO START THE SCRIPT ##################################################################
#########################################################################################################################################################

############################################################ RELEASE NOTES ##############################################################################

# Evaluation of both .czi and .zvi 12-bit images with 2 channels (DNA and foci)
# Automatic creation of Results and Analyzed_cells folders
# In the .csv file information about: 'Cell No', 'DAPI min', 'DAPI max', 'DAPI mean', 'DAPI homogeneity', 'GFP min', 'GFP max', 'GFP mean', 
#'GFP homogeneity', '% Area above threshold GFP', 'Foci size', 'Cell Area', 'N_Foci'
# ROI manager automatically saved in the folder where the images are - so that you can use them after for reevaluation

# Microscope calibration - the user can add the parameters from the microscope that he used
# The user can select if he wants to blur the foci channel in case that it is noisy
# Posibility of IN VIVO evaluations - vessel selection
# Can select if he wants to reevaluate cells if ROI manager saved
# Manual selection of noise tolerance 
# Manual or automatic selection of threshold
# User can estimate the cell sizes that he wants to include to the evaluation
# User can use exlusion criteria to exclude overexposed cells or include cells that are zero, but in case of normal evaluation would be considered as 1.

#########################################################################################################################################################
#########################################IMPORT OF ALL THE USED LIBRARIES AND SELECTION OF PARAMETERS ###################################################
#########################################################################################################################################################

from ij import IJ
from ij import ImageStack
from ij import ImagePlus
from ij.plugin.frame import RoiManager
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from ij.gui import WaitForUserDialog
from ij.plugin.filter import MaximumFinder
from ij.gui import GenericDialog  
from java.awt import Font
import csv, os, glob

# If set to true, windows with selected and analysed cells are left opened
verify_params = False

# If background set to a different colour than Black, then the script will not work correctly
IJ.run("Colors...", "foreground=magenta background=black selection=yellow")

#########################################################################################################################################################
#################################################### DEFINITION OF FUNCTIONS TO BE USED IN THE MAIN #####################################################
#########################################################################################################################################################

def opensavedROIman(ROIopenpath):
	#opens a saved ROImanager with the set of selected ROIs	
	rm = RoiManager.getInstance()
	if not rm:
  		rm = RoiManager()
	rm.reset()
	rm.runCommand("Open", ROIopenpath)	
	return rm

def ThresholdEst(DAPI, foci, thres_a, thres_b):
	#estimates the threshold for the MaximumFinder using the background without the cells which are excluded using a mask created from the DAPI channel
	DAPI_seg = DAPI.duplicate()
	DAPI_seg.show
	IJ.run(DAPI_seg, "8-bit", "")
	# run cell segmentation on DAPI channel. Method: Default
	IJ.run(DAPI_seg, "Subtract Background...", "rolling=100")
	IJ.run(DAPI_seg, "Gaussian Blur...", "sigma=5")
	IJ.run(DAPI_seg, "Enhance Contrast...", "saturated=10 normalize equalize")
	IJ.setAutoThreshold(DAPI, "Default dark no-reset")
	IJ.run(DAPI_seg, "Convert to Mask", "")
	IJ.run(DAPI_seg, "Fill Holes", "")
	IJ.run(DAPI_seg, "Create Selection", "")		
	roi1 = DAPI_seg.getRoi()	
	foci_img = foci.duplicate()
	foci_img.setRoi(roi1)
	IJ.run(foci_img, "Clear", "slice")
	IJ.run(foci_img, "Make Inverse", "")
	IJ.run(foci_img, "Measure", "") # Without this line it does not work correctly. Keep it!!!
	stat = foci_img.getProcessor().getStatistics()
	threshold = round(thres_a * stat.mean + thres_b) # Function obtained from threshold analysis 
	foci_img.hide()
	DAPI_seg.hide()	
	return threshold
	
def Select_Area(imp): # for the vessel selection in the in vivo samples from the DAPI channel
	img = imp.duplicate()	
	rm = RoiManager.getInstance()
	if not rm:
		rm = RoiManager()
	rm.runCommand("reset")	
	img.show()
	#IJ.setTool(Toolbar.FREEROI) # Maybe you will have to comment this line - it depends on the Java version you have on your PC
	WaitForUserDialog("Action required", "Select the vessel using the freehand selection tool."+
						" \nPlease press OK when done.").show()					
	roi1 = img.getRoi()
	img.close()			
	return roi1	
	
def Cell_Segmentation(imp):
	DAPI = imp.duplicate()
	IJ.run(DAPI, "8-bit", "")
	IJ.run(DAPI, "Subtract Background...", "rolling=50")
	IJ.run(DAPI, "Gaussian Blur...", "sigma=1")
	IJ.run(DAPI, "Enhance Contrast...", "saturated=10 normalize equalize")
	IJ.setAutoThreshold(DAPI, "Default dark no-reset")
	IJ.run(DAPI, "Convert to Mask", "")
	IJ.run(DAPI, "Fill Holes", "")
	IJ.run(DAPI, "Watershed", "")
	IJ.run(DAPI, "Create Selection", "")	
	IJ.run(DAPI, "Select None", "")
	DAPI.setTitle("Cell_Segmentation")
	return DAPI 

def Pick_Cells(Mask, DAPI, RoiManagerInstance, ROIsavepath): 
	# Pick the right cells ex vivo	
	DAPI_img = DAPI.duplicate()
	DAPI_img.show()
	IJ.run(Mask, "Create Selection", "")
	IJ.run(Mask, "RGB Color", "")
	Mask.show()	
	IJ.run("Add Image...", "image=" + DAPI_img.getTitle() + " x=0 y=0 opacity=50")
	WaitForUserDialog("Action required", "Use floodfill tool to pick cells with good segmentation"+
						" \nPlease press OK when done.").show()
	IJ.run(Mask, "Select None", "")
	IJ.run(Mask, "Remove Overlay", "")
	IJ.run(Mask, "8-bit", "")
	IJ.setThreshold(Mask, 32, 233)
	IJ.run(Mask, "Create Selection", "")
	IJ.run(Mask, "Clear Outside", "")
	IJ.run(Mask, "Set...", "value=255")
	RoiManagerInstance.addRoi(Mask.getRoi())
	RoiManagerInstance.runCommand("Split")
	RoiManagerInstance.select(0)
	RoiManagerInstance.runCommand(Mask,"Delete")
	RoiManagerInstance.runCommand("Save", ROIsavepath + ".zip")
	DAPI_img.close()
	Mask.hide()
	return Mask, RoiManagerInstance
	
def Pick_Cells_InVivo(Mask, DAPI, RoiManagerInstance, ROIsavepath, roi1): 
	# Pick the right cells in vivo
	# First create the area for the selection of the cells on the original DAPI image. 
	# Two lines are draw on the image. Between the lines we have 45 um.
	DAPI_img = DAPI.duplicate()	
	DAPI_img.show()	
	DAPI_img.setRoi(roi1)
	IJ.run(DAPI_img, "RGB Color", "")
	IJ.run("Colors...", "foreground=yellow background=black selection=yellow")
	IJ.run("Draw", "slice")
	IJ.run("Enlarge...", "enlarge=15")
	IJ.run("Enlarge...", "enlarge=15")
	IJ.run("Enlarge...", "enlarge=15")
	IJ.run("Draw", "slice")		
	Mask.show()
	IJ.run(Mask, "Create Selection", "")
	IJ.run(Mask, "RGB Color", "")			
	IJ.run("Add Image...", "image=" + DAPI_img.getTitle() + " x=0 y=0 opacity=50")	
	IJ.run("Colors...", "foreground=magenta background=black selection=yellow")
	WaitForUserDialog("Action required", "Use floodfill tool to pick cells with good segmentation"+
						" \nPlease press OK when done.").show()
	IJ.run(Mask, "Select None", "")
	IJ.run(Mask, "Remove Overlay", "")
	IJ.run(Mask, "8-bit", "")
	IJ.setThreshold(Mask, 32, 233)
	IJ.run(Mask, "Create Selection", "")
	IJ.run(Mask, "Clear Outside", "")
	IJ.run(Mask, "Set...", "value=255")			
	RoiManagerInstance.runCommand("reset")
	RoiManagerInstance.addRoi(Mask.getRoi())
	RoiManagerInstance.runCommand("Split")
	RoiManagerInstance.select(0)
	RoiManagerInstance.runCommand(Mask,"Delete")
	RoiManagerInstance.runCommand("Save", ROIsavepath + ".zip")
	DAPI_img.hide()
	Mask.hide()	
	return Mask, RoiManagerInstance
	
def excludeDAPI(Stack, cell_ROI, filename1):
	# This fuction gives information about the DAPI channel
	# The DAPI and GFP channel are in the opposite direction in the .czi and .zvi image formats 
	if filename1.endswith('.czi'):
		Stack.setSlice(2)
	else:
		Stack.setSlice(1)
	img = Stack.crop()
	ip = img.getProcessor()
	ip.setRoi(cell_ROI)
	stats = ip.getStatistics()
	DAminV = round(stats.min,2)
	DAmaxV = round(stats.max,2)
	DAmeanV = round(stats.mean,2)
	DAhomogen = round(DAmaxV/DAmeanV,3) 
	return DAmeanV, DAminV, DAmaxV, DAhomogen

def findFoci(Stack, cell_ROI, threshold, filename1):
	# Uses Find Maxima to detect foci in cell
	mf = MaximumFinder()
	# The DAPI and GFP channel are in the opposite direction in the .czi and .zvi image formats 
	if filename1.endswith('.czi'):
		Stack.setSlice(1)
	else:
		Stack.setSlice(2)
	img = Stack.crop()
	if Gaussian_blur_use == True:
		IJ.run(img, "Gaussian Blur...", "sigma="+sigma_foci+"")
	ip = img.getProcessor()
	ip2 = img.getProcessor().convertToFloat()
	ip.setRoi(cell_ROI)
	ip2.setRoi(cell_ROI)
	pixels = ip2.getPixels()
	stats = ip.getStatistics()
	minV = round(stats.min,2)
	maxV = round(stats.max,2)
	meanV = round(stats.mean,2)
	PecentAreaThres = round(reduce(lambda count, a: count + 1 if a > threshold else count, pixels, 0) / float(len(pixels)) * 100,2)
	GFPhomogen = round(maxV/meanV,3)
	maxima = mf.findMaxima(ip, noise_tolerance, threshold, MaximumFinder.SINGLE_POINTS, False, False)
	maxima2 =  mf.findMaxima(ip, noise_tolerance, threshold, MaximumFinder.SINGLE_POINTS, False, False)
	maxima2.invert()
	findmaxinvert = ImagePlus("Found Maxima", maxima2)
	findmaximashow = ImagePlus("Maxima", maxima)
	maximaip = findmaximashow.getProcessor()
	maximahist = maximaip.getHistogram()
	CountMethod = maximahist[255]
	Points = str(CountMethod)								
	return Points, findmaxinvert, meanV, minV, maxV, PecentAreaThres, GFPhomogen	

def Postprocess(Stack, cell_ROI):
	# Does a bit of postprocessing: adjust windowing for every channel and sets the foci as active ROI selection
	Stack.setRoi(cell_ROI);
	IJ.run(Stack, "Enhance Contrast...", "saturated=0.3 normalize process_all")
	IJ.run(Stack, "8-bit", "")		

def Cell_Area(IMPStack, filename1):
	# Gets Cell Area
	# Pick DAPI channel with black background (pixel value = 0)
	# The DAPI and GFP channel are in the opposite direction in the .czi and .zvi image formats 
	if filename1.endswith('.czi'):
		IMPStack.setSlice(1)
	else:
		IMPStack.setSlice(0)	
	# Invert to get correct selection
	cell_mask = IMPStack.crop()
	IJ.run(cell_mask, "Enhance Contrast", "saturated=0.35")
	IJ.run(cell_mask, "8-bit", "slice")
	cell_mask.show()
	IJ.setThreshold(1, 255)
	IJ.run(cell_mask, "Convert to Mask", "background=Dark calculate only black")
	IJ.run(cell_mask, "Create Selection", "")
	Area = round(cell_mask.getStatistics().area,2)
	ROI = cell_mask.getRoi()
	cell_mask.hide()
	return cell_mask, Area, ROI
	
##############################################################################################################################################
######################################################## MAIN FUNCTION #######################################################################
##############################################################################################################################################

def main():
	#################################### Welcome message and choise of parameters of analysis ################################################
	# Definition of global variables
	global threshold
	global noise_tolerance
	global TresSelect
	global RoiManHave
	global sigma_foci
	global Gaussian_blur_use
	
	#First script's dialog about microscope and pixel calibration
	font = Font("Arial", Font.BOLD, 14)
	font2 = Font("Arial", Font.BOLD, 12)
	gd = GenericDialog("Welcome to our Foci script! Insert the basic information.")
	gd.addMessage("This script evaluates foci from all the images (.czi or .zvi) with two chanels (DNA and foci) in a directory.",font)
	gd.addMessage("Firstly you have to define the pixel calibration. Each microscope/objective have a different micrometer value corresponding to a pixel." +
	"\nYou can find this value in the information of your image. Open it in ZEN and find this value in the information of the image. "+
	"\nHere the predefined value (0.161) correspond to AxioImager M1 with 40x objective.",font)
	gd.setInsets(10,40,10)
	gd.addNumericField("How many micrometers is 1 pixel in your case?", 0.161, 3)
	gd.addMessage("If your images are in their foci channel noisy then you can blur them if you want." +
	"\nIf you want to blur them then check the next field and add the sigma of the Gaussian blur that you want to use" ,font)
	gd.setInsets(5,40,5)
	gd.addCheckbox("I want to use the Gaussian blur in the foci channel", False)
	gd.setInsets(10,40,10)
	gd.addNumericField("What sigma you want to use:",1,1)
	gd.addMessage("If you want to select the vessels and then pick only cells 0.45 um from them, then check the next box.",font)
	gd.setInsets(5,40,5)
	gd.addCheckbox("I am analysing IN VIVO samples and I want to select the vessel on the image.", False)	
	gd.addMessage("If you already have the ROI manager then check the next box.",font)		
	gd.setInsets(5,40,5)
	gd.addCheckbox("I already have the ROI manager saved and I want to analyse these cells again.", False)	
	gd.setCancelLabel("Exit")
	gd.setOKLabel("Next step")
	gd.showDialog()
	if gd.wasCanceled(): 
		return 0   
	elif gd.wasOKed():
		pixel_size = gd.getNextNumber()
		pixel_cal = 1/pixel_size
		pixel_cal = str(pixel_cal)
		Gaussian_blur_use = gd.getNextBoolean()
		sigma_foci = gd.getNextNumber()
		sigma_foci = str(sigma_foci)
		TumorAnalysis = gd.getNextBoolean()	
		RoiManHave = gd.getNextBoolean()
		if (RoiManHave == True and TumorAnalysis == True):
			TumorAnalysis == False
	
	#Definition of noise tolerance and threshold
	gd = GenericDialog("Noise toleracne and threshold")
	gd.addMessage("Set the required noise tolerance - this value affects the detections of near by foci.",font)
	gd.setInsets(10,40,10)
	gd.addNumericField("Insert the noise tolerance", 100, 1)
	gd.addMessage("You can select either choosing your own threshold or it can be estimated from each image in a folder automatically." + 
	"\nIf you want the automatic threshold then do not check the next field.",font)
	gd.setInsets(5,40,5)
	gd.addCheckbox("I want to set the threshold manually using the numerical field below. ", False)
	gd.setInsets(10,40,10)
	gd.addNumericField("Manual Threshold", 500, 1)
	gd.addMessage("I want the automatic threshold estimation, so I am entering the parameters of my linear function." +
	"\nYou need to create a function: 'threshold = a * MEAN + b'. See the manual for further explanations." +
	"\nThe parameters given here are from a former calibration and they will probably not fit to your situation.",font)
	gd.setInsets(10,40,10)
	gd.addNumericField("Insert the first value for the threshold calibration a=",0.65,2)
	gd.setInsets(10,40,10)
	gd.addNumericField("Insert the second value for the threshold calibration b=",350,1)
	gd.setCancelLabel("Exit")
	gd.setOKLabel("Next step")
	gd.showDialog()
	if gd.wasCanceled(): 
		return 0   
	elif gd.wasOKed():
		noise_toler = gd.getNextNumber()
		TresSelect = gd.getNextBoolean()		
		manualTres = gd.getNextNumber()	
		thres_a = gd.getNextNumber() # The "a" parameter of the used calibration threshold = a * MEAN + b.
		thres_b = gd.getNextNumber() # The "b" parameter of the used calibration threshold = a * MEAN + b.		
				
	
	#Definition of cell parameters		
	gd = GenericDialog("Cell sizes to be considered")
	gd.addMessage("Here you can define the sizes of the cells you want to analyse."+
	"\nIn case that you will pick a cell than is outside these limits, then it will be excluded from the evaluation",font)
	gd.setInsets(10,40,10)
	gd.addNumericField("Smallest cell area you want to detect (in squared um)", 15, 1)
	gd.setInsets(10,40,10)
	gd.addNumericField("Largest cell area you want to detect (in squared um)", 300, 1)	
	gd.setCancelLabel("Exit")
	gd.setOKLabel("Next step")
	gd.showDialog()
	if gd.wasCanceled(): 
		return 0   
	elif gd.wasOKed():
		Min_cellarea = gd.getNextNumber()
		Max_cellarea = gd.getNextNumber()

	#Exclution criteria
	gd = GenericDialog("Definition of exclusion criteria")
	gd.addMessage("In case that you want to use exclution criteria, then you should firstly try this script without any and estimate your own criteria."+
	"\nOf course you can use the criteria given here, but it is not sure that they will fit to your data.",font)
	gd.addCheckbox("I want to use exclution criteria", False)
	gd.addMessage("Definition of zero cells",font2)
	gd.setInsets(10,40,10)
	gd.addNumericField("How many times higher the maximum value should not be than the threshold while the cell homogeneity is low to be considered as a zero cell",3,1)
	gd.setInsets(10,40,10)
	gd.addNumericField("To what value of homogeneity you want to consider the cell as zero if there is not any higher signal than that you defined in the field above",2,1)
	gd.addMessage("Definition of overexposed cells",font2)
	gd.setInsets(10,40,10)
	gd.addNumericField("How many times the average cell intensity value should be higher than the threshold to be considered as overexposed",4,1)
	gd.setInsets(10,40,10)
	gd.addNumericField("How much area of the cell above the threshold in % should be filled to be considered as overexposed",75,1)
	gd.setInsets(10,40,10)
	gd.addNumericField("From how big foci area the cell should be considered as strange and so overexposed",5,1)
	gd.addMessage("Sometimes you could want to calculate cells that where excluded, but you want to divide the area that was detected by a number to get foci number.", font)
	gd.setInsets(5,40,5)
	gd.addCheckbox("I want to calculate some cells that where due to the criteria above excluded using division",False)
	gd.setInsets(10,40,10)
	gd.addNumericField("From what homogeneity value the cell should be calculated by using division by the average foci size",2.3,1)
	gd.setInsets(10,40,10)
	gd.addNumericField("What is the average foci size you want to use for the division",3.5,1)
	gd.setCancelLabel("Exit")
	gd.setOKLabel("Start calculating!")
	gd.showDialog()
	if gd.wasCanceled(): 
		return 0   
	elif gd.wasOKed():
		ExclutionCriteria = gd.getNextBoolean()
		zeroMax = gd.getNextNumber()	
		zeroHom = gd.getNextNumber()
		ExcluMean = gd.getNextNumber()
		AreaAboveThres = gd.getNextNumber()
		MaxFS = gd.getNextNumber()
		DivideOEcells = gd.getNextBoolean()
		HomoToDiv = gd.getNextNumber()
		UserEstFociSize = gd.getNextNumber()
		 
	############################## Pick an image from a directory that you want to analyze ##################################################
	
	filename1 = IJ.getFilePath('Pick an image (.czi or .zvi) from the directory you want to analyze')	
	if not filename1:
		gd = GenericDialog("Error") 
		gd.addMessage("An image was not picked. Run the script again and pick an image.")	
		gd.hideCancelButton()	
		gd.showDialog()
		return 0
	if not (filename1.endswith('.czi') or filename1.endswith('.zvi')): 
		gd = GenericDialog("Error") 
		gd.addMessage("File is not .czi neither .zvi. Run the script again and pick an image .czi or .zvi.")
		gd.hideCancelButton()		
		gd.showDialog()
		return 0
		
	############################################ Definition of direcories ##################################################################
	
	# Τhe directory with images acquired with the same parameters	
	AnalysisDir = os.path.dirname(filename1) 
	# Creating directory to collect all the CSV files		
	resultpath2 = os.path.join(os.path.dirname(filename1), "Results") 
	resultpath2 = os.path.normpath(resultpath2)	
	if not os.path.isdir(resultpath2):
		try:
			os.mkdir(resultpath2)
		except:
			return 0
	# Creating directory for saving the individual imeges of the cells	
	Imagespath = os.path.join(os.path.dirname(filename1), "Analysed_cells") 
	Imagespath = os.path.normpath(Imagespath)	
	if not os.path.isdir(Imagespath):
		try:
			os.mkdir(Imagespath)
		except:
			return 0
	
	####################### Threshold establishment in case of manual and general threshold for one image in the directory #############
	
	if TresSelect == True:
		threshold = manualTres
			
    ##################################################### MAIN ANALYSIS ###############################################################
	
	# Definition of the directory with the images to be analysed
	if filename1.endswith('.czi'):
		workDir = glob.glob(os.path.join(AnalysisDir,'*.czi'))
		
	else:
		workDir = glob.glob(os.path.join(AnalysisDir,'*.zvi')) 
	
	NoImages = len(workDir) # number of images in the directory
	print 'There are ' + str(NoImages) + ' images for analysis in this folder.'
			
	for filename in workDir: # Analysis of all the images in the chosen directory	
		try:
			rm = RoiManager.getInstance()
			rm.reset()
		except:
			rm = RoiManager()
		try: 
			IJ.run("Close All", "")
		except:
			pass			
	
		# Open the image
		options = ImporterOptions()
		options.setId(filename)
		options.setSplitChannels(True)		
		imps = BF.openImagePlus(options)	
		if filename1.endswith('.czi'):	
			imps[1].setTitle("DAPI")
			DAPI = imps[1]
			imps[0].setTitle("GFP")
			GFP = imps[0]
		else: 
			imps[0].setTitle("DAPI")
			DAPI = imps[0]
			imps[1].setTitle("GFP")
			GFP = imps[1]

		# In case that the acquisition of the images was done with a low exposure, 
		# then I will try to change to noisetolerance accordingly.
		helpStat = GFP.getProcessor().getStatistics() 
		helpNoise =  4095/helpStat.max
		if helpStat.max < 4095:
			noise_tolerance = round(float(noise_toler)/float(helpNoise))
		else:
			noise_tolerance = noise_toler	
			
		# Setting the threshold from each image automatically			
		if TresSelect == False:	
			threshold = ThresholdEst(DAPI,GFP,thres_a,thres_b)			

		# Creation of folder with the analysed cells					
		resultpath = os.path.join(Imagespath, os.path.basename(filename)[:-4] + "_analyzed_cells_MaxFind_noiseTol_" + str(noise_tolerance)+ '_Threshold_' + str(threshold))
		resultpath = os.path.normpath(resultpath)	
		if not os.path.isdir(resultpath):
			try:
				os.mkdir(resultpath)
			except:
				return 0
				
		if RoiManHave == True: 
			resultpathROI = os.path.join(os.path.dirname(filename), os.path.basename(filename)[:-4] +"_ROI.zip")
			rm = opensavedROIman(resultpathROI)
		elif RoiManHave == False:
			resultpathROI = os.path.join(os.path.dirname(filename), os.path.basename(filename)[:-4] +"_ROI")				
		
			# Process DAPI channel - make the segmentation and select the cell to be analysed
			Cell_Map = Cell_Segmentation(DAPI)			
			if TumorAnalysis == True: 
				roi1 = Select_Area(DAPI)
				Cell_Map, rm = Pick_Cells_InVivo(Cell_Map, DAPI, rm, resultpathROI,roi1)
			else:	
				Cell_Map, rm = Pick_Cells(Cell_Map, DAPI, rm, resultpathROI)				

		######################################## Creation of the csv textbook for the results #################################################
		if RoiManHave == True:
			nameCSV = resultpath2 + '\\' + os.path.basename(filename)[:-4]+'_results_MaxFind_ROI_noiseTol_' + str(noise_tolerance) + '_UserThreshold_' + str(threshold) + '.csv'
		else:
			nameCSV = resultpath2 + '\\' + os.path.basename(filename)[:-4]+'_results_MaxFind_ORIG_noiseTol_' + str(noise_tolerance) + '_UserThreshold_' + str(threshold) + '.csv'
			
		with open(nameCSV, 'wb') as csvfile:
			spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			spamwriter.writerow(['Cell No', 'DAPI min', 'DAPI max', 'DAPI mean', 'DAPI homogeneity', 'GFP min', 'GFP max', 'GFP mean', 'GFP homogeneity', '% Area above threshold GFP', 'Foci size', 'Cell Area', 'N_Foci'])				

		######################################	Analysis of all selected cells ##############################################################
			Data = []
			N = rm.getCount()	
			ROIs = rm.getRoisAsArray()
			rm.reset()
			
			for ROI in ROIs: 	
				rm.addRoi(ROI)
				# copy cropped cell to separate image (all channels)
				Stack = ImageStack(int(ROI.getFloatWidth()),int(ROI.getFloatHeight()))
				
				# browse all channels
				for imp in imps:
					# duplicate image, clear area around ROI and add the remains to the Stack
					# Why? If cells are close to each other, a rectangular image celection around one cell 
					# might include areas of other cells and we do not want that
					imp_buffer = imp.duplicate()
					imp_buffer.setRoi(ROI)
					IJ.setAutoThreshold(imp_buffer, "Default dark reset")
					IJ.run(imp_buffer, "Clear Outside", "") 
					img = imp_buffer.crop()
					imp_buffer.close()	
					# add channels to Stack
					ip = img.getProcessor()
					Stack.addSlice(img.getTitle(), ip)
					
				# Turn Stack into ImagePlus
				Stack = ImagePlus(ROI.getName(), Stack)
				IJ.run(Stack, "Set Scale...", "distance="+pixel_cal+" known=1 pixel=1 unit=micron global")
	
				# Get Cell Area and ROI of Cell
				Stack.show()
				Cell_mask, Area_Cell, Cell_ROI = Cell_Area(Stack, filename1)

				# Get info from DAPI channell
				DAmeanV, DAminV, DAmaxV, DAhomogen = excludeDAPI(Stack, Cell_ROI, filename1)
				
				# Get Maxima/foci and several parameters of interest
				Points, findmaxinvert, meanV, minV, maxV, PecentAreaThres, GFPhomogen = findFoci(Stack, Cell_ROI, threshold, filename1)
				
				# Mean foci size
				if (int(Points) == 0):
					FS = "NA"
				else:	
					FS = round(float(Area_Cell) * float(PecentAreaThres) / 100 / int(Points),3)

				############################################### EXCLUSION CRITERIA ##############################################################
				if ExclutionCriteria == True:
					if float(Area_Cell) < Min_cellarea or float(Area_Cell) > Max_cellarea:
						Points = "OE"
						FS = 0
						FA = "NA"
					elif PecentAreaThres == 0:
						Points = 0
						FS = 0
						FA = 0	
					elif (float(maxV) < zeroMax * threshold and float(GFPhomogen) < zeroHom):
						Points = 0
						FS = 0
						FA = 0	
					elif ((float(meanV) > (ExcluMean * threshold)) or (float(PecentAreaThres) > AreaAboveThres) or (float(FS) >= MaxFS)):
						if DivideOEcells == True:
							if (float(GFPhomogen) > HomoToDiv): 
								help = round(float(Area_Cell)*float(PecentAreaThres)/100/UserEstFociSize)
								if int(Points) < float(help):
									Points = help
								else:
									Points = Points		
						else:
							Points = "OE"
							FA = "NA"
					elif Points == 0:
						FA = 0	
						FS = 0

				########################################## FINAL PROCESSING OF THE RESULTS ######################################################
						
				#Postprocess
				Postprocess(Stack, Cell_ROI)
	
				#Add another slice to the stack with the points that MaximumFinder took into account			
				Stack.getStack().addSlice('Found Maxima', findmaxinvert.getProcessor())			
											
				# Save image of cell
				Stack.show()
				IJ.saveAs("Tiff", resultpath + "\\" + Stack.getTitle())
							
				# If parameter verification is disabled: Close cell image.
				if verify_params:
					Stack.setSlice(5)
					IJ.run(Stack, "Options...", "black")
				else:
					Stack.close()
				# Write to results file
				spamwriter.writerow([Stack.getTitle()[:-4], DAminV, DAmaxV, DAmeanV, DAhomogen, minV, maxV, meanV, GFPhomogen, PecentAreaThres, FS, Area_Cell, Points])
			rm.runCommand("reset")
		NoImages = NoImages - 1
		print 'The are still ' + str(NoImages) + ' images for analysis in the folder.'					
	gd = GenericDialog("Progress") 
	gd.addMessage("All the images from the selected folder have been evaluated. Good job!")
	gd.showDialog()

main()


