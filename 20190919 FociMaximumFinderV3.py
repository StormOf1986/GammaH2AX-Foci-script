#########################################################################################################################################################
####################################################### PUSH "RUN" TO START THE SCRIPT ##################################################################
#########################################################################################################################################################

############################################################ RELEASE NOTES ##############################################################################

# Evaluation of both .czi and .zvi
# Manual selection of noise tolerance at the beggining
# Manual or from one image from a folder or from each image in a folder selection of threshold
# ROI manager saved in the folder where the images are - so that you can use them after for reevaluation
# Can reevaluate cells if ROI manager saved
# Automatic creation of Results and Analyzed_cells folders
# Overexposed cells are those that have a mean intensity of more than 3500 in the GFP channel
# In the .csv file information about: 'Cell No', 'Cell Area', '% Area above threshold', 'Min intensity', 'Max intensity', 'Mean intensity', 'N_Foci', 'Foci/Area'
# In 20190919 FociMaximumFinderV2 problem with FA when N_Foci was OE. Solved in this version.

#########################################################################################################################################################
#########################################IMPORT OF ALL THE USED LIBRARIES AND SELECTION OF PARAMETERS ###################################################
#########################################################################################################################################################

from ij import IJ
from ij import ImageStack
from ij import ImagePlus
from ij.process import FloatProcessor
from ij.plugin.frame import RoiManager
from loci.plugins import BF
from ij.plugin import Duplicator
from loci.plugins.in import ImporterOptions
from ij.gui import WaitForUserDialog
from ij.plugin.filter import MaximumFinder
from ij.plugin import RGBStackMerge
from jarray import array
from ij.gui import GenericDialog  
import csv, os, glob

# If set to true, windows with selected and analysed cells are left opened
verify_params = False; 

# If background set to a different colour than Black, then the script will not work correctly
IJ.run("Colors...", "foreground=magenta background=black selection=yellow"); 

#########################################################################################################################################################
#################################################### DEFINITION OF FUNCTIONS TO BE USED IN THE MAIN #####################################################
#########################################################################################################################################################

def opensavedROIman(ROIopenpath):
	#opens a saved ROImanager with the set of selected ROIs	
	rm = RoiManager.getInstance();
	if not rm:
  		rm = RoiManager()
	rm.reset();
	rm.runCommand("Open", ROIopenpath);
	
	return rm

def ThresholdEst(image):
	#estimates the threshold for the MaximumFinder using a user defined rectangular ROI
	#just select an area where is background
	imp = image.duplicate();	
	rm = RoiManager.getInstance()
	if not rm:
		rm = RoiManager()
	rm.runCommand("reset")	
	imp.show();	
	IJ.setTool(Toolbar.RECTANGLE);
	WaitForUserDialog("Action required", "Select an area of backgroud (cells without foci) to establish the threshold"+
						" \nPlease press OK when done.").show();					
	roi1 = imp.getRoi();
	imp.setRoi(roi1);
	stats = imp.getProcessor().getStatistics();
	threshold = round(stats.max)+100;
	imp.close();	

	return threshold

def Cell_Segmentation(imp):
	# duplicate input image
	DAPI = imp.duplicate();
	IJ.run(DAPI, "8-bit", "");

	# run cell segmentation on DAPI channel. Method: Default
	IJ.run(DAPI, "Subtract Background...", "rolling=50");
	IJ.run(DAPI, "Gaussian Blur...", "sigma=1");
	IJ.run(DAPI, "Enhance Contrast...", "saturated=10 normalize equalize");
	IJ.setAutoThreshold(DAPI, "Default dark no-reset");
	IJ.run(DAPI, "Convert to Mask", "");
	IJ.run(DAPI, "Fill Holes", "");
	IJ.run(DAPI, "Watershed", "");
	IJ.run(DAPI, "Create Selection", "");
	ROIs = DAPI.getRoi();
	IJ.run(DAPI, "Select None", "");
	DAPI.setTitle("Cell_Segmentation");	

	return DAPI , ROIs

def Pick_Cells(Mask, DAPI, RoiManagerInstance, ROIsavepath):
	# Pick the right cells	
	DAPI_img = DAPI.duplicate();
	DAPI_img.show();

	IJ.run(Mask, "Create Selection", "");
	IJ.run(Mask, "RGB Color", "");
	Mask.show();
	
	IJ.run("Add Image...", "image=" + DAPI_img.getTitle() + " x=0 y=0 opacity=10");
	WaitForUserDialog("Action required", "Use floodfill tool to pick cells with good segmentation"+
						" \nPlease press OK when done.").show();

	IJ.run(Mask, "Select None", "");
	IJ.run(Mask, "Remove Overlay", "");
	IJ.run(Mask, "8-bit", "");
	IJ.setThreshold(Mask, 32, 233);
	IJ.run(Mask, "Create Selection", "");
	IJ.run(Mask, "Clear Outside", "");
	IJ.run(Mask, "Set...", "value=255");
	
	# Add to Manager
	RoiManagerInstance.addRoi(Mask.getRoi());
	RoiManagerInstance.runCommand("Split");
	RoiManagerInstance.select(0);
	RoiManagerInstance.runCommand(Mask,"Delete");
	RoiManagerInstance.runCommand("Save", ROIsavepath + ".zip");
	DAPI_img.close();
	Mask.hide();

	return Mask, RoiManagerInstance

def findFoci(Stack, cell_ROI, threshold, filename1):
	# Uses Find Maxima to detect foci in cell
	mf = MaximumFinder();
	if filename1.endswith('.czi'):
		Stack.setSlice(1);
	else:
		Stack.setSlice(2);
	img = Stack.crop();
	ip = img.getProcessor();
	ip2 = img.getProcessor().convertToFloat();
	ip.setRoi(cell_ROI);
	ip2.setRoi(cell_ROI);	
	pixels = ip2.getPixels();
	stats = ip.getStatistics();
	minV = round(stats.min,2);
	maxV = round(stats.max,2);
	meanV = round(stats.mean,2);
	PecentAreaThres = round(reduce(lambda count, a: count + 1 if a > threshold else count, pixels, 0) / float(len(pixels)) * 100,2);
	if meanV >= 3500:
		Points = "OE"
		maxima2 =  mf.findMaxima(ip, noise_tolerance, threshold, MaximumFinder.SINGLE_POINTS, False, False);
		maxima2.invert();
		findmaxinvert = ImagePlus("Found Maxima", maxima2);
	else:
		maxima = mf.findMaxima(ip, noise_tolerance, threshold, MaximumFinder.SINGLE_POINTS, False, False);
		maxima2 =  mf.findMaxima(ip, noise_tolerance, threshold, MaximumFinder.SINGLE_POINTS, False, False);
		maxima2.invert();
		findmaxinvert = ImagePlus("Found Maxima", maxima2);
		findmaximashow = ImagePlus("Maxima", maxima);
		maximaip = findmaximashow.getProcessor();
		maximahist = maximaip.getHistogram();
		CountMethod = maximahist[255];
		Points = str(CountMethod);
		
	return Points, findmaxinvert, meanV, minV, maxV, PecentAreaThres	

def Postprocess(Stack, cell_ROI):
	# does a bit of postprocessing: adjust windowing for every channel and sets the foci as active ROI selection
	Stack.setRoi(cell_ROI);
	IJ.run(Stack, "Enhance Contrast...", "saturated=0.3 normalize process_all");
	IJ.run(Stack, "8-bit", "");		

def Cell_Area(IMPStack, filename1):
	# Get Cell Area
	# Pick DAPI channel with black background (pixel value = 0)
	if filename1.endswith('.czi'):
		IMPStack.setSlice(1);
	else:
		IMPStack.setSlice(0);
	
	# Invert to get correct selectio
	cell_mask = IMPStack.crop();
	IJ.run(cell_mask, "Enhance Contrast", "saturated=0.35");
	IJ.run(cell_mask, "8-bit", "slice");
	cell_mask.show();

	# make image binary and create selection
	IJ.setThreshold(1, 255);
	IJ.run(cell_mask, "Convert to Mask", "background=Dark calculate only black");
	IJ.run(cell_mask, "Create Selection", "");

	# Measure Area and get ROI
	Area = round(cell_mask.getStatistics().area,2);
	ROI = cell_mask.getRoi();
	cell_mask.hide();

	return cell_mask, Area, ROI

##############################################################################################################################################
######################################################## MAIN FUNCTION #######################################################################
##############################################################################################################################################

def main():
	#Instructions for users
	gd = GenericDialog("Welcome to our Foci script!")
	gd.addMessage("This script evaluates foci from all the images (.czi or .zvi) in a directory."+ 
	"\nThe first step is to estimate the noise tolerance - this value affects the detections of near by foci and the threshold - this will be the value under which foci will not be taken into account." +
	"\nYou can select either choosing your own threshold or it can be estimated from a selected image or from each image in a folder." +
	"\nIf you choose to estimate the threshold from an image, you will be able to select an area on the image from which the threshold will be estimated."+
	"\nJust follow the instructions given on each dialog and you will manage perfectly")
	gd.addNumericField("Insert the noise tolerance", 100, 1)
	gd.addChoice("Choose an option among the list to get information about the threshold", ["I want to insert the threshold here manually to be used for all the figures. Fill it below.", "I want to pick the threshold for EACH INDIVIDUAL figure using a ROI selection.", "I want to pick the threshold for ONE figure using a ROI selection and then use it for the other figures."], "I want to pick the threshold for ONE figure using a ROI selection and then use it for the other figures.") 
	gd.addNumericField("Manual Threshold", 500, 1)
	gd.addMessage("If you want to analyse the same cells that you have already analysed once and you have the ROI manager with the positions of these cells saved, then check the next field.")
	gd.addCheckbox("I already have the ROI manager saved and I want to analyse these cells again.", False)
	gd.showDialog()

	if gd.wasCanceled(): 
		gd = GenericDialog("Error") 
		gd.addMessage("You have to set a noise tolerance and threshold to run the script. Run the script again and set these parameters.")
		gd.hideCancelButton()		
		gd.showDialog()
		return 0   
	elif gd.wasOKed():
		noise_tolerance = gd.getNextNumber()
		global noise_tolerance
		selection = gd.getNextChoiceIndex()
		global selection
		manualTres = gd.getNextNumber()
		RoiManHave = gd.getNextBoolean()
		global RoiManHave
		
	#Pick the image to analyze
	filename1 = IJ.getFilePath('Pick an image (.czi or .zvi) from the directory you want to analyze');
	
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
		
	AnalysisDir = os.path.dirname(filename1); # the directory with images acquired with the same parameters	

	#Creating directory to collect all the CSV files
	resultpath2 = os.path.join(os.path.dirname(filename1), "Results")
	resultpath2 = os.path.normpath(resultpath2)	

	if not os.path.isdir(resultpath2):
		try:
			os.mkdir(resultpath2)
		except:
			gd = GenericDialog("Error") 
			gd.addMessage("Result directory could not be created. Check filename for special characters.")		
			gd.hideCancelButton()
			gd.showDialog()
			return 0	

	Imagespath = os.path.join(os.path.dirname(filename1), "Analysed_cells")
	Imagespath = os.path.normpath(Imagespath)	

	if not os.path.isdir(Imagespath):
		try:
			os.mkdir(Imagespath)
		except:
			gd = GenericDialog("Error") 
			gd.addMessage("Result directory could not be created. Check filename for special characters.")		
			gd.hideCancelButton()
			gd.showDialog()
			return 0	
	
	# Threshold establishment
	if selection == 0:
		threshold = manualTres
		global threshold	
	elif selection == 2:
		# Open the image
		options = ImporterOptions();
		options.setId(filename1);
		options.setSplitChannels(True);
			
		imps = BF.openImagePlus(options);
		if filename1.endswith('.czi'):	
			imps[1].setTitle("DAPI");
			imps[0].setTitle("GFP");
		else: 
			imps[0].setTitle("DAPI");
			imps[1].setTitle("GFP");
	
		# Establish the threshold 
		if filename1.endswith('.czi'):	
			threshold = ThresholdEst(imps[0]);
		else: 
			threshold = ThresholdEst(imps[1]);
		global threshold #make it global parameter
			

	#Analysis of the images	
	if filename1.endswith('.czi'):
		workDir = glob.glob(os.path.join(AnalysisDir,'*.czi'))
	else:
		workDir = glob.glob(os.path.join(AnalysisDir,'*.zvi'))
		
	for filename in workDir:
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
			imps[0].setTitle("GFP")
		else: 
			imps[0].setTitle("DAPI")
			imps[1].setTitle("GFP")

		if selection == 1:
			if filename1.endswith('.czi'):	
				threshold = ThresholdEst(imps[0]);
			else: 
				threshold = ThresholdEst(imps[1]);
			global threshold #make it global parameter	
	
		resultpath = os.path.join(Imagespath, os.path.basename(filename)[:-4] + "_analyzed_cells_MaxFind_noiseTol_" + str(noise_tolerance)+ '_Threshold_' + str(threshold)) #for images
		resultpath = os.path.normpath(resultpath)	
				
		if not os.path.isdir(resultpath):
			try:
				os.mkdir(resultpath)
			except:
				gd = GenericDialog("Error") 
				gd.addMessage("Result directory could not be created. Check filename for special characters.")		
				gd.hideCancelButton()
				gd.showDialog()
				return 0
		
		if RoiManHave == True:
			resultpathROI = os.path.join(os.path.dirname(filename), os.path.basename(filename)[:-4] +"_ROI.zip"); #ROI saved in the same directory as images
			rm = opensavedROIman(resultpathROI)
		elif RoiManHave == False:
			resultpathROI = os.path.join(os.path.dirname(filename), os.path.basename(filename)[:-4] +"_ROI"); #ROI saved in the same directory as images				
		# Process DAPI channel - make the segmentation and select the cell to be analysed
			if filename1.endswith('.czi'):	
				Cell_Map, Cell_ROIs = Cell_Segmentation(imps[1])
			else:
				Cell_Map, Cell_ROIs = Cell_Segmentation(imps[0])
		
			if filename1.endswith('.czi'):	
				Cell_Map, rm = Pick_Cells(Cell_Map, imps[1], rm, resultpathROI)
			else:
				Cell_Map, rm = Pick_Cells(Cell_Map, imps[0], rm, resultpathROI)

		# open a csv file to write all the results.
		with open(resultpath2 + '\\' + os.path.basename(filename)[:-4]+'_results_MaxFind_noiseTol_' + str(noise_tolerance) + '_UserThreshold_' + str(threshold) + '.csv', 'wb') as csvfile:
			spamwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL);
			spamwriter.writerow(['Cell No', 'Cell Area', '% Area above threshold', 'Min intensity', 'Max intensity', 'Mean intensity', 'N_Foci', 'Foci/Area']);
	
			Data = [];
			N = rm.getCount();
	
			ROIs = rm.getRoisAsArray();
			rm.reset();
			
			# browse all selected cells
			for ROI in ROIs:
	
				rm.addRoi(ROI);
				# copy cropped cell to separate image (all channels)
				Stack = ImageStack(int(ROI.getFloatWidth()),
								   int(ROI.getFloatHeight()));
				
				# browse all channels
				for imp in imps:
					# duplicate image, clear area around ROI and add the remains to the Stack
					# Why? If cells are close to each other, a rectangular image celection around one cell 
					# might include areas of other cells
					imp_buffer = imp.duplicate();
					imp_buffer.setRoi(ROI);
					IJ.setAutoThreshold(imp_buffer, "Default dark reset")
					IJ.run(imp_buffer, "Clear Outside", ""); 
					img = imp_buffer.duplicate();
					imp_buffer.close();
	
					# add channels to Stack
					ip = img.getProcessor();
					Stack.addSlice(img.getTitle(), ip);
					
				# Turn Stack into ImagePlus
				Stack = ImagePlus(ROI.getName(), Stack);
				IJ.run(Stack, "Set Scale...", "distance=6.2016 known=1 pixel=1 unit=micron global");
	
				# Get Cell Area and ROI of Cell
				Stack.show();
				Cell_mask, Area_Cell, Cell_ROI = Cell_Area(Stack, filename1);
				
				# Get Maxima/foci
				Points, findmaxinvert, meanV, minV, maxV, PecentAreaThres = findFoci(Stack, Cell_ROI, threshold, filename1);

				# Foci/Area
				if Points == "OE":
					FA = "OE"
				else:	
					FA = int(Points) / float(Area_Cell)
	
				#Postprocess
				Postprocess(Stack, Cell_ROI);
	
				#Add another slice to the stack with the points that MaximumFinder took into account			
				Stack.getStack().addSlice('Found Maxima', findmaxinvert.getProcessor());			
											
				# Save image of cell
				Stack.show();
				IJ.saveAs("Tiff", resultpath + "\\" + Stack.getTitle());
							
				# If parameter verification is disabled: Close cell image.
				if verify_params:
					Stack.setSlice(5);
					IJ.run(Stack, "Options...", "black");
				else:
					Stack.close();
				# Write to results file
				spamwriter.writerow([Stack.getTitle()[:-4], Area_Cell, PecentAreaThres, minV, maxV, meanV, Points, FA]);
			rm.runCommand("reset")
	gd = GenericDialog("Progress") 
	gd.addMessage("All the images from the selected folder have been evaluated. Good job!")
	gd.showDialog()

main();