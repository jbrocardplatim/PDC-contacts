//Macro used to measure distribution of infected cells in proximity of PDC cells
//Needs the Stardist plugin installed : https://github.com/stardist/stardist-imagej
//@Jacques Brocard for Garima Joshi & Marlène Dreux (CIRI), 2023

//STARTS with an open .czi file of 5 channels : 
//CTV = pDC ; BF ; AF488 = ORF2 ; CMPTx Red = All HepG2 cells ; AF647 = LFAI

//Thresholds used to detect infected cells; sizes in pixels
infected_min_size=250;
infected_max_size=2500;
infected_min_intensity=250;
pdc_min_intensity=800;
pdc_min_size=25;
prox=3; //nb microns considered as proximity between pdc and infected cells

//--- INITIALIZATION
list = getList("image.titles");
if (list.length==0){
	waitForUser("NO IMAGES OPEN!");
	exit;
}
t=getTitle();
if (substring(t,lengthOf(t)-4)!=".czi"){
	i=0;
	while (substring(t,i,i+4)!=".czi") i++;
	t=substring(t,0,i);
}else{t=substring(t,0,lengthOf(t)-4);}
dir=getDirectory("image");
//Saves a copy of MIP channels
run("Z Project...", "projection=[Max Intensity]");
saveAs("Tiff", dir+t+".tif");
rename(t);
close(t+".czi");
getPixelSize(unit, pixelWidth, pixelHeight);
cyc=floor(prox/pixelWidth);
run("Properties...", "pixel_width=1 pixel_height=1 voxel_depth=1 unit=pixel");
run("Set Measurements...", "area mean centroid redirect=None decimal=3");
print(t);

//--- SEGMENTATION of all cells from channel 4 = red
if (!(isOpen("ROI Manager"))){
	for (s=3;s<=4;s++){
		selectWindow(t);
		setSlice(s);
		run("Duplicate...", "title="+s);
		selectWindow(s);
		run("Smooth");
		run("Mean...", "radius=5");
		run("32-bit");
		//run("Log");
		run("Enhance Contrast", "saturated=0.35");
		setOption("ScaleConversions", true);
		run("8-bit");
	}
	imageCalculator("Average create", "3","4");
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'Result of 3', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.8', 'probThresh':'0.5', 'nmsThresh':'0.15', 'outputType':'ROI Manager', 'nTiles':'9', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	close();
	exclude_roi_size(infected_min_size, infected_max_size);
	close("3"); close("4");
}
//Macro stops if no cells have been detected
nROIs=roiManager("Count");
if (nROIs<2){
	waitForUser("No cells detected ; macro will stop");
	exit;
}

//--- SEGMENTATION of infected cells from channel 3 = green
selectWindow(t);
run("Select All");
setSlice(3);
run("Duplicate...", "title=green");
run("Smooth");
nb_infected_cells=0;
nb_cells=roiManager("Count");
for (r=0;r<nb_cells;r++){
	roiManager("Select",r);
	RoiManager.setGroup(0);
	getStatistics(area, mean);
	
	if (mean>infected_min_intensity){
		nb_infected_cells++;
		RoiManager.setGroup(1);
	}
}
roiManager("Deselect");
roiManager("Save",dir+t+"_all_cells.zip");
print("Total : " +nb_cells+" cells, including "+ nb_infected_cells+" infected cells");

//--- DISTANCE MAP of infected cells
selectWindow("green");
setForegroundColor(0,0,0);
run("Select All");
run("Fill");
run("8-bit");
run("Duplicate...", "title=dist");
setBatchMode(true);
setForegroundColor(255,255,255);
selectWindow("dist");
for (r=0;r<nb_cells;r++){
	roiManager("Select",0);
	if (Roi.getGroup()==1) roiManager("Add");
	roiManager("Fill");
	roiManager("Select",0);
	roiManager("Delete");
}
run("Select All");
run("Invert");
run("Distance Map");
if (nb_infected_cells>0) roiManager("Save",dir+t+"_infected.zip");

//--- SEGMENTATION of PDCs from channel 1 = CTV
roiManager("Reset");
selectWindow(t);
setSlice(1);
run("Duplicate...", "title=temp");
setThreshold(pdc_min_intensity,66666);
setOption("BlackBackground", true);
run("Convert to Mask");
run("Options...", "iterations=1 count=1 black pad do=Erode");
run("Options...", "iterations=1 count=1 black pad do=Open");
run("Fill Holes");
run("Watershed");
run("Analyze Particles...", "size="+pdc_min_size+"-Infinity show=Masks exclude clear add");
run("Grays");
rename("pdc");
close("temp");
roiManager("Deselect");
nPDC=roiManager("Count");

//--- COUNT PDCs in proximity of any cell
if (nPDC>0){
	nPDC_prox=0;
	selectWindow("dist");
	for (p=0;p<nPDC;p++){
		roiManager("Select",p);
		getStatistics(area, mean, min);
		if (min<=cyc) {
			nPDC_prox++;
			RoiManager.setGroup(1);
		}
	}
	roiManager("Save",dir+t+"_pdc.zip");
	close("dist");
	roiManager("Reset");
}
print("Total : " +nPDC+" PDCs, including "+ nPDC_prox+" within "+prox+" microns of any cell");

//--- FOR EACH INFECTED CELL, count PDCs in proximity
if (nb_infected_cells>0){
	roiManager("Open",dir+t+"_infected.zip");
	selectWindow("green");
	run("Select All");
	run("Copy");
	roiManager("Measure");
	centroX=Table.getColumn("X");
	centroY=Table.getColumn("Y");
	pdc_prox=newArray(nb_infected_cells);
	close("Results");
	
	setForegroundColor(255,255,255);
	for (r=0;r<nb_infected_cells;r++){
		selectWindow("green");
		run("Select All");
		run("Paste");
		roiManager("Select",r);
		roiManager("Fill");
		makeRectangle(centroX[r]-100, centroY[r]-100, 200, 200);
		run("Duplicate...", "title=temp");
		run("Options...", "iterations="+cyc+" count=1 black pad do=Dilate");
		run("Analyze Particles...", "add");
		close("temp");
		selectWindow("pdc");
		makeRectangle(centroX[r]-100, centroY[r]-100, 200, 200);
		run("Duplicate...", "title="+r);
		roiManager("Select",nb_infected_cells);
		run("Analyze Particles...", "summarize");
		pdc_prox[r]=Table.get("Count",r);
		roiManager("Select",nb_infected_cells);
		roiManager("Delete");
		close(r);
	}
	close("Summary");
}

//--- Color-Code EACH INFECTED CELL, based on #PDCs in proximity
pdc_prox1=0;
pdc_prox2=0;
pdc_prox3=0;
if (nb_infected_cells>0){
	for (n=0;n<nb_infected_cells;n++){
		if (pdc_prox[n]==1) pdc_prox1++;
		if (pdc_prox[n]==2) pdc_prox2++;
		if (pdc_prox[n]>=3){
			pdc_prox3++;
			pdc_prox[n]=3;
		}
		roiManager("Select",n)
		RoiManager.setGroup(pdc_prox[n]);
	}
}
close("pdc");
close("green");
setBatchMode("exit and display");
//print("# infected cells with PDC cells in proximity");
pdc_prox0=nb_infected_cells-pdc_prox1-pdc_prox2-pdc_prox3;
print("# infected cells : \t"+pdc_prox0+" \t"+pdc_prox1+" \t"+pdc_prox2+" \t"+pdc_prox3);
print("with #PDC cells in proximity : \t0 \t1 \t2 \t3+");
if (nb_infected_cells>0) roiManager("Save",dir+t+"_infected.zip");
close("ROI Manager");
selectWindow("Log");
saveAs("Txt", dir+t+".txt");
close();


function exclude_roi_size(min, max){
	nROIs=roiManager("Count");
	for (r=0;r<nROIs;r++){
		roiManager("Select",0);
		getStatistics(area);
		if (area>=min && area<=max) roiManager("Add");
		roiManager("Delete");
	}
	roiManager("Deselect");
	roiManager("Sort");
}
