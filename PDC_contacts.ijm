//Macro used to measure proximity and/or contact of PDC cells vis-à-vis infected cells
//Needs the Stardist plugin installed : https://github.com/stardist/stardist-imagej
//@Jacques Brocard for Garima Joshi & Marlène Dreux (CIRI), 2024

//STARTS with an open file of 5 channels already analyzed via one of the following macros:
//"PDC_PROX_PLC3" or "PDC_PROX_HepG2"
//CTV = pDC ; BF ; AF488 = ORF2 ; CMPTx Red = All HepG2 cells ; AF647 = LFAI

//---INITIALIZATION
roiManager("Reset");
t=getTitle();
t=substring(t,0,lengthOf(t)-4);
dir=getDirectory("image");

//---MAKE MASK OF INFECTED CELLS and record their centroids
run("Duplicate...", "title=temp duplicate channels=1");
getPixelSize(unit, pixelWidth, pixelHeight);
run("Properties...", "pixel_width=1 pixel_height=1 voxel_depth=1 unit=pixel");
run("Select All");
setForegroundColor(0,0,0);
run("Fill", "slice");

open(dir+t+"_infected.zip");
setForegroundColor(255,255,255);
roiManager("Fill");
nb_infected=roiManager("Count");
run("Set Measurements...", "centroid redirect=None decimal=3");
roiManager("Measure");
xcentro_infected=Table.getColumn("X");
ycentro_infected=Table.getColumn("Y");
selectWindow("Results");
run("Close");

//---OPEN ROIs of pdc and record (i) their interface with infected cells (=contacts)
//or (ii) their proximity to the closest infected cell if smaller than 1 micron
roiManager("Reset");
open(dir+t+"_pdc.zip");
nb_pdc=roiManager("Count");
nb_pdc_inter=0
interface=newArray();
nb_pdc_prox=0;
dist=newArray();
for (r=0;r<nb_pdc;r++){
	roiManager("Select",r);
	getStatistics(area, mean);
	if (mean>0){
		nb_pdc_inter++;
		interface[nb_pdc_inter-1]=mean/255;
	}else{
		setForegroundColor(255,255,255);
		roiManager("Fill");
		roiManager("Measure");
		xcentro=getResult("X",r);
		ycentro=getResult("Y",r);
		selectWindow("Results");
		run("Close");
		inter_min=1000;
		for (i=0;i<nb_infected;i++){
			makeLine(xcentro, ycentro, xcentro_infected[i], ycentro_infected[i]);
			prof=getProfile();
			inter=0;
			for (p=0;p<prof.length;p++){
				if (prof[p]==0) inter++;
			}
			if (inter<inter_min) inter_min=inter;
		}
		//print(r,inter_min);
		roiManager("Select",r);
		setForegroundColor(0,0,0);
		roiManager("Fill");	
		
		if (inter_min<1/pixelWidth){
			nb_pdc_prox++;
			dist[nb_pdc_prox-1]=inter_min*pixelWidth;
		}
	}
}
close("temp");

//---PRINT DETAILED RESULTS
print(t);
print(nb_pdc_prox, " pdc prox w/o contact");
Array.getStatistics(dist, min, max, mean);
if (nb_pdc_prox>0) print("dist moy (µm) =",mean);
print("");
print(nb_pdc_inter, "contacting pdcs (µm²):");
if (nb_pdc_inter>0) for (d=0;d<nb_pdc_inter;d++) print(interface[d]*pixelWidth*pixelWidth);
selectWindow("Log");
saveAs("Txt", dir+t+"_contacts.txt");
