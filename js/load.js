/*==========================================================
% MoDMaps3D - THREE DIMENSIONAL MOLECULAR DISTANCE MAPS 
% Version 3
% Coded by Rallis Karamichalis, 2016.
% ----------------------------------------------------------
% This is a vizualization tool for creating 
% 3D Molecular Distance Maps (MoDMaps3D). 
% ==========================================================*/


var MoDMaps3D = { 'version':'3.0', 
		'dbg':false, 
		'gglSearchEnabled':true, 
		'realtimeHighlight':true, 
		'highlightColor':parseInt("0x28EE2F") };
console.log(MoDMaps3D);

var dbg = MoDMaps3D['dbg'];                   
var gglSearchEnabled = MoDMaps3D['gglSearchEnabled']; 
var realtimeHighlight = MoDMaps3D['realtimeHighlight'];
var highlightColor = MoDMaps3D['highlightColor'];   
var mapsWithReadyFCGRs = {'mtDNA_Amphibians.txt':'mtDNA', 'mtDNA_Hominidae_Hylobatidae.txt':'mtDNA', 'mtDNA_InsecMamAmph.txt':'mtDNA', 'mtDNA_Insects.txt':'mtDNA', 'mtDNA_Mammals.txt':'mtDNA', 'mtDNA_Primates.txt':'mtDNA', 'mtDNA_Vertebrata.txt':'mtDNA'};
var mapsWithReadyDistMatrix = {'mapfilename':'distMatrixFilename'};
var mapid, dim1, dim2, dim3, radius, alldata, distMatrix = [];
var setOfPoints, colors, numberOfLabels, namesOfLabels, legendColors, legendLabels; 
var globalScaledPointsCoord, globalPointsLabels;      
var intersectedPoint, geometries, materials, meshes, offsets, scene;         
var renderer, camera, mouse, raycaster, meshesRotX, meshesRotY, cameraX, cameraY, cameraZ; 
var selectedIndex, minDistance, allHits;   // FOR RAYCASTING aka 'selection' BY USER 
var selectedGeometry, selectedMaterial, selectedMesh, selectedPoint, protoSelectedSphereGeometry, protoVertexCount, searchOffsets, axisUnit;                   
var linkNCBI = 'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=nuccore&id=';
var linkAssemblyNCBI = 'http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=assembly&id=';
var allAccessionNums, allNames, allBioInfo, allSearchHits;    
var fromIndex="", toIndex="", dists, howmany;
var dxTimer, dyTimer, dzTimer;   


// CREATE VARIOUS INFO DIVS AND SET THEIR IDs
var pointInfoDiv, fcgrDiv, cgrInfoDiv, searchDiv, distPointsDiv, staticInfoDiv, fcgrInfo;
var separator = document.createElement("hr");  
separator.style.color = "white";
separator.style.width = "60%";
separator.appendChild(document.createTextNode(""));

// ON WINDOW RESIZE, RESIZE MAP
window.addEventListener('resize', onWindowResize, false);
function onWindowResize() {
	camera.aspect = window.innerWidth / window.innerHeight;
	camera.updateProjectionMatrix();
	renderer.setSize(window.innerWidth, window.innerHeight);
}

// GET URL ARGUMENT VALUE
function geturlparamvalue(name) {
	name = name.replace(/[\[]/, "\\\[").replace(/[\]]/, "\\\]");
	var regexS = "[\\?&]" + name + "=([^&#]*)";
	var regex = new RegExp(regexS);
	var results = regex.exec(window.location.href);
	if(results!=null){
		return results[1];
	}else{
		return -1;
	}
}

// SEARCH FUNCTION 
function startSearch(){
	allSearchHits = [];
	var query = $("#tosearch").val();
	var queryLen = query.length;
	if(4-queryLen > 0){
		// IF LESS THAN 4 CHARACTERS HAVE BEEN ENTERED
		scene.remove(selectedMesh);
		var missingLen=4-queryLen;
		if(missingLen==1){ 
			$("#searchstatus").html('Enter 1 more character..'); 
		}else{
			$("#searchstatus").html("Enter "+missingLen+" more characters.."); 
		}
	}else{
		// IF >= 4 CHARACTERS HAVE BEEN ENTERED
		$("#searchstatus").html('');
		var searchOutput="";
		var patt,res1,res2,res3,ind1,ind2,ind3,redcol1,redcol2,redcol3,bfr1,bfr2,bfr3,after1,after2,after3;
		var howManySearchResults=0;
		for(var i=0;i<allAccessionNums.length;i++){
			//console.log("q=", query);
			patt=new RegExp(query,'gi');
			// console.log("patt=", patt);
			res1=allAccessionNums[i].match(patt);
			res2=allNames[i].match(patt);
			if(allBioInfo.length>0){
				res3=allBioInfo[i].match(patt);
			}else{
				res3=null;
			}
			
			if((res1!=null)||(res2!=null)||(res3!=null)){
				allSearchHits.push(i);
				howManySearchResults++;
				if(res1!=null){
					ind1=allAccessionNums[i].toLowerCase().indexOf(query.toLowerCase());
					bfr1=allAccessionNums[i].substr(0,ind1);
					redcol1=allAccessionNums[i].substr(ind1,queryLen);
					after1=allAccessionNums[i].substr(ind1+queryLen);
					searchOutput=searchOutput+'<input type="radio" name="pointFoundFromSearch" value="'+i+'" onclick="selectAndFill('+i+');">['+bfr1+'<strong style="color:red">'+redcol1+'</strong>'+after1+'] <br/> ['+allNames[i]+']<br/>';
				}else if(res2!=null){
					ind2=allNames[i].toLowerCase().indexOf(query.toLowerCase());
					bfr2=allNames[i].substr(0,ind2);
					redcol2=allNames[i].substr(ind2,queryLen);
					after2=allNames[i].substr(ind2+queryLen);
					searchOutput=searchOutput+'<input type="radio" name="pointFoundFromSearch" value="'+i+'" onclick="selectAndFill('+i+');">['+allAccessionNums[i]+'] <br/> ['+bfr2+'<strong style="color:red">'+redcol2+'</strong>'+after2+']<br/>';
				}else if(res3!=null){ 
					ind3=allBioInfo[i].toLowerCase().indexOf(query.toLowerCase());
					bfr3=allBioInfo[i].substr(0,ind3);
					redcol3=allBioInfo[i].substr(ind3,queryLen);
					after3=allBioInfo[i].substr(ind3+queryLen);
					searchOutput=searchOutput+'<input type="radio" name="pointFoundFromSearch" value="'+i+'" onclick="selectAndFill('+i+');">['+allAccessionNums[i]+'] <br/> ['+bfr3+'<strong style="color:red">'+redcol3+'</strong>'+after3+']<br/>';
				} 
			}
		}
		
		if(howManySearchResults==0){
			$("#searchstatus").html('Nothing matches your search.');
		}else if(howManySearchResults<=10){
			$("#searchstatus").html(searchOutput);
		}else{
			$("#searchstatus").html('Found '+howManySearchResults+' matches. To select a specific point, please narrow it down to less than 10 results. <input type="button" value="Plot this!" onclick="saveAndPlot()">');
		}            

		if(realtimeHighlight){
			updateRealTimeSelectedMesh();
		}     
	}
}

// SAVE LOCALLY SELECTED DATASET AND PLOT ITFROM SCRATCH
function saveAndPlot(){
	var out='';
	for(var i=0; i<allSearchHits.length; i++){
		out += allAccessionNums[allSearchHits[i]]+',';
	}
	out = out.substring(0,out.length-1);
	//console.log(out);
	
	var dataTimeStamp = new Date().getTime();
	if(typeof(Storage) !== "undefined") {
		localStorage.setItem("dataset"+dataTimeStamp, out);
	} else {
		console.log("ERROR! NO LOCAL STORAGE AVAILABLE!?");
		alert("ERROR! NO LOCAL STORAGE AVAILABLE!?");
	}
	window.location = "index.html?dataset=local"+dataTimeStamp;
}

// SELECTION OF SEARCH RESULTS
function selectAndFill(id){
	$("#searchstatus").html("");
	$("#tosearch").val(allAccessionNums[id]);
	selectedIndex = id;
	intersectedPoint = {x:offsets[selectedIndex].x, y:offsets[selectedIndex].y, z:offsets[selectedIndex].z};
	console.log("SELECTED point is= "+selectedIndex);
	selectedPoint = intersectedPoint;
	$("#searchstatus").html("Take a closer look, point is selected!!");
	updateSelectedMesh();
	updateInfoDiv();
}

// COMPUTE DIST 
function computeDist(preload){
	if(dbg){console.log("preload= "+preload);}

	// Load distance matrix if it's first time
	if(preload==true){ 
		if(distMatrix.length>0){
			console.log('distMatrix.length > 0');
			howmany=Math.sqrt(distMatrix.length);
			console.log("howmany= "+howmany+" x "+howmany);
			dists = [];
			for(var i=0; i<distMatrix.length; i++){
				if(i%(howmany)==0){ dists.push([]); }
				dists[dists.length-1].push(distMatrix[i]);
			}
			// temp dbg
			// var toprint="[";
			// for(var i=0; i<dists.length; i++){
			// 	toprint+='';
			// 	for (var j=0; j<dists.length; j++){
			// 		toprint+=dists[i][j]+',';
			// 	}
			// 	toprint = toprint.substring(0, toprint.length - 1);
			// 	toprint+="\n"
			// }
			// toprint = toprint.substring(0, toprint.length - 1);
			// toprint+=']';
			// console.log("matrix= ",toprint);
			$("#computeDist").html('');
			distPointsDiv.innerHTML = '<em><strong><font color="yellow" size="4">Distance between Points</font></strong></em><br>\
					<table border="0">\
					<tr><td>From:</td><td><input id="fromHere" value="'+fromIndex+'" maxlength="7" size="7" disabled /></td>\
					<td><input type="button" value="Add" onclick="add(\'from\');"></td></tr>\
					<tr><td>To:</td><td><input id="toHere" value="'+toIndex+'" maxlength="7" size="7" disabled /></td>\
					<td><input type="button" value="Add" onclick="add(\'to\');"></td></tr>\
					<tr><td>Distance:</td><td><input id="outputDist" type="text" value="" maxlength="10" size="7" disabled /></td>\
					<td><div id="computeDist"></div></td></tr></table>\
					<hr color="white" width="60%">';
		}else{
			var xmlhttp;
			xmlhttp=new XMLHttpRequest();
			$("#computeDist").html('Loading, please wait.. <img src="img/loading.gif" height="30">');
			xmlhttp.onreadystatechange=function(){
				if (xmlhttp.readyState==4 && xmlhttp.status==200){
					dists=xmlhttp.responseText.split("\n");
					howmany=dists.length;
					for(var i=0; i<dists.length; i++){
						dists[i] = dists[i].split(',');
					}
					$("#computeDist").html('');
					distPointsDiv.innerHTML = '<em><strong><font color="yellow" size="4">Distance between Points</font></strong></em><br>\
							<table border="0">\
							<tr><td>From:</td><td><input id="fromHere" value="'+fromIndex+'" maxlength="7" size="7" disabled /></td>\
							<td><input type="button" value="Add" onclick="add(\'from\');"></td></tr>\
							<tr><td>To:</td><td><input id="toHere" value="'+toIndex+'" maxlength="7" size="7" disabled /></td>\
							<td><input type="button" value="Add" onclick="add(\'to\');"></td></tr>\
							<tr><td>Distance:</td><td><input id="outputDist" type="text" value="" maxlength="10" size="7" disabled /></td>\
							<td><div id="computeDist"></div></td></tr></table>\
							<hr color="white" width="60%">';
				}else if(xmlhttp.readyState==4 && xmlhttp.status==404){
					console.log("Distance file is not available..");
					distPointsDiv.innerHTML='<em><strong><font color="yellow" size="4">Distance between Points</font></strong></em><br>Distance Matrix has not been made available for this map.<hr color="white" width="60%">';
				}       
			}
			xmlhttp.open("GET","dists/"+mapsWithReadyDistMatrix[mapid],true);
			xmlhttp.send();	
		}
	}else if(preload==false){
		var idFrom=$("#fromHere").val();
		var idTo=$("#toHere").val();
		if((idFrom!="")&&(idTo!="")){
			$("#outputDist").val('');
			$("#computeDist").html('<img src="img/loading.gif" height="30">');
			idFrom=parseFloat(globalPointsLabels[idFrom][0]);
			idTo=parseFloat(globalPointsLabels[idTo][0]);
			// old method, only upper triangular
			// var res=dists[(idFrom-1)*howmany+idTo-1];
			var res=dists[idFrom][idTo];
			if(preload!=true){$("#outputDist").val(res);}
			$("#computeDist").html('');
		}
	}else{
		console.log("Weird call. Please debug it...");
	}
}

// ADD FROM/TO INDICES
function add(place){
	if(selectedIndex!=undefined){
		if(place=='from'){
			$('#fromHere').val(selectedIndex);
			fromIndex=selectedIndex;
		}
		if(place=='to'){
			$('#toHere').val(selectedIndex);
			toIndex=selectedIndex;
		}
	}
	computeDist(false);
}

// UPDATES ALL DIVS THAT REQUIRE DYNAMIC DATA
function updateInfoDiv() {
	
	// POINTINFO DIV
	var infoTable='';
	for(var indLabel=0; indLabel<namesOfLabels.length; indLabel++){
		if(selectedIndex!=undefined){
			if(namesOfLabels[indLabel]=="Acc"){
				var prefix=globalPointsLabels[selectedIndex][indLabel].substring(0,3);
				var link;
				if((prefix=="GCA")||(prefix=="GCF")){
					link=linkAssemblyNCBI+globalPointsLabels[selectedIndex][indLabel];
				}else{
					link=linkNCBI+globalPointsLabels[selectedIndex][indLabel];  
				}
				infoTable=infoTable+'<tr><td>'+namesOfLabels[indLabel]+'</td><td>'+globalPointsLabels[selectedIndex][indLabel]+' (<a href="'+link+'" target="_blank"><font color="yellow">NCBI</font></a>)</td></tr>';
				
				if(mapid in mapsWithReadyFCGRs){
					console.log("readyFCGRS!!");
					fcgrInfo='<a href="fcgrs/'+mapsWithReadyFCGRs[mapid]+'/'+globalPointsLabels[selectedIndex][indLabel]+'.png" target="_blank"><img src="fcgrs/'+mapsWithReadyFCGRs[mapid]+'/'+globalPointsLabels[selectedIndex][indLabel]+'.png" height="200px" width="200px" alt="FCGR IMAGE NOT AVAILABLE" title="Click here to Zoom In"></a>';
				}else{
					fcgrInfo="Not Available.";			 
				}
			}else if(namesOfLabels[indLabel]=="AsmAcc"){
				var link=linkNCBI+globalPointsLabels[selectedIndex][indLabel];
				infoTable=infoTable+'<tr><td>'+namesOfLabels[indLabel]+'</td><td>'+globalPointsLabels[selectedIndex][indLabel]+' (<a href="'+link+'" target="_blank"><font color="yellow">NCBI</font></a>)</td></tr>';
				if(mapid in mapsWithReadyFCGRs){
					console.log("readyFCGRS!!");
					fcgrInfo='<a href="fcgrs/'+mapsWithReadyFCGRs[mapid]+'/'+globalPointsLabels[selectedIndex][indLabel]+'.png" target="_blank"><img src="fcgrs/'+mapsWithReadyFCGRs[mapid]+'/'+globalPointsLabels[selectedIndex][indLabel]+'.png" height="200px" width="200px" alt="FCGR IMAGE NOT AVAILABLE" title="Click here to Zoom In"></a>';
				}else{
					fcgrInfo="Not Available.";				 
				}
			}else if(namesOfLabels[indLabel]=="Name"){
				if(globalPointsLabels[selectedIndex][indLabel].substring(0,13)=='mitochondrion'){
					var link='https://www.google.com/search?hl=en&site=imghp&tbm=isch&source=hp&q='+globalPointsLabels[selectedIndex][indLabel].substring(14,globalPointsLabels[selectedIndex][indLabel].length);
				}else{
					var link='https://www.google.com/search?hl=en&site=imghp&tbm=isch&source=hp&q='+globalPointsLabels[selectedIndex][indLabel];  
				}
				if(gglSearchEnabled){
					infoTable=infoTable+'<tr><td>'+namesOfLabels[indLabel]+'</td><td>'+globalPointsLabels[selectedIndex][indLabel]+' (<a href="'+link+'" target="_blank"><font color="yellow">Google it!</font></a>)</td></tr>';   
				}else{
					infoTable=infoTable+'<tr><td>'+namesOfLabels[indLabel]+'</td><td>'+globalPointsLabels[selectedIndex][indLabel]+'</td></tr>';
				}	
			}else if(namesOfLabels[indLabel]=="Taxa"){
				infoTable=infoTable+'<tr><td>'+namesOfLabels[indLabel]+'</td><td><div id="taxa">'+globalPointsLabels[selectedIndex][indLabel]+'</div></td></tr>'; 
			}else if(namesOfLabels[indLabel]=="Index"){
				infoTable=infoTable+'<tr><td>'+namesOfLabels[indLabel]+'</td><td>'+globalPointsLabels[selectedIndex][indLabel]+' ('+selectedIndex+')</td></tr>';
			}else{
				infoTable=infoTable+'<tr><td>'+namesOfLabels[indLabel]+'</td><td>'+globalPointsLabels[selectedIndex][indLabel]+'</td></tr>';        
			}
		}else{
			infoTable=infoTable+'<tr><td>'+namesOfLabels[indLabel]+'</td><td>Undefined</td></tr>';
			fcgrInfo='Unavailable for this map!';
		}
	}
	pointInfoDiv.innerHTML= '<em><strong><font color="yellow" size="4">Point Info</font></strong></em><table border="3">'+infoTable+'</table>';
	if(document.getElementById('taxa')!=undefined){
		document.getElementById('taxa').style.display='block';    
	}
	
	// FCGRINFO DIV
	fcgrInfo+='<hr color="white" width="60%">';
	fcgrDiv.innerHTML=fcgrInfo;

	// DISTPOINTS DIV
	if(document.getElementById("fromHere")!=undefined){$("#fromHere").val(fromIndex);}
	if(document.getElementById("toHere")!=undefined){$("#toHere").val(toIndex);}
}

// CLONE GEOMETRY
function cloneGeometry(from, to, vertexOffsets, vertexCountOffset) {
	for (var i in from.vertices) {
		var v = from.vertices[i];
		to.vertices.push(new THREE.Vector3(v.x + vertexOffsets.x, v.y + vertexOffsets.y, v.z + vertexOffsets.z));  
	}
	for (var i in from.faces) {
		var f = from.faces[i];
		to.faces.push(new THREE.Face3(f.a + vertexCountOffset, f.b + vertexCountOffset, f.c + vertexCountOffset));
	}
}

// UPDATE SELECTED MESH
function updateSelectedMesh() {
	scene.remove(selectedMesh);
	selectedGeometry = new THREE.Geometry();
	var selectedCount = 0;
	if (selectedPoint) {
		cloneGeometry(protoSelectedSphereGeometry, selectedGeometry, selectedPoint, protoSelectedSphereGeometry.vertices.length * selectedCount);
		selectedCount++;
	}
	selectedGeometry.mergeVertices();
	selectedGeometry.computeFaceNormals();
	selectedGeometry.computeVertexNormals();

	//selectedGeometry.computeFaceNormals();
	selectedMaterial = new THREE.MeshLambertMaterial({color: highlightColor, wireframe:false, shading: THREE.SmoothShading});
	selectedMesh = new THREE.Mesh(selectedGeometry, selectedMaterial);
	selectedMesh.rotation.x = meshes[0].rotation.x;
	selectedMesh.rotation.y = meshes[0].rotation.y;
	scene.add(selectedMesh);
}

// HIGHLIGHT SEARCH RESULTS REAL TIME
function updateRealTimeSelectedMesh(){
	scene.remove(selectedMesh);
	selectedGeometry = new THREE.Geometry();
	searchOffsets=[];
	allSearchHits
	for(var i = 0; i<allSearchHits.length; i++){
		searchOffsets.push({
			x:20*globalScaledPointsCoord[allSearchHits[i]][0]*axisUnit, 
			y:20*globalScaledPointsCoord[allSearchHits[i]][1]*axisUnit, 
			z:20*globalScaledPointsCoord[allSearchHits[i]][2]*axisUnit
		});
	}
	var tmpInd=0;
	for (var r = 0; r<allSearchHits.length; r++) {
		cloneGeometry(protoSphereGeometry, selectedGeometry, searchOffsets[r], protoSelectedSphereGeometry.vertices.length * tmpInd);
		tmpInd++;
	}
	selectedGeometry.mergeVertices();
	selectedGeometry.computeFaceNormals();
	selectedGeometry.computeVertexNormals();

	selectedMaterial = new THREE.MeshLambertMaterial({color: highlightColor, wireframe:false, shading: THREE.SmoothShading});
	selectedMesh = new THREE.Mesh(selectedGeometry, selectedMaterial);
	selectedMesh.rotation.x = meshes[0].rotation.x;
	selectedMesh.rotation.y = meshes[0].rotation.y;
	scene.add(selectedMesh);

	selectedIndex=undefined;
	updateInfoDiv();
}

// TOGGLE VISIBILITY OF AN ELEMENT
function toggle(id){
	$("#"+id).slideToggle();
}

// TOGGLE LEFT PANEL
function toggleleft(){
	$("#leftmenu").slideToggle();
	$("#leftminimized").slideToggle();
}

// TOGGLE RIGHT PANEL
function toggleright(){
	$("#rightmenu").slideToggle();
	$("#rightminimized").slideToggle();
}

// CHANGE HIGHLIGHT COLOR
function changeHighlightColor(){
	highlightColor= parseInt('0x'+ $("#pickbtn").val());
	MoDMaps3D['highlightColor'] = highlightColor;
	alert("Highlight color has changed!\nClick on any point to see it in action!");
}

// REDRAW MAP WITH NEW SETTINGS AS SET BY USER
function redraw(){
	dim1selected=$("#dim1").val();
	dim2selected=$("#dim2").val();
	dim3selected=$("#dim3").val();
	radius=$("#radius").val();
	if(dbg){console.log("You selected: "+dim1selected+"-"+dim2selected+"-"+dim3selected+" and radius="+radius);}
	var camX = camera.position.x;
	var camY = camera.position.y;
	var camZ = camera.position.z;
	var meshesX= meshes[0].rotation.x;
	var meshesY= meshes[0].rotation.y; 

	// AUTONAVIGATE ONLY IF NOT CHANGING THE DIMENSIONS USED
	var redirectTo;
	if((dim1==dim1selected) && (dim2==dim2selected) && (dim3==dim3selected)){
		redirectTo="load.html?mapid="+mapid +"&dim1="+dim1+"&dim2="+dim2+"&dim3="+dim3+"&radius="+radius+"&cameraX="+camX+"&cameraY="+camY+"&cameraZ="+camZ+"&meshesRotX="+meshesX+"&meshesRotY="+meshesY+"&autonavigate=true";
	}else{
		redirectTo="load.html?mapid="+mapid +"&dim1="+dim1selected+"&dim2="+dim2selected+"&dim3="+dim3selected+"&radius="+radius;
	} 
		
	window.location = redirectTo;
}

// AUTONAVIGATION OF CAMERA
function myCameraMove(dim, goal) {
	var d=0.1;
	if(dim=='x'){
		dif = camera.position.x - goal;
		if(dif>=0){camera.position.x -=d;}else{camera.position.x +=d;}
		if(Math.abs(dif)<d){clearInterval(dxTimer);}
	}
	if(dim=='y'){
		dif = camera.position.y - goal;
		if(dif>=0){camera.position.y -=d;}else{camera.position.y +=d;}
		if(Math.abs(dif)<d){clearInterval(dyTimer);}
	}
	if(dim=='z'){
		dif = camera.position.z - goal;
		if(dif>=0){camera.position.z -=d;}else{camera.position.z +=d;}
		if(Math.abs(dif)<d){clearInterval(dzTimer);}
	}
	if(dim=='xRot'){
		dif = meshes[0].rotation.x - goal;
		for(var meshInd=0; meshInd<meshes.length; meshInd++){
			if(dif>=0){meshes[meshInd].rotation.x -=d;}else{meshes[meshInd].rotation.x +=d;}  
		}
		if(dif>=0){selectedMesh.rotation.x -=d;}else{selectedMesh.rotation.x +=d;}
		if(Math.abs(dif)<d){clearInterval(dxRotTimer);}
	}
	if(dim=='yRot'){
		dif = meshes[0].rotation.y - goal;
		for(var meshInd=0; meshInd<meshes.length; meshInd++){
			if(dif>=0){meshes[meshInd].rotation.y -=d;}else{meshes[meshInd].rotation.y +=d;}  
		}
		if(dif>=0){selectedMesh.rotation.y -=d;}else{selectedMesh.rotation.y +=d;}
		if(Math.abs(dif)<d){clearInterval(dyRotTimer);}
	}
}

// DRAW 3D GRAPHICS - STEP 3
function initGraphics(){
	if(dbg){console.log("begin initGraphics()..");}

	// remove loading image 
	$('#showload').remove();
	
	scene = new THREE.Scene();
	camera = new THREE.PerspectiveCamera(45, window.innerWidth/window.innerHeight, 0.1, 1000);
	mouse = new THREE.Vector2();
	raycaster = new THREE.Raycaster();
	renderer = new THREE.WebGLRenderer();
	renderer.setSize(window.innerWidth, window.innerHeight);
	document.body.appendChild(renderer.domElement);
	renderer.setClearColor( 0xC7C7C7, 1); //change color  E9E9E9   A1A1A1  C7C7C7

	selectedPoint = null;
	axisUnit = 0.15; 

	protoSphereGeometry = new THREE.SphereGeometry(radius);
	protoSelectedSphereGeometry = new THREE.SphereGeometry(1.7*radius);
	protoVertexCount = protoSphereGeometry.vertices.length;

	geometries = [];
	materials= [];
	meshes = []; 
	offsets = [];
	var s=0;
	var startInd;
	for(var i=0; i<setOfPoints.length; i++){
		var tmpGeometry = new THREE.Geometry();
		geometries.push(tmpGeometry);
		startInd=s;
		for(var j=0; j<setOfPoints[i]; j++){
			offsets.push({
				x:20*globalScaledPointsCoord[s][0]*axisUnit, 
				y:20*globalScaledPointsCoord[s][1]*axisUnit, 
				z:20*globalScaledPointsCoord[s][2]*axisUnit}); 
			s++;
		}

		var tmpInd=0;
		for (var r = startInd; r<s; r++) {
				cloneGeometry(protoSphereGeometry, geometries[i], offsets[r], protoVertexCount * tmpInd);
				tmpInd++;
		}
		geometries[i].mergeVertices();
		geometries[i].computeFaceNormals();
		geometries[i].computeVertexNormals();

		var tmpMaterial = new THREE.MeshLambertMaterial({
			color: colors[i], 
			wireframe:false, 
			shading: THREE.SmoothShading
		});
		materials.push(tmpMaterial);
		var tmpMesh = new THREE.Mesh(geometries[i], materials[i]);
		meshes.push(tmpMesh);
		scene.add(meshes[i]);
	}

	var facedownGeom = new THREE.Geometry();
	var faceupGeom = new THREE.Geometry();
	var line1Geom = new THREE.Geometry();
	var line2Geom = new THREE.Geometry();
	var line3Geom = new THREE.Geometry();
	var line4Geom = new THREE.Geometry();
	var cubeMaterial = new THREE.LineBasicMaterial({color: 0x000000});
		
	facedownGeom.vertices.push(new THREE.Vector3( -3.5, -3.5, -3.5 ),new THREE.Vector3( 3.5, -3.5, -3.5 ),new THREE.Vector3( 3.5, 3.5, -3.5 ),new THREE.Vector3( -3.5, 3.5, -3.5 ),new THREE.Vector3( -3.5, -3.5, -3.5 ));
	faceupGeom.vertices.push(new THREE.Vector3( -3.5, -3.5, 3.5 ),new THREE.Vector3( 3.5, -3.5, 3.5 ), new THREE.Vector3( 3.5, 3.5, 3.5 ),new THREE.Vector3( -3.5, 3.5, 3.5 ),new THREE.Vector3( -3.5, -3.5, 3.5 ));
	line1Geom.vertices.push(new THREE.Vector3( -3.5, -3.5, -3.5 ),new THREE.Vector3( -3.5, -3.5, 3.5 ));
	line2Geom.vertices.push(new THREE.Vector3( -3.5, 3.5, -3.5 ),new THREE.Vector3( -3.5, 3.5, 3.5 ));
	line3Geom.vertices.push(new THREE.Vector3( 3.5, -3.5, -3.5 ),new THREE.Vector3( 3.5, -3.5, 3.5 ));
	line4Geom.vertices.push(new THREE.Vector3( 3.5, 3.5, -3.5 ),new THREE.Vector3( 3.5, 3.5, 3.5 ));

	var facedown = new THREE.Line( facedownGeom, cubeMaterial );
	var faceup = new THREE.Line( faceupGeom, cubeMaterial );
	var line1 = new THREE.Line( line1Geom, cubeMaterial );
	var line2= new THREE.Line( line2Geom, cubeMaterial );
	var line3 = new THREE.Line( line3Geom, cubeMaterial );
	var line4 = new THREE.Line( line4Geom, cubeMaterial );

	meshes.push(facedown); 
	meshes.push(faceup); 
	meshes.push(line1); 
	meshes.push(line2); 
	meshes.push(line3); 
	meshes.push(line4);
	scene.add(facedown); 
	scene.add(faceup);
	scene.add(line1); 
	scene.add(line2); 
	scene.add(line3); 
	scene.add(line4);

	selectedGeometry = new THREE.Geometry();  
	selectedMaterial = new THREE.MeshLambertMaterial({
		color: highlightColor, 
		wireframe:false, shading: THREE.SmoothShading
	}); 
	selectedMesh = new THREE.Mesh(selectedGeometry, selectedMaterial);
	scene.add(selectedMesh);

	camera.position.x = 0;
	camera.position.y = 0;
	camera.position.z = 15;

	var ambientLight = new THREE.AmbientLight(0x444444);
	scene.add(ambientLight);
	var directionalLight = new THREE.DirectionalLight( 0xffffff, 1 );
	directionalLight.position.set( 1, 1, 1 );
	scene.add( directionalLight );

	var render = function () {
		requestAnimationFrame(render);
		renderer.render(scene, camera);
	};
		
	render();
	pointInfoDiv = document.getElementById('pointInfoDiv');
	fcgrDiv = document.getElementById('fcgrDiv');
	cgrInfoDiv = document.getElementById('cgrInfoDiv');
	searchDiv = document.getElementById('searchDiv');
	distPointsDiv = document.getElementById('distPointsDiv');

	$("#leftmenu").css({
		"display": "block",
		"position": "absolute", 
		"padding": "5px 8px 5px 8px",
		"background-color": "black", 
		"color": "white",
		"top": "0px",
		"bottom": "0px",
		"overflow-y": "auto",
		//"overflow-x": "auto",
		"max-width": "220px"
	});
    
	$("#leftminimized").css({
		"display": "none",
		"position": "absolute", 
		"padding": "5px 8px 5px 8px",
		"background-color": "black", 
		"color": "white",
		"top": "0px",
		"left": "0px"
	});
	$("#leftminimized").html('<input type="button" value="Click to Maximize Left Panel" onclick="toggleleft();"/>');

	$("#rightmenu").css({
		"display": "block",
		"position": "absolute", 
		"padding": "5px 8px 5px 8px",
		"background-color": "black", 
		"color": "white",
		"top": "0px",
		"right": "0px",
		"bottom": "0px",
		"overflow-y": "auto",
		//"overflow-x": "auto",
		"max-width": "220px" //'30%'; 270 before//
	});	
		
	$("#rightminimized").css({
		"display": "none",
		"position": "absolute", 
		"padding": "5px 8px 5px 8px",
		"background-color": "black", 
		"color": "white",
		"top": "0px",
		"right": "0px",
		"overflow-y": "auto"
	});	
	$("#rightminimized").html('<input type="button" value="Click to Maximize Right Panel" onclick="toggleright();"/>');

	$("#searchDiv").css({
		"display": "block",
		"padding": "0px 5px 5px 5px",
		"background-color": "black", 
		"color": "white",
		"overflow-y": "auto",
		"max-width": "270px",
		"max-height": "300px"
	});
	document.getElementById("rightmenu").appendChild(searchDiv);
	document.getElementById("rightmenu").appendChild(separator);
		
	$("#fcgrDiv").css({
		"display": "block",
		"padding": "5px 5px 5px 5px",
		"background-color": "black", 
		"color": "white",
		"max-width": "230px"
	});


	var cgrTitle = document.createElement("font");
	cgrTitle.style.color = "yellow";
	cgrTitle.appendChild(document.createTextNode("CGR image"));
	var cgrTitle2 = document.createElement("strong");
	var cgrTitle3 = document.createElement("em");
	cgrTitle2.appendChild(cgrTitle);
	cgrTitle3.appendChild(cgrTitle2);
	cgrInfoDiv.appendChild(cgrTitle3);
	cgrInfoDiv.appendChild(fcgrDiv);
	document.getElementById("rightmenu").appendChild(cgrInfoDiv);
	cgrInfoDiv.style.display='none';
	
	var colorTitle = document.createElement("font");
	colorTitle.style.color = "yellow";
	colorTitle.appendChild(document.createTextNode("Pick Highlight Color"));
	var colorTitle2 = document.createElement("strong");
	var colorTitle3 = document.createElement("em");
	colorTitle2.appendChild(colorTitle);
	colorTitle3.appendChild(colorTitle2);
	document.getElementById("pickcolor").insertBefore(colorTitle3,document.getElementById("pickbtn"));
	document.getElementById("pickcolor").appendChild(separator);
	document.getElementById("rightmenu").appendChild(document.getElementById("pickcolor"));
	
	document.getElementById("rightmenu").appendChild(distPointsDiv);

	$("#pointInfoDiv").css({
		"display": "block",
		"padding": "5px 5px 5px 5px",
		"background-color": "black", 
		"color": "white",
		"max-width": "250px"
	});

	document.getElementById("rightmenu").appendChild(pointInfoDiv);


	// INFO DIV
	staticInfoDiv = '<input type="button" value="Click to Minimize this Panel" onclick="toggleleft();"/><br>\
	<em><strong><font color="yellow" size="4">How to Navigate</font></strong></em>\
	<br>Left | Right = A | D\
	<br>Up | Down = W | S\
	<br>Zoom In | Out = E | Q\
	<br>Rotate = Left Click & Drag\
	<br><hr color="white" width="60%"><em><strong><font color="yellow" size="4">Advanced Settings</font></strong></em>\
	<div id="advsettings" style="display:block"><table>\
	<tr><td>Dim1</td><td><select id="dim1"><option value="0">SELECT</option><option value="1">1</option><option value="2">2</option><option value="3">3</option><option value="4">4</option><option value="5">5</option></select></td></tr>\
	<tr><td>Dim2</td><td><select id="dim2"><option value="0">SELECT</option><option value="1">1</option><option value="2">2</option><option value="3">3</option><option value="4">4</option><option value="5">5</option></select></td></tr>\
	<tr><td>Dim3</td><td><select id="dim3"><option value="0">SELECT</option><option value="1">1</option><option value="2">2</option><option value="3">3</option><option value="4">4</option><option value="5">5</option></select></td></tr>\
	<tr><td>Radius</td><td><input type="number" id="radius" min="0.005" max="0.1" step="0.005" value='+radius+'></td></tr>\
	</table>\
	<input type="submit" value="Re-Draw!" onclick="redraw();">\
	<br>\
	<br>Show CGR image: <input type="checkbox" onchange="toggle(\'cgrInfoDiv\');" unchecked>\
	<br>Show highlight color: <input type="checkbox" onchange="toggle(\'pickcolor\');" unchecked>\
	<br>Show distances: <input type="checkbox" onchange="toggle(\'distPointsDiv\');" unchecked>\
	</div>';
		
	staticInfoDiv += '<hr color="white" width="60%"><em><strong><font color="yellow" size="4">Legend</font></strong></em><br>'+mapCaption+'<br><br>';
	//<hr color="white" width="60%">';   

	for (var styleInd=0; styleInd<legendColors.length; styleInd++) {
		staticInfoDiv = staticInfoDiv + '<div>'+'<img src="img/'+legendColors[styleInd]+'.png" height="20" style="border:0px solid white"> '+legendLabels[styleInd].charAt(0).toUpperCase()+legendLabels[styleInd].slice(1)+'</div>';
	}
	$("#leftmenu").html(staticInfoDiv);

	// SETUP SEARCH DIV
	searchDiv.innerHTML='<input type="button" value="Click to Minimize this Panel" onclick="toggleright();"/>\
	<br><em><strong><font color="yellow" size="4">Search</font></strong></em> \
	<br><input id="tosearch" type="text" onkeyup="startSearch();" onfocus="startSearch();" value="" maxlength="15" >\
	<div id="searchstatus"></div><hr color="white" width="60%">';

	// DISTPOINTS DIV
	distPointsDiv.innerHTML='<em><strong><font color="yellow" size="4">Distance between Points</font></strong></em><br>\
	<input type="button" value="Load Map\'s Distance Matrix" onclick="computeDist(true);">\
	<div id="computeDist"></div>\
	<hr color="white" width="60%">';

	// SET DIMENSIONS IN DROPDOWN MENU
	document.getElementById("dim1").selectedIndex=dim1;
	document.getElementById("dim2").selectedIndex=dim2;
	document.getElementById("dim3").selectedIndex=dim3;

	updateInfoDiv();

	function getDistanceSquare(p1, p2) {
		if (!p1 || !p2) return null;
		var dx = p1.x - p2.x;
		var dy = p1.y - p2.y;
		var dz = p1.z - p2.z;
		return dx * dx + dy * dy + dz * dz;
	}

	var canvas = document.getElementsByTagName('canvas')[0];
	var downx = 0, downxConst = 0;
	var downy = 0, downyConst = 0;
	mouseIsDown = false;

	canvas.onmousedown = function(e){
		if (e.which != 1) return;
		downx = downxConst = e.clientX;
		downy = downyConst = e.clientY;
		mouseIsDown = true;
	}
		
	canvas.onmouseup = function(e){
		if (e.which != 1) return;
		mouseIsDown = false;
		mouse.x = ( e.clientX / canvas.width ) * 2 - 1;
		mouse.y = - ( e.clientY / canvas.height) * 2 + 1;
		
		var vector = new THREE.Vector3(mouse.x, mouse.y, 1);
		vector.unproject(camera);
		raycaster.set(camera.position, vector.sub(camera.position).normalize());
		var intersects = raycaster.intersectObjects(scene.children);
				
		if(dbg){console.log("Found="+intersects.length+" intersections..");}

		if (intersects.length > 0 && downxConst == e.clientX && downyConst == e.clientY) {
			allHits = new Array();
			minDistance=1000;
			var minIndex;
			for(var myInd=0; myInd <intersects.length; myInd++){
				intersectedPoint= intersects[myInd].point;
				var op = [intersectedPoint.x, intersectedPoint.y, intersectedPoint.z];
				var im = new THREE.Matrix4();
				im.getInverse(meshes[0].matrixWorld);
				im.applyToVector3Array(op);
				op = {x:op[0], y:op[1], z:op[2]};
				if(dbg){console.log(op);}

				// loop through all of them to find the point closer to the ray
				//minDistance = getDistanceSquare(op, offsets[0]);
				//var minIndex = 0;
				for (var i = 0; i < offsets.length; i ++) {
					var dist = getDistanceSquare(op, offsets[i]);
					if (dist < minDistance) {
						minDistance = dist;
						minIndex = i;
					}
				}
				
				allHits[myInd]=[];
				allHits[myInd].push(minIndex);
				allHits[myInd].push(minDistance);
			}
			
			if(minDistance<3){
				for(var hitInd=0; hitInd<allHits.length; hitInd++){
					if(100*(allHits[hitInd][1]-minDistance)<1){
						minIndex=hitInd;
						break;
					}  
				}
				//console.log("hitInd="+hitInd);
				selectedIndex=allHits[minIndex][0];
				op = {x:offsets[selectedIndex].x, y:offsets[selectedIndex].y, z:offsets[selectedIndex].z};
				console.log("FINAL Closest point is= "+selectedIndex);
				
				selectedPoint = op;
				//// for disabling selection 
				updateSelectedMesh();
				updateInfoDiv();      
			}         
		}
	}

	canvas.onmousemove = function(e){
		if(!mouseIsDown) return;
		var dx = e.clientX - downx;
		var dy = e.clientY - downy;
		downx = e.clientX;
		downy = e.clientY;
		
		for(var i1=0; i1<meshes.length;i1++){
			meshes[i1].rotation.y += 0.005 * dx;
			while (meshes[i1].rotation.y > 2 * Math.PI){meshes[i1].rotation.y -= 2 * Math.PI;}
			while (meshes[i1].rotation.y < -2 * Math.PI){meshes[i1].rotation.y += 2 * Math.PI;}
			if ( (dy > 0 && meshes[i1].rotation.x < Math.PI/2) || (dy < 0 && meshes[i1].rotation.x > -Math.PI/2) ) {
				meshes[i1].rotation.x += 0.005 * dy;
			}
			selectedMesh.rotation.x = meshes[i1].rotation.x;
			selectedMesh.rotation.y = meshes[i1].rotation.y;
		}
		if(dbg){console.log('rotX='+meshes[0].rotation.x+'|rotY='+meshes[0].rotation.y);}   
		return false;
	}
	
	canvas.onmouseleave = function(e){
		mouseIsDown = false;
	}
		
	document.onkeydown = function(e) {
		if(dbg){console.log("keypressed="+e.which)};
		
		// toggle left menu with ctrl+alt+m
		if (e.ctrlKey && e.altKey && e.which == 77) {toggleright();}
		// toggle right menu with ctrl+alt+n
		if (e.ctrlKey && e.altKey && e.which == 78) {toggleleft();}

		if(document.activeElement.id!="tosearch"){
			var d =0.2;
			// ctrl + alt + z 90 => go to search field
			if (e.ctrlKey && e.altKey && e.which == 90) {$("#tosearch").focus();}
			// ctrl + alt + r 82 => reset camera position, keep mesh rotation
			if (e.ctrlKey && e.altKey && e.which == 82) {
				camera.position.x=0;
				camera.position.y=0;
				camera.position.z=15;
			}			
			// +- numpad  or eq
			if (e.which == 109 || e.which == 81) {camera.position.z += d;}
			if (e.which == 107 || e.which == 69) {camera.position.z -= d;}
			//wasd 87 65 83 68
			if (e.which == 87) {camera.position.y -= d;}
			if (e.which == 65) {camera.position.x += d;}
			if (e.which == 83) {camera.position.y += d;}
			if (e.which == 68) {camera.position.x -= d;}
			if(dbg){console.log('x='+camera.position.x+'|y='+camera.position.y+'|z='+camera.position.z);}
		}
	} 

	if(dbg){console.log("endGraphics!");}

	if(dbg){console.log("check for wrong coordinates!");}
	for(var indtmp=0;indtmp<globalScaledPointsCoord.length;indtmp++){
		if(Math.abs(globalScaledPointsCoord[indtmp][0])>3.5){
			console.log(indtmp+" = "+globalScaledPointsCoord[indtmp][0]+" WHAAAAAAAAAAAAATTTT??");
		}
	}

	if(dbg){console.log("autonavigation mode is ON!");}
		
	if(geturlparamvalue("search")!=-1){
		document.getElementById("tosearch").value=decodeURIComponent(geturlparamvalue("search"));
		if(dbg){console.log("we got filter value! Now we start the search!");}
		startSearch();
	}


	if(geturlparamvalue('autonavigate')=='true'){
		camera.position.x = 0;
		camera.position.y = 0;
		camera.position.z = 10;
		dxTimer = setInterval(function(){myCameraMove('x',cameraX)},50);
		dyTimer = setInterval(function(){myCameraMove('y',cameraY)},50);
		dzTimer = setInterval(function(){myCameraMove('z',cameraZ)},50);
		dxRotTimer = setInterval(function(){myCameraMove('xRot',meshesRotX)},50);
		dyRotTimer = setInterval(function(){myCameraMove('yRot',meshesRotY)},50);       
	}
}

// MAP CONTENT PARSING - STEP 2
function mapContentParsing() {
	if(dbg){console.log("begin mapContentParsing()..");}
	
	var res=alldata.split("\n");
	setOfPoints=res[0].split(','); 
	colors=res[1].split(',');  
	numberOfLabels=parseFloat(res[2]);
	namesOfLabels=res[3].split(',');
	legendColors=res[4].split(',');
	legendLabels=res[5].split(',');
	mapCaption=res[6];

	var pointDimensions=5;
	var startInd=7;
	if(res[7]=="pointDimensions"){
		pointDimensions=parseInt(res[8]);
		startInd=9;
	}
	console.log("Each Point's Dimensions="+pointDimensions);
	alldata="";  //free memory

	allAccessionNums = [];
	allNames = [];
	allBioInfo = [];
	globalScaledPointsCoord = [];
	globalPointsLabels = [];
	var tmpCoordsAllDims,centerX,centerY,centerZ;
	var accInd = namesOfLabels.indexOf('Acc');
	var nameInd = namesOfLabels.indexOf('Name');
	var taxaInd = namesOfLabels.indexOf('Taxa');
	var bioInd = namesOfLabels.indexOf('BioInfo');
	console.log("accInd="+accInd);
	console.log("nameInd="+nameInd);
	console.log("taxaInd="+taxaInd);
	console.log("bioInd="+bioInd);
	var pt_counter=-1;
	for(var ind=startInd; ind<res.length;ind=ind+pointDimensions+numberOfLabels){
		pt_counter++;
		tmpCoordsAllDims=res.slice(ind,ind+pointDimensions); 
		centerX=parseFloat(tmpCoordsAllDims[dim1-1]);
		centerY=parseFloat(tmpCoordsAllDims[dim2-1]);
		centerZ=parseFloat(tmpCoordsAllDims[dim3-1]);
		globalScaledPointsCoord[pt_counter]=[centerX,centerY,centerZ]; 

		globalPointsLabels[pt_counter]=[];
		for(var comm_ind=0;comm_ind<numberOfLabels;comm_ind++){
			globalPointsLabels[pt_counter].push(res[ind+pointDimensions+comm_ind]);
			// POPULATE allAccessionNums AND allNames BY TAKING THE RIGHT INDEX EACH TIME
			if(comm_ind==accInd){allAccessionNums.push(res[ind+pointDimensions+comm_ind]);}
			if(comm_ind==nameInd){allNames.push(res[ind+pointDimensions+comm_ind]);}
			if(comm_ind==taxaInd){allBioInfo.push(res[ind+pointDimensions+comm_ind]);}
			if(comm_ind==bioInd){allBioInfo.push(res[ind+pointDimensions+comm_ind]);}
		}
	}

	if(dbg){
		console.log("setofpoints="+setOfPoints);
		console.log("colors="+colors);
		console.log("#labels="+numberOfLabels);
		console.log("namesLabels="+namesOfLabels);
		console.log("legend_Colors="+legendColors);
		console.log("legend_labels="+legendLabels);
		console.log("globalPointsLabels="+globalPointsLabels);
	}

	if(dbg){console.log("end mapContentParsing()..");}
	if(dbg){console.log("to call initGraphics()..");}
	initGraphics();
}

// READ PARAMETER VALUES + LOAD MAP FILE - STEP 1
$(document).ready(function(){
	
	$("#cgrInfoDiv").hide();
	$("#distPointsDiv").hide();

	if(geturlparamvalue('dbg')=='true'){dbg=true;}    
	mapid=geturlparamvalue('mapid');
	dim1=geturlparamvalue('dim1');
	dim2=geturlparamvalue('dim2');
	dim3=geturlparamvalue('dim3');
	if((dim1<0)||(dim1>5)){dim1=1;}else{dim1=parseInt(dim1);}
	if((dim2<0)||(dim2>5)){dim2=2;}else{dim2=parseInt(dim2);}
	if((dim3<0)||(dim3>5)){dim3=3;}else{dim3=parseInt(dim3);}
	radius=geturlparamvalue('radius');
	if((radius<=0)||(radius>0.1)){radius=0.04;}
	
	if(geturlparamvalue('autonavigate')=='true'){
		meshesRotX = parseFloat(geturlparamvalue('meshesRotX'));
		meshesRotY = parseFloat(geturlparamvalue('meshesRotY'));
		cameraX = parseFloat(geturlparamvalue('cameraX'));
		cameraY = parseFloat(geturlparamvalue('cameraY'));
		cameraZ = parseFloat(geturlparamvalue('cameraZ'));
	}else{
		meshesRotX = 0;
		meshesRotY = 0;
		cameraX = 0;
		cameraY = 0;
		cameraZ = 15;
	}
	
	if(dbg){console.log("x="+cameraX+"|y="+cameraY+"|z="+cameraZ+"|rotX="+meshesRotX+"|rotY="+meshesRotY);}
	if(dbg){console.log("radius="+radius);}
	if(dbg){console.log("mapid="+mapid);}
	if(mapid.length>0){
		var tmpMapId = geturlparamvalue('mapid');
		if(tmpMapId.slice(0,5)=='local'){
			alldata = localStorage.getItem("mapfile"+tmpMapId.slice(5));
			distMatrix = localStorage.getItem("distMatrix"+tmpMapId.slice(5)).split(',');
			console.log(distMatrix.length);
			//console.log(distMatrix);
			//console.log(alldata);
			if(alldata!=null){
				mapContentParsing();
			}else{
				document.write('<font color="red" size="4"><b>The mapid you provided ['+mapid+'] cannot be found. <br>Available in localStorage ['+Object.keys(localStorage)+'].</b></font>');
			}
		}else{
			$.ajax({
				type: 'GET',
				url: "maps/"+mapid, 
				error: function(xhr, status, error) {
				  console.log('STATUS=',status,'ERROR=',error);
				  document.write('<font color="red" size="4"><b>The mapid you provided ['+mapid+'] has status ['+status+'] and error ['+error+'].</b></font>');
				},
				success: function(result){
					alldata=result;
					if(dbg){console.log("to call mapContentParsing()..");}
					mapContentParsing();
				}
			});	
		}
	}else{
		document.write('<br><center><font color="red" size="4"><b>Please specify a map to load from <a href="index.html">here</a></b></font></center>');
	}
});

