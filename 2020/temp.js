// this.namesets = namesets;
        // this.allSequences = [];
        // this.error = false;
        // this.errorReason = '';
        // this.calls2NCBI = 0;
        // this.totalNumOfSeq = 0;
        // this.totalNumOfSeqEachGroup = [];
        // var noNameSet = false;
        // if(this.namesets.length == 0){ noNameSet = true; }
        // for(var i=0; i<this.accession_numbersSets.length; i++){
        //     this.totalNumOfSeqEachGroup.push(this.accession_numbersSets[i].length);
        //     if(noNameSet ){ this.namesets.push("NA"); }
        //     for(var j=0; j<this.accession_numbersSets[i].length; j++){
        //         this.totalNumOfSeq += 1;
        //         this.allSequences.push(["", "", "NA", "", "" ]);
        //         this.accession_numbers.push(this.accession_numbersSets[i][j]);
        //     }
        // }
        // this.distMatrix = [];
        // this.finalEigenvec = [];
        // this.finalEigenval = [];
        // this.finalPoints = [];
        // this.labelColors = ['blue', 'red', 'green', 'orange', 'magenta', 'yellow', 'brown', 'lime', 'gray', 'pink','cyan','seagreen', 'olive', 'purple', 'black'];
        // this.loadSequencesStart = '';
        // this.loadSequencesEnd = '';
        // this.clusters_found = [];
        // this.clustering_assignment = [];
        // this.clusters_composition = [];


        			// if(ui_element_for_progress_update != null){
			// 	let progress = index*100.0/accession_numbers.length;
			// 	progress = progress.toFixed(2);
			// 	$('#' + ui_element_for_progress_update).html(progress+'%');
            // }
            

            loadFastaFromNCBICallback (ind, accID, output, header, localDBG = false){
                this.allSequences[ind][0] = accID;
                this.allSequences[ind][1] = output.length;
                this.allSequences[ind][3] = buildFCGR(output, this.kmer, accID, localDBG);
                this.allSequences[ind][4] = header;
                if(this.taxa_info && !runLocally){ this.calls2NCBI += 0.5; }else{ this.calls2NCBI += 1; }
                if(localDBG){console.log("loadFastaCB: ",[ind,accID,output.length,'DONE']);}
                if(this.calls2NCBI%100 == 0){ console.log('Completed '+this.calls2NCBI+' of '+this.totalNumOfSeq+' calls ('+(100*this.calls2NCBI/this.totalNumOfSeq).toFixed(2)+'%)'); }
                if(this.calls2NCBI == this.totalNumOfSeq && !this.error){ 
                    this.loadSequencesEnd = + new Date();
                    console.log("loadSequences for ["+this.id+"] COMPLETED in ["+((this.loadSequencesEnd-this.loadSequencesStart)/1000)+"] sec");
                    this.computeDistMatrix(localDBG); 
                }
                // console.log("getFastaCB-NCBI: ",accID, this.calls2NCBI, this.totalNumOfSeq);
            }	
        
            loadFastaFromNCBI(ind, accID, callback, callbackObj, localDBG = false){
                var urlToLoad;
                // console.log("runLocally=",runLocally);
                if(runLocally){ 
                    urlToLoad = "http://localhost/haplomaps3d/allfiles/"+accID+".fasta"; 
                }else{
                    urlToLoad = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="+accID;
                }
                $.ajax({
                    url: urlToLoad,
                    success: function(result){
                        var data = result.split(/\n\r?/gi);
                        var header;
                        if( data.length && data[0][0]==='>' ){ header = data[0]; }
                        while (data.length && data[0][0] === '>') {data.shift();}
                        var outputFasta=data.join('');
                        if(localDBG){console.log('loadFasta: ['+accID+'] len= ['+outputFasta.length+']');}
                        callback.apply(callbackObj, [ind, accID, outputFasta, header, localDBG]);
                    },
                    error: function(xhr, status, error) {
                        this.error = true;
                        this.errorReason = 'On "loadFastaFromNCBI" ['+accID+'], Status= ['+status+'], Error= ['+error+']';
                    }
                });
            }
        
            getTaxaFromNCBICallback (ind, accID, taxa, localDBG = false){
                this.allSequences[ind][2] = taxa;
                if(this.taxa_info){ this.calls2NCBI += 0.5; }
                if(localDBG){console.log("getTaxaCB: ["+accID+"]");}
                if(this.calls2NCBI%100 == 0){ console.log('Completed '+this.calls2NCBI+' of '+this.totalNumOfSeq+' calls ('+(100*this.calls2NCBI/this.totalNumOfSeq).toFixed(2)+'%)'); }
                if(this.calls2NCBI == this.totalNumOfSeq && !this.error){ 
                    this.loadSequencesEnd = + new Date();
                    console.log("loadSequences for ["+this.id+"] COMPLETED in ["+((this.loadSequencesEnd-this.loadSequencesStart)/1000)+"] sec");
                    this.computeDistMatrix(localDBG); 
                }
                // console.log("getTaxaCB-NCBI: ",accID, this.calls2NCBI, this.totalNumOfSeq);
            }
        
            getTaxaFromNCBI(ind, accID, callback, callbackObj, localDBG = false){
                $.ajax({
                    url: "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gb&retmode=xml&id="+accID,
                    success: function(result){
                        var taxaString = new XMLSerializer().serializeToString(result);
                        var taxa = /<GBSeq_taxonomy>(.*)<\/GBSeq_taxonomy>/g.exec(taxaString);
                        if(localDBG){console.log("getTaxa: ["+accID+"] taxa: ["+taxa[1]+"]");}
                        callback.apply(callbackObj, [ind, accID, taxa[1].trim(), localDBG]);
                    },
                    error: function(xhr, status, error) {
                        this.error = true;
                        this.errorReason = 'On "getTaxaFromNCBI" ['+accID+'], Status= ['+status+'], Error= ['+error+']';
                    }
                });		
            }
        
            loadSequences(localDBG = false){
                console.log("loadSequences for ["+this.id+"]");
                if(this.accession_numbers.length < 5){
                    this.error = true;
                    this.errorReason = 'Found '+this.accession_numbers.length+' sequences. Please enter at least 5 sequences and try again.';
                }
                this.loadSequencesStart = + new Date();
                for(var i=0; i<this.accession_numbers.length; i++){
                    if(!this.error){
                        this.loadFastaFromNCBI(i, this.accession_numbers[i].trim(), this.loadFastaFromNCBICallback, this, localDBG);
                        if(this.taxa_info && !runLocally){ this.getTaxaFromNCBI(i, this.accession_numbers[i].trim(), this.getTaxaFromNCBICallback, this, localDBG); }
                    }
                }
            }


            // if(dbg){console.log(distMatrix);}
			// var toprint="{";
			// for(var i=0; i<distMatrix.length; i++){
			// 	toprint+='{';
			// 	for (var j=0; j<distMatrix.length; j++){
			// 		toprint+=distMatrix[i][j]+',';
			// 	}
			// 	toprint = toprint.substring(0, toprint.length - 1);
			// 	toprint+='},'
			// }
			// toprint = toprint.substring(0, toprint.length - 1);
			// toprint+='}';
			// console.log('COPY/PASTE MATHEMATICA');
			// console.log("Eigensystem[",toprint,"]");

            time = new Date();
			step2B = time.getTime();		
			$("#stepstatus").html('\
				<p><strong>Step 1 (NCBI + FCGRs) DONE in ~'+(step1B-step1A)/1000+' sec</p>\
				<p>Step 2 (AID) DONE in ~'+(step2B-step2A)/1000+' sec</p></strong>\
				<p>Step 3 of 3: (MDS)</p>');
			$('#progress').html('Please wait.. (browser may freeze for a while) <img src="img/loading.gif" width="30" height="30"/>');
            
            

            // function assembleChunks(data){
		// 	for(let i=0; i<data.length; i++){
		// 		let tmpInfo = data[i];
		// 		// console.log("tmpInfo=",i,tmpInfo);
		// 		let bgX = tmpInfo[0], bgY = tmpInfo[2];
		// 		let endX = tmpInfo[1], endY = tmpInfo[3];
		// 		for(let indexRow=0; indexRow<= endX - bgX; indexRow++){
		// 			for(let indexCol=0; indexCol<= endY -bgY; indexCol++){
		// 				// console.log(bgX + indexRow, bgY + indexCol, tmpInfo[4][indexRow][indexCol],distMatrix[bgX + indexRow][bgY + indexCol]);
		// 				distMatrix[bgX + indexRow][bgY + indexCol] = tmpInfo[4][indexRow][indexCol];
		// 				distMatrix[bgY + indexCol][bgX + indexRow] = tmpInfo[4][indexRow][indexCol];
		// 			}
		// 		}
        //     }
        // }
        

        // computeDistMatrix(localDBG = false){
    //     console.log("computeDistMatrix for ["+this.id+"]");	
    //     var t1 = + new Date(), t2;

    //     var dim = this.totalNumOfSeq;
    //     for(var i=0; i<dim; i++){
    //         this.distMatrix.push([]);
    //         for(var j=0; j<dim; j++){
    //             this.distMatrix[i].push(0);
    //         }
    //     }
        
    //     var splitBy = Math.min( Math.max(Math.ceil(dim/5), 50) , 100);
    //     if(localDBG){ console.log("splitBy= ",splitBy); }
    //     var input = []
    //     for(var i=1; i<= Math.ceil(dim/splitBy); i++){
    //         for(var j=i; j<= Math.ceil(dim/splitBy); j++){
    //             var minX, maxX, minY, maxY;
    //             minX = (i-1)*splitBy;
    //             maxX = Math.min(i*splitBy-1,dim-1);
    //             minY = (j-1)*splitBy;
    //             maxY = Math.min(j*splitBy-1,dim-1);
    //             var allSequencesX = [] , allSequencesY = [];
    //             for(var i00=minX; i00<=maxX; i00++){
    //                 allSequencesX.push(this.allSequences[i00]);
    //             }
    //             for(var i00=minY; i00<=maxY; i00++){
    //                 allSequencesY.push(this.allSequences[i00]);
    //             }
    //             input.push([minX, maxX, minY, maxY, allSequencesX, allSequencesY, this.kmer]);
    //         }
    //     }
    //     console.log("input= ",input);
    //     console.log('Parallel computation of '+input.length+' matrices of ~ '+Math.min(splitBy,this.totalNumOfSeq)+'x'+Math.min(splitBy,this.totalNumOfSeq)+' each. Please wait..');

    //     var p = new Parallel( input );
        
    //     var computeChunks = function (chunk) {
    //         var bgX = chunk[0], endX = chunk[1], bgY = chunk[2], endY = chunk[3];
    //         var allSequencesX = chunk[4], allSequencesY = chunk[5], fcgrRes = chunk[6];
    //         var res=[], tmpRow, numerator, denominator ;
    //         var timeBegin = new Date();
    //         // console.log("Worker ["+bgX+","+endX+" - "+bgY+","+endY+"]");
            
    //         for(var i = 0; i<allSequencesX.length; i++){
    //             tmpRow=[];
    //             for(var j = 0; j<allSequencesY.length; j++){
    //                 numerator = allSequencesX[i][3]['size'] + allSequencesY[j][3]['size'];
    //                 denominator=0;
    //                 for(var unionInd=0; unionInd<Math.pow(4,fcgrRes); unionInd++){
    //                     if( (allSequencesX[i][3]['flatFCGR'][unionInd]==1) || (allSequencesY[j][3]['flatFCGR'][unionInd]==1) ){
    //                         denominator++;
    //                     } 
    //                 }
    //                 // console.log(i,j,'=',numerator,denominator,2 - numerator*1.0/denominator);
    //                 tmpRow.push(2 - numerator*1.0/denominator);
    //             }
    //             res.push(tmpRow);
    //         }
            
    //         var timeEnd = new Date();
    //         console.log("Worker ["+bgX+","+endX+" - "+bgY+","+endY+"] in ["+String((timeEnd-timeBegin)/1000)+"] sec");
    //         return [bgX, endX, bgY, endY, res];
    //     };

    //     p.map(computeChunks).then(data => {
            
    //         for(var i=0; i<data.length; i++){
    //             var tmpInfo = data[i];
    //             // console.log("tmpInfo=",i,tmpInfo);
    //             var bgX = tmpInfo[0], bgY = tmpInfo[2];
    //             var endX = tmpInfo[1], endY = tmpInfo[3];	
    //             for(var indexRow=0; indexRow<= endX - bgX; indexRow++){
    //                 for(var indexCol=0; indexCol<= endY -bgY; indexCol++){
    //                     this.distMatrix[bgX + indexRow][bgY + indexCol] = tmpInfo[4][indexRow][indexCol];
    //                     this.distMatrix[bgY + indexCol][bgX + indexRow] = tmpInfo[4][indexRow][indexCol];
    //                 }
    //             }

    //         }

    //         if(localDBG){ console.log(this.distMatrix); }
    //         t2 = + new Date();
    //         console.log("computeDistMatrix for ["+this.id+"]  COMPLETED in ["+((t2-t1)/1000)+"] sec");
    //         this.computeMDS(localDBG);
    //       });
    // }





    
	// function computeDistMatrix(sequences_info_list){
	// 	time = new Date();
	// 	step2A = time.getTime();
	// 	console.log("BEGIN computeDistMatrix()");	
	// 	$("#stepstatus").html('\
	// 		<p><strong>Step 1 (NCBI + FCGRs) DONE in ~'+(step1B-step1A)/1000+' sec</strong></p>\
	// 		<p>Step 2 of 3: (AID)</p>');
	// 	$('#progress').html('0%');

	// 	var dim=accIDs.length;
	// 	distMatrix=[];
	// 	for(var i=0; i<dim; i++){
	// 		distMatrix.push([]);
	// 		for(var j=0; j<dim; j++){
	// 			distMatrix[i].push(0);
	// 		}
	// 	}
		
	// 	splitBy = Math.max(Math.ceil(dim/5), 50);
	// 	console.log("splitBy= ",splitBy);
	// 	var input = []
	// 	for(var i=1; i<= Math.ceil(dim/splitBy); i++){
	// 		for(var j=i; j<= Math.ceil(dim/splitBy); j++){
	// 			var minX, maxX, minY, maxY;
	// 			minX = (i-1)*splitBy;
	// 			maxX = Math.min(i*splitBy-1,dim-1);
	// 			minY = (j-1)*splitBy;
	// 			maxY = Math.min(j*splitBy-1,dim-1);
	// 			var allSequencesX = [] , allSequencesY = [];
	// 			for(var i00=minX; i00<=maxX; i00++){
	// 				allSequencesX.push(sequences_info_list[i00]);
	// 			}
	// 			for(var i00=minY; i00<=maxY; i00++){
	// 				allSequencesY.push(sequences_info_list[i00]);
	// 			}
	// 			input.push([minX, maxX, minY, maxY, allSequencesX, allSequencesY, fcgr_resolution]);
	// 		}
	// 	}
	// 	console.log("input= ",input);
	// 	$('#progress').html('Parallel computation of '+input.length+' matrices of ~ '+splitBy+'x'+splitBy+' each. Please wait..');
	// 	var p = new Parallel( input );
		
	// 	var computeChunks = function (chunk) {
	// 		var bgX = chunk[0], endX = chunk[1], bgY = chunk[2], endY = chunk[3];
	// 		var allSequencesX = chunk[4], allSequencesY = chunk[5], fcgrRes = chunk[6];
	// 		var res=[], tmpRow, numerator, denominator ;
	// 		var timeBegin = new Date();
	// 		console.log("Worker ["+bgX+","+endX+" - "+bgY+","+endY+"]");

	// 		for(var i = 0; i<allSequencesX.length; i++){
	// 			tmpRow=[];
	// 			for(var j = 0; j<allSequencesY.length; j++){
	// 				numerator = allSequencesX[i]['fcgr']['size'] + allSequencesY[j]['fcgr']['size'];
	// 				denominator=0;
	// 				for(var unionInd=0; unionInd<Math.pow(4,fcgrRes); unionInd++){
	// 					if( (allSequencesX[i]['fcgr']['flatFCGR'][unionInd]==1) || (allSequencesY[j]['fcgr']['flatFCGR'][unionInd]==1) ){
	// 						denominator++;
	// 					} 
	// 				}
	// 				// console.log(i,j,'=',numerator,denominator,2 - numerator*1.0/denominator);
	// 				tmpRow.push(2 - numerator*1.0/denominator);
	// 			}
	// 			res.push(tmpRow);
	// 		}

	// 		var timeEnd = new Date();
	// 		console.log("Worker ["+bgX+","+endX+" - "+bgY+","+endY+"] execution time: "+String((timeEnd-timeBegin)/1000)+" sec");
	// 		return [bgX, endX, bgY, endY, res];
	// 	};
 

	// 	function assembleChunks(data){
	// 		// console.log("data=",data,data.length);
	// 		// console.log("dMat=",distMatrix.length);

	// 		for(var i=0; i<data.length; i++){
	// 			var tmpInfo = data[i];
	// 			// console.log("tmpInfo=",i,tmpInfo);
	// 			var bgX = tmpInfo[0], bgY = tmpInfo[2];
	// 			var endX = tmpInfo[1], endY = tmpInfo[3];
	// 			for(var indexRow=0; indexRow<= endX - bgX; indexRow++){
	// 				for(var indexCol=0; indexCol<= endY -bgY; indexCol++){
	// 					// console.log(bgX + indexRow, bgY + indexCol, tmpInfo[4][indexRow][indexCol],distMatrix[bgX + indexRow][bgY + indexCol]);
	// 					distMatrix[bgX + indexRow][bgY + indexCol] = tmpInfo[4][indexRow][indexCol];
	// 					distMatrix[bgY + indexCol][bgX + indexRow] = tmpInfo[4][indexRow][indexCol];
						
	// 				}
	// 			}
	// 		}	
			
	// 		// if(dbg){console.log(distMatrix);}
	// 		// var toprint="{";
	// 		// for(var i=0; i<distMatrix.length; i++){
	// 		// 	toprint+='{';
	// 		// 	for (var j=0; j<distMatrix.length; j++){
	// 		// 		toprint+=distMatrix[i][j]+',';
	// 		// 	}
	// 		// 	toprint = toprint.substring(0, toprint.length - 1);
	// 		// 	toprint+='},'
	// 		// }
	// 		// toprint = toprint.substring(0, toprint.length - 1);
	// 		// toprint+='}';
	// 		// console.log('COPY/PASTE MATHEMATICA');
	// 		// console.log("Eigensystem[",toprint,"]");


	// 		time = new Date();
	// 		step2B = time.getTime();		
	// 		$("#stepstatus").html('\
	// 			<p><strong>Step 1 (NCBI + FCGRs) DONE in ~'+(step1B-step1A)/1000+' sec</p>\
	// 			<p>Step 2 (AID) DONE in ~'+(step2B-step2A)/1000+' sec</p></strong>\
	// 			<p>Step 3 of 3: (MDS)</p>');
	// 		$('#progress').html('Please wait.. (browser may freeze for a while) <img src="img/loading.gif" width="30" height="30"/>');
			

	// 	}

	// 	return p.map(computeChunks).then(assembleChunks);		
	// }

	// function buildFCGR(seq, k){
	// 	var numOfKmers=0, curmer, curmerInd, kmers={}, mapACGT={"A":0, "C":1, "G":2, "T":3};
	// 	var newseq = seq.replace(/B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|U|V|W|X|Y|Z/gi, "");
	// 	var kmersInds = new Uint8Array(Math.pow(4, k));
	// 	if(dbg){console.log((seq.length - newseq.length)+' lost bp');}
	// 	for(var i=0; i<newseq.length-k+1; i++){
	// 		curmer = newseq.slice(i,i+k);
	// 		curmerInd = 0;
	// 		for(var j=0; j<curmer.length; j++){ 
	// 			curmerInd += Math.pow(4, k-1-j)*mapACGT[curmer[j]]; 
	// 		}
			
	// 		if(kmersInds[curmerInd]==0){
	// 			numOfKmers++;
	// 			kmersInds[curmerInd]=1;
	// 		}

	// 	}
	// 	if(dbg){console.log(numOfKmers+' #kmers '+kmersInds.length);}
	// 	return {"size":numOfKmers, "flatFCGR":kmersInds};
	// }

	// async function fetch_fasta_from_ncbi(idx, accession_number){
	// 	// let fasta_seq = fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="+accession_number, {method: 'GET'})
	// 	let fasta_seq = fetch("http://localhost:8081/2020/fasta/"+accession_number , {method: 'GET'})
	// 	.then(response => {
	// 		// console.log(response);
	// 		const reader = response.body.getReader();
	// 		return new ReadableStream({
	// 			start(controller) {
	// 			return pump();
	// 			function pump() {
	// 				return reader.read().then(({ done, value }) => {
	// 				// When no more data needs to be consumed, close the stream
	// 				if (done) {
	// 					controller.close();
	// 					return;
	// 				}
	// 				// Enqueue the next data chunk into our target stream
	// 				controller.enqueue(value);
	// 				return pump();
	// 				});
	// 			}
	// 			}  
	// 		})
	// 	})
	// 	.then(stream => new Response(stream))
	// 	.then(response => response.text())
	// 	.then(function(result){
	// 		let data = result.split(/\n\r?/gi);
	// 		let header;
	// 		if( data.length && data[0][0]==='>' ){ header = data[0]; }
	// 		while (data.length && data[0][0] === '>') {data.shift();}
	// 		let fasta_contents=data.join('');
	// 		console.info('[' + idx + '] [' + accession_number + '] len= ' + fasta_contents.length + ' bp');
			
	// 		let seq_info = {"index": idx, 
	// 						"accession_number": accession_number, 
	// 						"sequence_length": fasta_contents.length, 
	// 						"taxa_info": null, 
	// 						"fcgr": fasta_contents,   // this will be updated with fcgr later
	// 						"header": header};
	// 		return seq_info;
	// 	})
	// 	.catch(err => console.error(err));

	// 	return fasta_seq;
	// }

	// async function fetch_taxa_from_ncbi(idx, accession_number){
	// 	let taxa_info = fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gb&retmode=xml&id="+accession_number, {method: 'GET'})
	// 	.then(response => {
	// 		// console.log(response);
	// 		const reader = response.body.getReader();
	// 		return new ReadableStream({
	// 			start(controller) {
	// 			return pump();
	// 			function pump() {
	// 				return reader.read().then(({ done, value }) => {
	// 				// When no more data needs to be consumed, close the stream
	// 				if (done) {
	// 					controller.close();
	// 					return;
	// 				}
	// 				// Enqueue the next data chunk into our target stream
	// 				controller.enqueue(value);
	// 				return pump();
	// 				});
	// 			}
	// 			}  
	// 		})
	// 	})
	// 	.then(stream => new Response(stream))
	// 	.then(response => response.text())
	// 	.then(function(response) {
	// 		let taxa = /<GBSeq_taxonomy>(.*)<\/GBSeq_taxonomy>/g.exec(response);
	// 		console.debug("taxa=", taxa[1]);
	// 		return taxa[1];
	// 	});

	// 	return taxa_info;
	// }

	// async function fetch_sequences_and_build_fcgrs(accession_ids_list, enable_taxa_info=false, fcgr_resolution, ui_element_for_progress_update=null){
	// 	console.debug('Start fetch_sequences');
	// 	let all_sequences = [];

	// 	for (let index = 0; index < accession_ids_list.length; index++) {
	// 		// const tmp =  await new Promise(resolve => setTimeout(resolve, 1200));
	// 		// console.log(index);
			
	// 		let cur_accession = accession_ids_list[index];

	// 		let seq_info = await fetch_fasta_from_ncbi(index, cur_accession);
	// 		console.info(seq_info);

	// 		if(enable_taxa_info){
	// 			let taxa_info = await fetch_taxa_from_ncbi(index, cur_accession)
	// 			console.info(taxa_info);
	// 			seq_info["taxa_info"] = taxa_info;
	// 		}

	// 		seq_info["fcgr"] = buildFCGR(seq_info["fcgr"], fcgr_resolution);
	// 		all_sequences.push(seq_info);

	// 		if(ui_element_for_progress_update != null){
	// 			let progress = index*100.0/accession_ids_list.length;
	// 			progress = progress.toFixed(2);
	// 			$('#' + ui_element_for_progress_update).html(progress+'%');
	// 		}
	// 	}

	// 	console.debug('End fetch_sequences');
	// 	return all_sequences;
	// }
	// function computeDistMatrix(sequences_info_list){
	// 	time = new Date();
	// 	step2A = time.getTime();
	// 	console.log("BEGIN computeDistMatrix()");	
	// 	$("#stepstatus").html('\
	// 		<p><strong>Step 1 (NCBI + FCGRs) DONE in ~'+(step1B-step1A)/1000+' sec</strong></p>\
	// 		<p>Step 2 of 3: (AID)</p>');
	// 	$('#progress').html('0%');

	// 	var dim=accIDs.length;
	// 	distMatrix=[];
	// 	for(var i=0; i<dim; i++){
	// 		distMatrix.push([]);
	// 		for(var j=0; j<dim; j++){
	// 			distMatrix[i].push(0);
	// 		}
	// 	}
		
	// 	splitBy = Math.max(Math.ceil(dim/5), 50);
	// 	console.log("splitBy= ",splitBy);
	// 	var input = []
	// 	for(var i=1; i<= Math.ceil(dim/splitBy); i++){
	// 		for(var j=i; j<= Math.ceil(dim/splitBy); j++){
	// 			var minX, maxX, minY, maxY;
	// 			minX = (i-1)*splitBy;
	// 			maxX = Math.min(i*splitBy-1,dim-1);
	// 			minY = (j-1)*splitBy;
	// 			maxY = Math.min(j*splitBy-1,dim-1);
	// 			var allSequencesX = [] , allSequencesY = [];
	// 			for(var i00=minX; i00<=maxX; i00++){
	// 				allSequencesX.push(sequences_info_list[i00]);
	// 			}
	// 			for(var i00=minY; i00<=maxY; i00++){
	// 				allSequencesY.push(sequences_info_list[i00]);
	// 			}
	// 			input.push([minX, maxX, minY, maxY, allSequencesX, allSequencesY, fcgr_resolution]);
	// 		}
	// 	}
	// 	console.log("input= ",input);
	// 	$('#progress').html('Parallel computation of '+input.length+' matrices of ~ '+splitBy+'x'+splitBy+' each. Please wait..');
	// 	var p = new Parallel( input );
		
	// 	var computeChunks = function (chunk) {
	// 		var bgX = chunk[0], endX = chunk[1], bgY = chunk[2], endY = chunk[3];
	// 		var allSequencesX = chunk[4], allSequencesY = chunk[5], fcgrRes = chunk[6];
	// 		var res=[], tmpRow, numerator, denominator ;
	// 		var timeBegin = new Date();
	// 		console.log("Worker ["+bgX+","+endX+" - "+bgY+","+endY+"]");

	// 		for(var i = 0; i<allSequencesX.length; i++){
	// 			tmpRow=[];
	// 			for(var j = 0; j<allSequencesY.length; j++){
	// 				numerator = allSequencesX[i]['fcgr']['size'] + allSequencesY[j]['fcgr']['size'];
	// 				denominator=0;
	// 				for(var unionInd=0; unionInd<Math.pow(4,fcgrRes); unionInd++){
	// 					if( (allSequencesX[i]['fcgr']['flatFCGR'][unionInd]==1) || (allSequencesY[j]['fcgr']['flatFCGR'][unionInd]==1) ){
	// 						denominator++;
	// 					} 
	// 				}
	// 				// console.log(i,j,'=',numerator,denominator,2 - numerator*1.0/denominator);
	// 				tmpRow.push(2 - numerator*1.0/denominator);
	// 			}
	// 			res.push(tmpRow);
	// 		}

	// 		var timeEnd = new Date();
	// 		console.log("Worker ["+bgX+","+endX+" - "+bgY+","+endY+"] execution time: "+String((timeEnd-timeBegin)/1000)+" sec");
	// 		return [bgX, endX, bgY, endY, res];
	// 	};
 

	// 	function assembleChunks(data){
	// 		// console.log("data=",data,data.length);
	// 		// console.log("dMat=",distMatrix.length);

	// 		for(var i=0; i<data.length; i++){
	// 			var tmpInfo = data[i];
	// 			// console.log("tmpInfo=",i,tmpInfo);
	// 			var bgX = tmpInfo[0], bgY = tmpInfo[2];
	// 			var endX = tmpInfo[1], endY = tmpInfo[3];
	// 			for(var indexRow=0; indexRow<= endX - bgX; indexRow++){
	// 				for(var indexCol=0; indexCol<= endY -bgY; indexCol++){
	// 					// console.log(bgX + indexRow, bgY + indexCol, tmpInfo[4][indexRow][indexCol],distMatrix[bgX + indexRow][bgY + indexCol]);
	// 					distMatrix[bgX + indexRow][bgY + indexCol] = tmpInfo[4][indexRow][indexCol];
	// 					distMatrix[bgY + indexCol][bgX + indexRow] = tmpInfo[4][indexRow][indexCol];
						
	// 				}
	// 			}
	// 		}	
			
	// 		// if(dbg){console.log(distMatrix);}
	// 		// var toprint="{";
	// 		// for(var i=0; i<distMatrix.length; i++){
	// 		// 	toprint+='{';
	// 		// 	for (var j=0; j<distMatrix.length; j++){
	// 		// 		toprint+=distMatrix[i][j]+',';
	// 		// 	}
	// 		// 	toprint = toprint.substring(0, toprint.length - 1);
	// 		// 	toprint+='},'
	// 		// }
	// 		// toprint = toprint.substring(0, toprint.length - 1);
	// 		// toprint+='}';
	// 		// console.log('COPY/PASTE MATHEMATICA');
	// 		// console.log("Eigensystem[",toprint,"]");


	// 		time = new Date();
	// 		step2B = time.getTime();		
	// 		$("#stepstatus").html('\
	// 			<p><strong>Step 1 (NCBI + FCGRs) DONE in ~'+(step1B-step1A)/1000+' sec</p>\
	// 			<p>Step 2 (AID) DONE in ~'+(step2B-step2A)/1000+' sec</p></strong>\
	// 			<p>Step 3 of 3: (MDS)</p>');
	// 		$('#progress').html('Please wait.. (browser may freeze for a while) <img src="img/loading.gif" width="30" height="30"/>');
			

	// 	}

	// 	return p.map(computeChunks).then(assembleChunks);		
	// }

	// function buildFCGR(seq, k){
	// 	var numOfKmers=0, curmer, curmerInd, kmers={}, mapACGT={"A":0, "C":1, "G":2, "T":3};
	// 	var newseq = seq.replace(/B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|U|V|W|X|Y|Z/gi, "");
	// 	var kmersInds = new Uint8Array(Math.pow(4, k));
	// 	if(dbg){console.log((seq.length - newseq.length)+' lost bp');}
	// 	for(var i=0; i<newseq.length-k+1; i++){
	// 		curmer = newseq.slice(i,i+k);
	// 		curmerInd = 0;
	// 		for(var j=0; j<curmer.length; j++){ 
	// 			curmerInd += Math.pow(4, k-1-j)*mapACGT[curmer[j]]; 
	// 		}
			
	// 		if(kmersInds[curmerInd]==0){
	// 			numOfKmers++;
	// 			kmersInds[curmerInd]=1;
	// 		}

	// 	}
	// 	if(dbg){console.log(numOfKmers+' #kmers '+kmersInds.length);}
	// 	return {"size":numOfKmers, "flatFCGR":kmersInds};
	// }

	// async function fetch_fasta_from_ncbi(idx, accession_number){
	// 	// let fasta_seq = fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="+accession_number, {method: 'GET'})
	// 	let fasta_seq = fetch("http://localhost:8081/2020/fasta/"+accession_number , {method: 'GET'})
	// 	.then(response => {
	// 		// console.log(response);
	// 		const reader = response.body.getReader();
	// 		return new ReadableStream({
	// 			start(controller) {
	// 			return pump();
	// 			function pump() {
	// 				return reader.read().then(({ done, value }) => {
	// 				// When no more data needs to be consumed, close the stream
	// 				if (done) {
	// 					controller.close();
	// 					return;
	// 				}
	// 				// Enqueue the next data chunk into our target stream
	// 				controller.enqueue(value);
	// 				return pump();
	// 				});
	// 			}
	// 			}  
	// 		})
	// 	})
	// 	.then(stream => new Response(stream))
	// 	.then(response => response.text())
	// 	.then(function(result){
	// 		let data = result.split(/\n\r?/gi);
	// 		let header;
	// 		if( data.length && data[0][0]==='>' ){ header = data[0]; }
	// 		while (data.length && data[0][0] === '>') {data.shift();}
	// 		let fasta_contents=data.join('');
	// 		console.info('[' + idx + '] [' + accession_number + '] len= ' + fasta_contents.length + ' bp');
			
	// 		let seq_info = {"index": idx, 
	// 						"accession_number": accession_number, 
	// 						"sequence_length": fasta_contents.length, 
	// 						"taxa_info": null, 
	// 						"fcgr": fasta_contents,   // this will be updated with fcgr later
	// 						"header": header};
	// 		return seq_info;
	// 	})
	// 	.catch(err => console.error(err));

	// 	return fasta_seq;
	// }

	// async function fetch_taxa_from_ncbi(idx, accession_number){
	// 	let taxa_info = fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gb&retmode=xml&id="+accession_number, {method: 'GET'})
	// 	.then(response => {
	// 		// console.log(response);
	// 		const reader = response.body.getReader();
	// 		return new ReadableStream({
	// 			start(controller) {
	// 			return pump();
	// 			function pump() {
	// 				return reader.read().then(({ done, value }) => {
	// 				// When no more data needs to be consumed, close the stream
	// 				if (done) {
	// 					controller.close();
	// 					return;
	// 				}
	// 				// Enqueue the next data chunk into our target stream
	// 				controller.enqueue(value);
	// 				return pump();
	// 				});
	// 			}
	// 			}  
	// 		})
	// 	})
	// 	.then(stream => new Response(stream))
	// 	.then(response => response.text())
	// 	.then(function(response) {
	// 		let taxa = /<GBSeq_taxonomy>(.*)<\/GBSeq_taxonomy>/g.exec(response);
	// 		console.debug("taxa=", taxa[1]);
	// 		return taxa[1];
	// 	});

	// 	return taxa_info;
	// }

	// async function fetch_sequences_and_build_fcgrs(accession_ids_list, enable_taxa_info=false, fcgr_resolution, ui_element_for_progress_update=null){
	// 	console.debug('Start fetch_sequences');
	// 	let all_sequences = [];

	// 	for (let index = 0; index < accession_ids_list.length; index++) {
	// 		// const tmp =  await new Promise(resolve => setTimeout(resolve, 1200));
	// 		// console.log(index);
			
	// 		let cur_accession = accession_ids_list[index];

	// 		let seq_info = await fetch_fasta_from_ncbi(index, cur_accession);
	// 		console.info(seq_info);

	// 		if(enable_taxa_info){
	// 			let taxa_info = await fetch_taxa_from_ncbi(index, cur_accession)
	// 			console.info(taxa_info);
	// 			seq_info["taxa_info"] = taxa_info;
	// 		}

	// 		seq_info["fcgr"] = buildFCGR(seq_info["fcgr"], fcgr_resolution);
	// 		all_sequences.push(seq_info);

	// 		if(ui_element_for_progress_update != null){
	// 			let progress = index*100.0/accession_ids_list.length;
	// 			progress = progress.toFixed(2);
	// 			$('#' + ui_element_for_progress_update).html(progress+'%');
	// 		}
	// 	}

	// 	console.debug('End fetch_sequences');
	// 	return all_sequences;
    // }
    



    endOfMoDMap(localDBG = false){
        console.log(this);
        this.toTextFile(localDBG);
        proceedToNextMap = true;
        if(this.clustering_assignment[this.clustering_assignment.length - 1] == -1){
            proceedToNextMap = false;
            clearInterval(myTimer);
        }
    }

    function computeMDS(distMatrix, allSequences){
		time = new Date();
		step3A = time.getTime();
		console.log("BEGIN computeMDS()");	
		
		// pre-processign MDS

		// part 1
		var numOfSeq = distMatrix.length;
		var distMatrixSqTotals = [];
		for(var i=0; i<numOfSeq; i++){
			var tmpColSum = 0;
			for(var j=0; j<numOfSeq; j++){
				tmpColSum += Math.pow(distMatrix[i][j], 2);
			}
			distMatrixSqTotals.push(tmpColSum);
		}
		if(dbg){console.log('sqtotals=',distMatrixSqTotals);}
		
		// part 2
		var sumOfSqTotals = 0;
		for(var i=0; i<numOfSeq; i++){
			sumOfSqTotals += distMatrixSqTotals[i];
		}
		if(dbg){console.log('sumofsqtotals=',sumOfSqTotals);}
		
		// part 3
		var distMatrixCentered = [];
		for(var i=0; i<numOfSeq; i++){
			distMatrixCentered.push([]);
			for(var j=0; j<numOfSeq; j++){
				var tmpVal = -(1/2)*(Math.pow(distMatrix[i][j], 2) - (1/numOfSeq)*distMatrixSqTotals[i] - (1/numOfSeq)*distMatrixSqTotals[j] + (1/(Math.pow(numOfSeq, 2)))* sumOfSqTotals);
				distMatrixCentered[i].push(tmpVal);
			}
		}
		
		//DEBUGING FOR PRINTING dMatrix
		if(noMDS){
			if(dbg){console.log(distMatrixCentered);}
			
			$('#progress').html('<p>Please continue computation offline. <br> <strong>Instructions:</strong>\
				<ul><li>Download zip file (It is downloaded automatically)</li>\
				<li>Extract zip contents to a directory of your choice</li>\
				<li>Download either <a href="./extra/mds.nb" target="_blank">MDS-Mathematica</a> or <a href="./extra/mds.py" target="_blank">MDS-Python</a> code and place it in the same directory</li>\
				<li>Run the code, you should get an output file [map_output.txt]</li>\
				<li>Open [map_output.txt], copy its content and paste it <a href="./extra/input.html" target="_blank">here</a></li>\
				</ul></p>');
			console.log('COPY/PASTE MATHEMATICA');
			
			var toprint="{";
			for(var i=0; i<distMatrixCentered.length; i++){
				toprint+='{';
				for (var j=0; j<distMatrixCentered.length; j++){
					toprint+=distMatrixCentered[i][j].toFixed(10)+',';
				}
				toprint = toprint.substring(0, toprint.length - 1);
				toprint+='},'
			}
			toprint = toprint.substring(0, toprint.length - 1);
			toprint+='}';
			
			var toprint2='{';
			for(var i=0; i<accIDsSets.length; i++){
				toprint2 += '{';
				for(var j=0; j<accIDsSets[i].length; j++){
					toprint2 += '"'+accIDsSets[i][j]+'",';
				}
				toprint2 = toprint2.substring(0, toprint2.length - 1);
				toprint2 +='},'
			}
			toprint2 = toprint2.substring(0, toprint2.length - 1);
			toprint2 +='}';
			
			var toprint3='{';
			var totalNumOfSeq = Object.keys(allSequences).length;
			for(var i=0; i<totalNumOfSeq; i++){
				toprint3 += '{';
				toprint3 += '"'+allSequences[i][0]+'",';
				toprint3 += '"'+allSequences[i][1]+'",';
				toprint3 += '"'+allSequences[i][2]+'",';
				toprint3 += '"'+allSequences[i][4]+'",';
				toprint3 = toprint3.substring(0, toprint3.length - 1);
				toprint3 +='},'
			}
			toprint3 = toprint3.substring(0, toprint3.length - 1);
			toprint3 +='}';

			var toprint4 = '';
			for(var i=0; i<accIDsSets.length; i++){	toprint4 += accIDsSets[i].length + ",";	}
			toprint4 = toprint4.slice(0,-1);
			toprint4 += "\n";
			for(var i=0; i<accIDsSets.length; i++){ toprint4 += $("#colorset"+i).val() + ","; }
			toprint4 = toprint4.slice(0,-1);
			toprint4 += "\n";
			toprint4 += "5\nIndex,Acc,Name,Length,Taxa\n";
			for(var i=0; i<accIDsSets.length; i++){ toprint4 += $("#colorset"+i).val() + ","; }
			toprint4 = toprint4.slice(0,-1);
			toprint4 += "\n";
			for(var i=0; i<accIDsSets.length; i++){toprint4 +=$("#nameset"+i).val()+" ("+accIDsSets[i].length+"),";}
			toprint4 = toprint4.slice(0,-1);
			toprint4 += "\n";

			for(var i=0; i<$("#numofsets").val(); i++){ toprint4 += $("#nameset"+i).val() + ', '; }
			toprint4 = toprint4.slice(0,-2);  //because of space above!
			// OLD: toprint4 += '. Number of sequences is '+accIDs.length+ "\n";
			toprint4 += '#User-input sequences#'+accIDs.length+"#Approx.Inf.Dist#NA\n";
			
			
			// completely debugging part
			var toprint5="{";
			for(var i=0; i<distMatrix.length; i++){
				toprint5+='{';
				for (var j=0; j<distMatrix.length; j++){
					toprint5+=distMatrix[i][j].toFixed(10)+',';
				}
				toprint5 = toprint5.substring(0, toprint5.length - 1);
				toprint5+='},'
			}
			toprint5 = toprint5.substring(0, toprint5.length - 1);
			toprint5+='}';

			var zip = new JSZip();
			zip.file("dMatrix.txt", toprint);
			zip.file("accIDs.txt", toprint2);
			zip.file("allSequences.txt", toprint3);
			zip.file("mapheader.txt", toprint4);
			zip.file("debug.txt", toprint5);
			var content = zip.generate({type:"blob"});
			// see FileSaver.js
			saveAs(content, "dataForMMA.zip");
			return;
		}
		

		// EIGENVALUES - EIGENVECTORS
		try{
			var res2;
			// part 4 NEW MDS
			var data_type = jsfeat.F64_t | jsfeat.C1_t;
			var columns = numOfSeq, rows = numOfSeq;
			var my_matrix = new jsfeat.matrix_t(columns, rows, data_type);
			my_matrix.data = [].concat.apply([], distMatrixCentered);    
			var my_eigenvec = new jsfeat.matrix_t(columns, rows, data_type);
			var my_eigenval = new jsfeat.matrix_t(columns, rows, data_type);
			jsfeat.linalg.eigenVV(my_matrix, my_eigenvec, my_eigenval);
			
			// var eigenval, eigenvec, eigensystem;
			eigenval = my_eigenval.data;
			eigenvec = [];
			for( var mdsInd = 0; mdsInd < numOfSeq; mdsInd++){
				eigenvec.push(my_eigenvec.data.slice( mdsInd*numOfSeq, (mdsInd+1)*numOfSeq ));
			}
		}catch(err){
			console.log(err.message);
			errorHasOccured = true;
			$('#progress').html('<div class="alert alert-danger" role="alert"><span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span><span class="sr-only">Error:</span>'+' ERROR [<strong>'+err.message+'</strong>] while computing MDS. Internal errors originate from numerical errors in distance matrices. Consider changing/increasing value of length of k-mers.</div>');
			return;
		}


		// try{
		// 	res2 = numeric.eig(distMatrixCentered);	
		// }catch(err){
		// 	console.log(err.message);
		// 	errorHasOccured = true;
		// 	$('#progress').html('<div class="alert alert-danger" role="alert"><span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span><span class="sr-only">Error:</span>'+' ERROR [<strong>'+err.message+'</strong>] while computing MDS. Internal errors originate from numerical errors in distance matrices. Consider changing/increasing value of length of k-mers.</div>');
		// 	return;
		// }
		// if(dbg){console.log('numeric.js=',res2['lambda']['x']);}

		// eigenval = res2['lambda']['x'];
		// eigenvec = numeric.transpose(res2['E']['x']);   /// here transpose!!!

		
		eigensystem =[];
		for(var i=0; i<numOfSeq; i++){
			eigensystem.push([eigenval[i], eigenvec[i]]);
		}
		if(dbg){console.log("eigensystem=",eigensystem);}
		eigensystem.sort(function (a,b){ return Math.abs(b[0]) - Math.abs(a[0]) });
		if(dbg){console.log("eigensystem=",eigensystem);}
		

		// PICK 5 LARGEST AND POSITIVE
		finalEigenvec=[];
		finalEigenval=[];
		for(var i=0; i<numOfSeq; i++){
			if(eigensystem[i][0]>0){
				finalEigenvec.push(eigensystem[i][1]);
				finalEigenval.push(eigensystem[i][0]);
			}
			if(dbg){console.log("i-th eigeinvalue= [",i,"] [",eigensystem[i][0],"]");}
			if(finalEigenvec.length === 5){break;}
		}
		if(dbg){console.log("finalEigenval=",finalEigenval);}
		if(dbg){console.log("finalEigenvec=",finalEigenvec);}
		
		// POINTS COORDINATES
		try{
			points = numeric.dot(numeric.transpose(finalEigenvec), numeric.sqrt(numeric.T.diag(finalEigenval)['x']));
		}catch(err){
			console.log(err.message);
			errorHasOccured = true;
			$('#progress').html('<div class="alert alert-danger" role="alert"><span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span><span class="sr-only">Error:</span>'+' ERROR [<strong>'+err.message+'</strong>] while computing MDS. Internal errors originate from numerical errors in distance matrices. Consider changing/increasing value of length of k-mers.</div>');
			return;
		}
		console.log("points=",points);
		
		// MIN-MAX PER DIMENSION (FOR LATER SCALING)
		var allMinMax =[], tmpMinMax = numeric.transpose(points);
		for(var i=0; i<5; i++){
			allMinMax.push( [ Math.min.apply(Math,tmpMinMax[i]) , Math.max.apply(Math,tmpMinMax[i]) ] );
		}
		if(dbg){console.log(allMinMax);}
		
		
		// PRODUCING MAPFILE
		mapfile = '';
		// // line 1
		// for(var i=0; i<accIDsSets.length; i++){
		// 	mapfile += accIDsSets[i].length + ",";
		// }
		// mapfile = mapfile.slice(0,-1);
		// mapfile += "\n";
		// // line 2
		// for(var i=0; i<accIDsSets.length; i++){
		// 	mapfile += $("#colorset"+i).val() + ",";
		// }
		// mapfile = mapfile.slice(0,-1);
		// mapfile += "\n";
		// // line 3+4
		// mapfile += "5\nIndex,Acc,Name,Length,Taxa\n";
		// // line 5
		// for(var i=0; i<accIDsSets.length; i++){
		// 	mapfile += $("#colorset"+i).val() + ",";
		// }
		// mapfile = mapfile.slice(0,-1);
		// mapfile += "\n";
		// // line 6
		// for(var i=0; i<accIDsSets.length; i++){
		// 	mapfile += $("#nameset"+i).val()+" ("+accIDsSets[i].length+"),";
		// }
		// mapfile = mapfile.slice(0,-1);
		// mapfile += "\n";
		
		// line 7
		// OLD mapfile += "No Description given\n";
		for(var i=0; i<$("#numofsets").val(); i++){
			mapfile += $("#nameset"+i).val() + ', ';
		}
		mapfile = mapfile.slice(0,-2);  //because of space above!
		// old mapfile
		// mapfile += '. Number of sequences is '+accIDs.length+ "\n";
		// mapfile += '#User-input sequences#'+accIDs.length+"#Approx.Inf.Dist#NA\n";

		// rest of the lines
		for(var i=0; i<numOfSeq; i++){
			for(var j=0; j<points[i].length; j++){
				mapfile += ((2*points[i][j] - allMinMax[j][0] - allMinMax[j][1]) / ( allMinMax[j][1] - allMinMax[j][0] ) ) +"\n";	
			}
			mapfile += i+"\n";                      // INDEX
			mapfile += allSequences[i]['accession_number']+"\n";     // ACC
			mapfile += allSequences[i]['header'].split("|").pop().trim()+"\n";     // NAME
			mapfile += allSequences[i]['sequence_length']+"\n";     // LENGTH
			mapfile += allSequences[i]['taxa_info']+"\n";     // TAXA
		}
		mapfile = mapfile.slice(0,-1);
		if(dbg){console.log("c/p=\n",mapfile);}

		var mapTimeStamp = new Date().getTime();
		if(typeof(Storage) !== "undefined") {
			try{
				localStorage.setItem("mapfile"+mapTimeStamp, mapfile);
				console.log(mapfile);
				// localStorage.setItem("distMatrix"+mapTimeStamp, distMatrix);	
			}catch(err){
				console.log(err.message);
				errorHasOccured = true;
				$('#progress').html('<div class="alert alert-danger" role="alert"><span class="glyphicon glyphicon-exclamation-sign" aria-hidden="true"></span><span class="sr-only">Error:</span>'+' ERROR [<strong>'+err.message+'</strong>]. Consider emptying localStorage.</div>');
				return;
			}
		} else {
			console.log("ERROR! NO LOCAL STORAGE AVAILABLE!?");
			alert("ERROR! NO LOCAL STORAGE AVAILABLE!?");
		}
		
		$('#progress').html('100%');
		$("#file_contents").html('<pre>'+mapfile+'</pre>');
		$("#file_contents").hide();
		// END OF MDS
		

		time = new Date();
		step3B = time.getTime();
		$("#stepstatus").html('\
			<p><strong>Step 1 (NCBI + FCGRs) DONE in ~'+(step1B-step1A)/1000+' sec</p>\
			<p>Step 2 (AID) DONE in ~'+(step2B-step2A)/1000+' sec</p>\
			<p>Step 3 (MDS) DONE in ~'+(step3B-step3A)/1000+' sec</p></strong>');
		$("#stepstatus").html($("#stepstatus").html()+'<p><a href="load.html?mapid=local'+mapTimeStamp+'" target="_blank">Show MoDMap</a>');

		$("#stepstatus").html($("#stepstatus").html()+' or <a href="#viewfile" id="viewfile">view contents</a>');

		$("#stepstatus").html($("#stepstatus").html()+'</p>');
		$("#mainprogressbar").hide();
		$("#viewfile").click(function(){
			$("#file_contents").toggle();
		});

		// list local storage for debuging
		for (var a in localStorage) {console.log(a);}
		console.log(Object.keys(localStorage));
	}


    -------------------------


    // ( async() => {
			// var map_data = {
			// 	"description": "No_Description",
			// 	"kmer": 9,
			// 	"enable_taxa_info": false,
			// 	"groups": [
			// 		{
			// 			"name": "Strepsirrhini",
			// 			"color": "green",
			// 			"accession_numbers": ["NC_014453","NC_012761","NC_012762","NC_012763","NC_012764","NC_012766","NC_012769","NC_012771","NC_012773","NC_011053","NC_010299","NC_010300","NC_004025","NC_002765"]
			// 		},
			// 		{
			// 			"name": "Haplorrhini",
			// 			"color": "red",
			// 			"accession_numbers": ["NC_018115","NC_018116","NC_018096","NC_018057","NC_018058","NC_018059","NC_018060","NC_018061","NC_018062","NC_018063","NC_016666","NC_015485","NC_015486","NC_014042","NC_014045","NC_014047","NC_014051","NC_013993","NC_012920","NC_012774","NC_012775","NC_012670","NC_011519","NC_011137","NC_011120","NC_009747","NC_009748","NC_008215","NC_008216","NC_008217","NC_008218","NC_008219","NC_008220","NC_008066","NC_007009","NC_006900","NC_006901","NC_005943","NC_002811","NC_002763","NC_002764","NC_001643","NC_001644","NC_001645","NC_001646","NC_001992","NC_002082","NC_002083"]
			// 		},
			// 	]
			// }
		
		// 	// var curHaploMap = new MoDMap3D(map_data);
		// 	// curHaploMap.print();
		// 	// await curHaploMap.fetch_sequences_and_build_fcgrs();
		// 	// curHaploMap.print();

		// 	// await curHaploMap.compute_distance_matrix();
		// 	// curHaploMap.print();

		// 	// await curHaploMap.multidimensional_scaling();
		// 	// curHaploMap.print();

		// 	// await curHaploMap.toTextFile(true);
		// 	// curHaploMap.print();

		// })();


        -----------------


        				/// convertion utils tests

				// let json_mapfile = cur_modmap.toJSON();
				// console.log("json_mapfile=", json_mapfile);

				// let text_mapfile = new ConvertionUtils().fromJSONToText(json_mapfile);
				// console.log(text_mapfile);

				// let json_mapfile_2 = new ConvertionUtils().fromTextToJSON(text_mapfile);
				// console.log(json_mapfile_2);