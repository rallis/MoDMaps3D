"use strict";

class MoDMap3D {

    constructor (map_data) {
        this.id = + new Date();
        this.groups = map_data["groups"];
        this.kmer = map_data["kmer"];
        this.enable_taxa_info = map_data["enable_taxa_info"];
        this.description = map_data["description"];

        this.accession_numbers = [];
        this.group_names = [];
        this.group_colors = [];
        this.group_sizes = [];
        for(let i = 0; i <this.groups.length; i++){
            let cur_group = this.groups[i];
            this.group_names.push(cur_group.name);
            this.group_colors.push(cur_group.color);
            this.group_sizes.push(cur_group.accession_numbers.length);
            for(let j = 0; j < cur_group.accession_numbers.length; j++){
                this.accession_numbers.push(cur_group.accession_numbers[j]);
            }
        }

        this.all_sequences_info = [];
        this.dist_matrix = [];
        this.finalEigenvec = [];
        this.finalEigenval = [];
        this.finalPoints = [];
    }

    toString (){ 
        return "MoDMap ID= ["+this.id+"], Groups= ["+this.groups.length + "] sizes= ["+this.group_sizes+"] names= ["+this.group_names+"] total= ["+this.accession_numbers.length+"], k= [" + this.kmer + "], Taxa= [" + this.enable_taxa_info+"], description= ["+this.description+"]"; 
    }

    print (){ 
        console.log( this.toString() ); 
        console.log(this.all_sequences_info);
        console.log(this.dist_matrix);
    }


    ///////////////////////
    // FETCHING SEQUENCES
    ///////////////////////

    _buildFCGR(seq, kmer_length){
		let numOfKmers=0, curmer, curmerInd, kmers={}, mapACGT={"A":0, "C":1, "G":2, "T":3};
		let newseq = seq.replace(/B|D|E|F|H|I|J|K|L|M|N|O|P|Q|R|S|U|V|W|X|Y|Z/gi, "");
		let kmersInds = new Uint8Array(Math.pow(4, kmer_length));
		console.log((seq.length - newseq.length)+' lost bp');
		for(let i=0; i<newseq.length-kmer_length+1; i++){
			curmer = newseq.slice(i,i+kmer_length);
			curmerInd = 0;
			for(let j=0; j<curmer.length; j++){ 
				curmerInd += Math.pow(4, kmer_length-1-j)*mapACGT[curmer[j]]; 
			}
			
			if(kmersInds[curmerInd]==0){
				numOfKmers++;
				kmersInds[curmerInd]=1;
			}

		}
		console.log(numOfKmers+' #kmers '+kmersInds.length);
		return {"size":numOfKmers, "flatFCGR":kmersInds};
	}

    async _fetch_fasta_from_file(accession_number, manual_fasta_file){
        console.log(accession_number, manual_fasta_file);

        return new Promise((resolve, reject) => {
            let fr = new FileReader();  
            fr.onload = () => {
                console.log(fr.result);
                
                let data = fr.result.split(/\n\r?/gi);
                let header;
                if( data.length && data[0][0]==='>' ){ header = data[0]; }
                while (data.length && data[0][0] === '>') {data.shift();}
                let fasta_contents=data.join('');
                console.info('[' + accession_number + '] len= ' + fasta_contents.length + ' bp');
			
                let seq_info = {"accession_number": accession_number, 
							"sequence_length": fasta_contents.length, 
							"taxa_info": null, 
							"fcgr": fasta_contents,   // this will be updated with fcgr later
							"header": header};

                resolve(seq_info );
            };
            fr.readAsText(manual_fasta_file);
        });
	}

	async _fetch_fasta_from_ncbi(accession_number){
		// let fasta_seq = fetch("http://localhost:8081/2020/fasta/"+accession_number , {method: 'GET'})
		let fasta_seq = fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id="+accession_number, {method: 'GET'})
		.then(response => {
			// console.log(response);
			const reader = response.body.getReader();
			return new ReadableStream({
				start(controller) {
				return pump();
				function pump() {
					return reader.read().then(({ done, value }) => {
					// When no more data needs to be consumed, close the stream
					if (done) {
						controller.close();
						return;
					}
					// Enqueue the next data chunk into our target stream
					controller.enqueue(value);
					return pump();
					});
				}
				}  
			})
		})
		.then(stream => new Response(stream))
		.then(response => response.text())
		.then(function(result){
			let data = result.split(/\n\r?/gi);
			let header;
			if( data.length && data[0][0]==='>' ){ header = data[0]; }
			while (data.length && data[0][0] === '>') {data.shift();}
			let fasta_contents=data.join('');
			console.info('[' + accession_number + '] len= ' + fasta_contents.length + ' bp');
			
			let seq_info = {"accession_number": accession_number, 
							"sequence_length": fasta_contents.length, 
							"taxa_info": null, 
							"fcgr": fasta_contents,   // this will be updated with fcgr later
							"header": header};
            return seq_info;
		})
		.catch(function(error){
		    console.error(error);
		    return Promise.reject("UNABLE TO FETCH ");	
		});

		return fasta_seq;
	}

	async _fetch_taxa_from_ncbi(accession_number){
		let taxa_info = fetch("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=gb&retmode=xml&id="+accession_number, {method: 'GET'})
		.then(response => {
			// console.log(response);
			const reader = response.body.getReader();
			return new ReadableStream({
				start(controller) {
				return pump();
				function pump() {
					return reader.read().then(({ done, value }) => {
					// When no more data needs to be consumed, close the stream
					if (done) {
						controller.close();
						return;
					}
					// Enqueue the next data chunk into our target stream
					controller.enqueue(value);
					return pump();
					});
				}
				}  
			})
		})
		.then(stream => new Response(stream))
		.then(response => response.text())
		.then(function(response) {
			let taxa = /<GBSeq_taxonomy>(.*)<\/GBSeq_taxonomy>/g.exec(response);
			console.debug("taxa=", taxa[1]);
			return taxa[1];
		});

		return taxa_info;
	}

    async fetch_sequences_and_build_fcgrs(manual_fasta_files, ui_progress_element=null){
        console.debug('Start fetch_sequences');
        let all_sequences = [];
        let manual_fasta_files_idx = 0;
        let progress;

		for (let index = 0; index < this.accession_numbers.length; index++) {
            // artificial delay, because ncbi throttles requests higher than 3 tps
            // 1 tps when fetching fasta only, 2 tps when fetching taxa as well
            // 1200 ms delay since 1000 ms delay may stil get throttled after some iterations
			await new Promise(resolve => setTimeout(resolve, 1200));
			console.log(index);
			
			let cur_accession = this.accession_numbers[index];
            let seq_info;
            if(cur_accession.trim().substr(0,9) == "fastaFile"){
                console.log("Reading from user-uploaded fastaFile..");
                seq_info = await this._fetch_fasta_from_file(cur_accession, manual_fasta_files[manual_fasta_files_idx]);
                manual_fasta_files_idx += 1;
            }else{
                seq_info = await this._fetch_fasta_from_ncbi(cur_accession);
            }
            console.info(seq_info);

			if(this.enable_taxa_info){
				let taxa_info = await this._fetch_taxa_from_ncbi(cur_accession)
				console.info(taxa_info);
				seq_info["taxa_info"] = taxa_info;
			}

			seq_info["fcgr"] = this._buildFCGR(seq_info["fcgr"], this.kmer);
            
            if(seq_info['fcgr']['size'] == 0){
                throw Error("Unable to build FCGR for " + this.accession_numbers[index]);
            }
            all_sequences.push(seq_info);

            if(ui_progress_element != null){
                progress = (index+1)*100.0/this.accession_numbers.length;
                progress = progress.toFixed(2);
                $('#' + ui_progress_element).html(progress+'%');
            }
		}

		console.debug('End fetch_sequences');
        this.all_sequences_info = all_sequences;
        return this.all_sequences_info;
    }


    //////////////////////////
    // COMPUTING DIST MATRIX
    //////////////////////////
    
    async _parallel_compute_dist_matrix_chunks(){
		let dim = this.accession_numbers.length;
		// let splitBy = Math.max(Math.ceil(dim/5), 50);
		let splitBy = 50;  // seems to be optimal in terms of performance
		console.log("splitBy= ",splitBy);
        
        let input = []
		for(let i=1; i<= Math.ceil(dim/splitBy); i++){
			for(let j=i; j<= Math.ceil(dim/splitBy); j++){
				let minX, maxX, minY, maxY;
				minX = (i-1)*splitBy;
				maxX = Math.min(i*splitBy-1,dim-1);
				minY = (j-1)*splitBy;
				maxY = Math.min(j*splitBy-1,dim-1);
				let allSequencesX = [] , allSequencesY = [];
				for(let i00=minX; i00<=maxX; i00++){
					allSequencesX.push(this.all_sequences_info[i00]);
				}
				for(let i00=minY; i00<=maxY; i00++){
					allSequencesY.push(this.all_sequences_info[i00]);
				}
				input.push([minX, maxX, minY, maxY, allSequencesX, allSequencesY, this.kmer]);
			}
		}
		console.log("Parallel Computation input= ",input);
		let p = new Parallel( input );
		
		let computeChunks = function (chunk) {
			let bgX = chunk[0], endX = chunk[1], bgY = chunk[2], endY = chunk[3];
			let allSequencesX = chunk[4], allSequencesY = chunk[5], fcgrRes = chunk[6];
			let res=[], tmpRow, numerator, denominator ;
			let timeBegin = new Date();
			console.log("Worker ["+bgX+","+endX+" - "+bgY+","+endY+"]");

			for(let i = 0; i<allSequencesX.length; i++){
				tmpRow=[];
				for(let j = 0; j<allSequencesY.length; j++){
                    // TODO                    
                    // this is not aware of _compute_AID_distance method, so it re-implements it
                    // would prefer to simply call _compute_AID_distance(allSequencesX[i]['fcgr'], allSequencesY[j]['fcgr'])
                    numerator = allSequencesX[i]['fcgr']['size'] + allSequencesY[j]['fcgr']['size'];
					denominator=0;
					for(let unionInd=0; unionInd<Math.pow(4,fcgrRes); unionInd++){
						if( (allSequencesX[i]['fcgr']['flatFCGR'][unionInd]==1) || (allSequencesY[j]['fcgr']['flatFCGR'][unionInd]==1) ){
							denominator++;
						} 
                    }
					// console.log(i,j,'=',numerator,denominator,2 - numerator*1.0/denominator);
					tmpRow.push(2 - numerator*1.0/denominator);
				}
				res.push(tmpRow);
			}

			let timeEnd = new Date();
			console.log("Worker ["+bgX+","+endX+" - "+bgY+","+endY+"] execution time: "+String((timeEnd-timeBegin)/1000)+" sec");
			return [bgX, endX, bgY, endY, res];
		};

        return p.map(computeChunks);
	}

    _compute_AID_distance(fcgr_1, fcgr_2, kmer_length){
		console.log("custom = ", fcgr_1, fcgr_2, kmer_length);
        
		let aidnumerator = fcgr_1['size'] + fcgr_2['size'];
		let aiddenominator = 0;
		for(let unionInd = 0; unionInd < Math.pow(4,kmer_length); unionInd++){
			if( (fcgr_1['flatFCGR'][unionInd]==1) || (fcgr_2['flatFCGR'][unionInd]==1) ){
				aiddenominator++;
			} 
		}
		console.debug(aidnumerator,aiddenominator,2 - aidnumerator*1.0/aiddenominator);

		return 2 - aidnumerator*1.0/aiddenominator;
	}


    async compute_distance_matrix(){
        let dim = this.accession_numbers.length;
		for(let i=0; i<dim; i++){
			this.dist_matrix.push([]);
			for(let j=0; j<dim; j++){
				this.dist_matrix[i].push(0);
			}
		}
        
        let dist_matrix_chunks = await this._parallel_compute_dist_matrix_chunks();
        
        for(let i=0; i<dist_matrix_chunks.length; i++){
            let tmpInfo = dist_matrix_chunks[i];
            // console.log("tmpInfo=",i,tmpInfo);
            let bgX = tmpInfo[0], bgY = tmpInfo[2];
            let endX = tmpInfo[1], endY = tmpInfo[3];
            for(let indexRow=0; indexRow<= endX - bgX; indexRow++){
                for(let indexCol=0; indexCol<= endY - bgY; indexCol++){
                    // console.log(bgX + indexRow, bgY + indexCol, tmpInfo[4][indexRow][indexCol],distMatrix[bgX + indexRow][bgY + indexCol]);
                    this.dist_matrix[bgX + indexRow][bgY + indexCol] = tmpInfo[4][indexRow][indexCol];
                    this.dist_matrix[bgY + indexCol][bgX + indexRow] = tmpInfo[4][indexRow][indexCol];
                }
            }
        }
        
        // console.log("dist=",dist_matrix);
        return this.dist_matrix;
    }

    
    /////////////////////////////
    // MULTIDIMENSIONAL SCALING
    /////////////////////////////
    
    multidimensional_scaling(){
        console.log("computeMDS for ["+this.id+"]"); 	
        let t1 = + new Date();
        
        // pre-processign MDS
        // part 1
        let numOfSeq = this.accession_numbers.length;
        let distMatrixSqTotals = [];
        for(let i=0; i<numOfSeq; i++){
            let tmpColSum = 0;
            for(let j=0; j<numOfSeq; j++){
                tmpColSum += Math.pow(this.dist_matrix[i][j], 2);
            }
            distMatrixSqTotals.push(tmpColSum);
        }
        console.debug('sqtotals=',distMatrixSqTotals);
        
        // part 2
        let sumOfSqTotals = 0;
        for(let i=0; i<numOfSeq; i++){
            sumOfSqTotals += distMatrixSqTotals[i];
        }
        console.debug('sumofsqtotals=',sumOfSqTotals);
        
        // part 3
        let distMatrixCentered = [];
        for(let i=0; i<numOfSeq; i++){
            distMatrixCentered.push([]);
            for(let j=0; j<numOfSeq; j++){
                let tmpVal = -(1/2)*(Math.pow(this.dist_matrix[i][j], 2) - (1/numOfSeq)*distMatrixSqTotals[i] - (1/numOfSeq)*distMatrixSqTotals[j] + (1/(Math.pow(numOfSeq, 2)))* sumOfSqTotals);
                distMatrixCentered[i].push(tmpVal);
            }
        }

        // part 4 NEW MDS
        let data_type = jsfeat.F64_t | jsfeat.C1_t;
        let columns = numOfSeq, rows = numOfSeq;
        let my_matrix = new jsfeat.matrix_t(columns, rows, data_type);
        my_matrix.data = [].concat.apply([], distMatrixCentered);    
        let my_eigenvec = new jsfeat.matrix_t(columns, rows, data_type);
        let my_eigenval = new jsfeat.matrix_t(columns, rows, data_type);
        jsfeat.linalg.eigenVV(my_matrix, my_eigenvec, my_eigenval);
        
        let eigenval, eigenvec, eigensystem;
        eigenval = my_eigenval.data;
        eigenvec = [];
        for( let mdsInd = 0; mdsInd < numOfSeq; mdsInd++){
            eigenvec.push(my_eigenvec.data.slice( mdsInd*numOfSeq, (mdsInd+1)*numOfSeq ));
        }
        
        if( eigenvec.length == 0  || eigenval.length == 0){
            console.error("ERRROR with Eigenvec or Eigenval", eigenvec, eigenval);
        }
        
        eigensystem =[];
        for(let i=0; i<numOfSeq; i++){
            eigensystem.push([eigenval[i], eigenvec[i]]);
        }
        console.debug("eigensystem=",eigensystem);
        eigensystem.sort(function (a,b){ return Math.abs(b[0]) - Math.abs(a[0]) });
        console.debug("eigensystem=",eigensystem);
    
        // part 5
        // PICK 5 LARGEST AND POSITIVE
        for(let i=0; i<numOfSeq; i++){
            if(eigensystem[i][0]>0){
                this.finalEigenvec.push(eigensystem[i][1]);
                this.finalEigenval.push(eigensystem[i][0]);
            }
            console.debug("i-th eigenvalue= [",i,"] [",eigensystem[i][0],"]");
            if(this.finalEigenvec.length === 5){break;}
        }
        console.debug("finalEigenval=",this.finalEigenval);
        console.debug("finalEigenvec=",this.finalEigenvec);
        
        // part 6
        // POINTS COORDINATES
        let points = numeric.dot(numeric.transpose(this.finalEigenvec), numeric.sqrt(numeric.T.diag(this.finalEigenval)['x']));
        console.debug("points=",points);
        
        // part 7
        // MIN-MAX PER DIMENSION (FOR LATER SCALING)
        let allMinMax =[], tmpMinMax = numeric.transpose(points);
        for(let i=0; i<5; i++){
            allMinMax.push( [ Math.min.apply(Math,tmpMinMax[i]) , Math.max.apply(Math,tmpMinMax[i]) ] );
        }
        console.debug("allMinMax=",allMinMax);

        // part 8
        // SCALE POINTS
        for(let i=0; i<numOfSeq; i++){
            this.finalPoints.push([]);
            for(let j=0; j<points[i].length; j++){
                this.finalPoints[i][j] = ((2*points[i][j] - allMinMax[j][0] - allMinMax[j][1]) / ( allMinMax[j][1] - allMinMax[j][0] ) );
            }
        }
        console.debug("finalPoints=",this.finalPoints);

        let t2 = + new Date();
        console.log("computeMDS for ["+this.id+"] COMPLETED in ["+((t2-t1)/1000)+"] sec");
    }	

    
    /////////////////////////////
    // BACKWORD compatible export
    /////////////////////////////
    
    toJSON(){
        let dataset = {
            "description": this.description,
            "kmer": this.kmer,
            "enable_taxa_info": this.enable_taxa_info,
            "groups": this.groups,
            "points": this.finalPoints,
            "all_sequences_info": this.all_sequences_info,
        }
        return dataset;
    }

    toTextFile(){
        console.log("toTextFile for ["+this.id+"]");

        // PRODUCING MAPFILE
        let mapfile = '';
        // line 1
        for(let i=0; i<this.group_sizes.length; i++){
            mapfile += this.group_sizes[i] + ",";
        }
        mapfile = mapfile.slice(0,-1);
        mapfile += "\n";
        // line 2 
        for(let i=0; i<this.group_colors.length; i++){
            mapfile += this.group_colors[i] + ",";
        }
        mapfile = mapfile.slice(0,-1);
        mapfile += "\n";
        // line 3+4
        mapfile += "5\nIndex,Acc,Name,Length,Taxa\n";
        // line 5
        for(let i=0; i<this.group_colors.length; i++){
            mapfile += this.group_colors[i] + ",";
        }
        mapfile = mapfile.slice(0,-1);
        mapfile += "\n";
        // line 6
        for(let i=0; i<this.group_names.length; i++){
            mapfile += this.group_names[i] +" ("+this.group_sizes[i]+"),";
        }
        mapfile = mapfile.slice(0,-1);
        mapfile += "\n";
        // line 7
        for(let i=0; i<this.group_names.length; i++){
            mapfile += this.group_names[i] + ', ';
        }
        mapfile = mapfile.slice(0,-2);  //because of space above!
        
        // old mapfile
        mapfile += '#'+this.description+'#'+this.accession_numbers.length+"#Approx.Inf.Dist#NA\n";
        // rest of the lines
        for(let i=0; i<this.accession_numbers.length; i++){
            for(let j=0; j<this.finalPoints[i].length; j++){
                // mapfile += ((2*points[i][j] - allMinMax[j][0] - allMinMax[j][1]) / ( allMinMax[j][1] - allMinMax[j][0] ) ) +"\n";	
                mapfile += this.finalPoints[i][j] +"\n";	
            }
            mapfile += i+"\n";                      // INDEX
            mapfile += this.all_sequences_info[i]['accession_number']+"\n";     // ACC
            mapfile += this.all_sequences_info[i]['header'].split("|").pop().trim()+"\n";     // NAME
            mapfile += this.all_sequences_info[i]['sequence_length']+"\n";     // LENGTH
            mapfile += this.all_sequences_info[i]['taxa']+"\n";     // TAXA
        }
        mapfile = mapfile.slice(0,-1);
        // console.log("c/p=\n",mapfile);
        return mapfile;
    }

    toOfflineComputation(){
		console.log("toOfflineComputation for ["+this.id+"]"); 	

        // part 1
        let numOfSeq = this.accession_numbers.length;
        let distMatrixSqTotals = [];
        for(let i=0; i<numOfSeq; i++){
            let tmpColSum = 0;
            for(let j=0; j<numOfSeq; j++){
                tmpColSum += Math.pow(this.dist_matrix[i][j], 2);
            }
            distMatrixSqTotals.push(tmpColSum);
        }
        console.debug('sqtotals=',distMatrixSqTotals);
        
        // part 2
        let sumOfSqTotals = 0;
        for(let i=0; i<numOfSeq; i++){
            sumOfSqTotals += distMatrixSqTotals[i];
        }
        console.debug('sumofsqtotals=',sumOfSqTotals);
        
        // part 3
        let distMatrixCentered = [];
        for(let i=0; i<numOfSeq; i++){
            distMatrixCentered.push([]);
            for(let j=0; j<numOfSeq; j++){
                let tmpVal = -(1/2)*(Math.pow(this.dist_matrix[i][j], 2) - (1/numOfSeq)*distMatrixSqTotals[i] - (1/numOfSeq)*distMatrixSqTotals[j] + (1/(Math.pow(numOfSeq, 2)))* sumOfSqTotals);
                distMatrixCentered[i].push(tmpVal);
            }
        }
        
        let cur_dataset = {
            "description": this.description,
            "kmer": this.kmer,
            "enable_taxa_info": this.enable_taxa_info,
            "groups": this.groups,
        }

        let partial_sequences_info = [];
        for(let i=0; i<this.all_sequences_info.length; i++){
            partial_sequences_info.push({
                "accession_number": this.all_sequences_info[i]['accession_number'],
                "header": this.all_sequences_info[i]["header"],
                "sequence_length": this.all_sequences_info[i]["sequence_length"],
                "taxa_info": this.all_sequences_info[i]["taxa_info"],
            })
        }

        let zip = new JSZip();
        zip.file("dMatrix.json", JSON.stringify(distMatrixCentered));
        zip.file("dataset.json", JSON.stringify(cur_dataset));
        zip.file("all_sequences_info.json", JSON.stringify(partial_sequences_info));
        let content = zip.generate({type:"blob"});
        // see FileSaver.js
        saveAs(content, "dataForMMA.zip");
        return;
    }

    // applyClustering(typeOfClustering = "kmeans", localDBG = false){
    //     console.log("start CLUSTERING",typeOfClustering);
        
    //     for (let i = 0; i < this.accession_numbers.length; i++) { this.clustering_assignment.push(-1); }
        
    //     // select the points' coordinates to apply clustering
    //     let datasetOfPoints = [];
    //     for(let i=0; i < this.finalPoints.length; i++){ datasetOfPoints.push(this.finalPoints[i].slice(0,3)); } 
    //     if(localDBG){ console.log("dataset has been formed"); }

    //     // dbscan
    //     if(typeOfClustering == "dbscan"){
    //         let dbscan = new DBSCAN();
    //         let clusters_dbscan = dbscan.run(datasetOfPoints, 0.2, 5);
    //         let noise_dbscan  = dbscan.noise;
    //         console.log("clusters-DBSCAN=",clusters_dbscan);
    //         console.log("noise-DBSCAN=",noise_dbscan);
    //         if(localDBG){
    //             for(let i=0; i<noise_dbscan.length; i++){
    //                 console.log(this.allSequences[noise_dbscan[i]]);
    //             }	
    //         }
    //         this.clusters_found = clusters_dbscan;	
    //     }
        
    //     // optics
    //     if( typeOfClustering == "optics"){
    //         let optics = new OPTICS();
    //         let clusters_optics = optics.run(datasetOfPoints, 0.2, 5);
    //         console.log("clusters-OPTICS=",clusters_optics);
    //         if(localDBG){
    //             let plot_optics = optics.getReachabilityPlot();
    //             console.log("reachability-plot=",plot_optics);
    //         }
    //         this.clusters_found = clusters_optics;
    //     }

    //     // kmeans
    //     if(typeOfClustering == "kmeans"){
    //         let kmeans = new KMEANS();
    //         let clusters_kmeans = kmeans.run(datasetOfPoints, 3);
    //         let noise_kmeans  = kmeans.noise;
    //         console.log("clusters-KMEANS=",clusters_kmeans);
    //         this.clusters_found = clusters_kmeans;
    //     }

    //     let clusters_size = this.totalNumOfSeqEachGroup;
    //     let clusters_from_clustering = this.clusters_found;			
        
    //     // print and store cluster composition for each cluster returned by clustering method above
    //     let cur_clust;
    //     let cur_clust_res;
    //     for (let i = 0; i < clusters_from_clustering.length; i++) {
    //         this.clusters_composition.push([]);
    //         cur_clust = clusters_from_clustering[i];
    //         cur_clust_res = {"tot": 0};
    //         for (let j = 0; j < cur_clust.length; j++) {
    //             // console.log(cur_clust[j]);
    //             let left_int = 0;
    //             for (let k = 0; k < clusters_size.length; k++) {
    //                 left_int += clusters_size[k];
    //                 if(cur_clust[j] < left_int){
    //                     cur_clust_res["tot"] += 1;
    //                     if( k in cur_clust_res){
    //                         cur_clust_res[k] += 1;
    //                     }else{
    //                         cur_clust_res[k] = 1;
    //                     }
    //                     // console.log(k);
    //                     break;
    //                 }
    //             }		
    //         }
    //         let out = "Tot: "+cur_clust_res["tot"]+" = ";
    //         for (let m = 0; m < clusters_size.length; m++) {
    //             // out += m+": ";
    //             if(cur_clust_res[m]){
    //                 out += cur_clust_res[m]+" ";
    //                 this.clusters_composition[i].push(cur_clust_res[m]);
    //             }else{
    //                 out += "NA ";
    //                 this.clusters_composition[i].push(0);
    //             }
    //         }
    //         if(/*localDBG*/ true){ console.log(out); }
    //     }

    //     // compute clustering assignment
    //     for (let i = 0; i < clusters_from_clustering.length; i++) {
    //         if(localDBG){ console.log("current cluster is = ", clusters_from_clustering[i]); }
            
    //         // detect in which cluster each element of this cluster falls
    //         cur_clust_res = {};	
    //         for (let j = 0; j < clusters_from_clustering[i].length; j++) {
    //             let left_int = 0;
    //             for (let k = 0; k < clusters_size.length; k++) {
    //                 left_int += clusters_size[k];
    //                 if(clusters_from_clustering[i][j] < left_int){
    //                     if( k in cur_clust_res){
    //                         cur_clust_res[k] += 1;
    //                     }else{
    //                         cur_clust_res[k] = 1;
    //                     }
    //                     break;
    //                 }
    //             }
    //         }
    //         if(localDBG){ console.log(cur_clust_res); }

    //         // detect the majority
    //         let final_cluster = -1;
    //         let cur_max = -1;
    //         for (let j = 0; j < this.accession_numbersSets.length; j++) {
    //             if (j in cur_clust_res){
    //                 if(cur_clust_res[j] > cur_max){
    //                     final_cluster = j;
    //                     cur_max = cur_clust_res[j];
    //                 }
    //             }
    //         }
    //         if(localDBG){ console.log(final_cluster); }
        
    //         // assign the correct cluster
    //         for (let j = 0; j < clusters_from_clustering[i].length; j++) {
    //             this.clustering_assignment[ clusters_from_clustering[i][j] ] = final_cluster;
    //         }
    //         if(localDBG){ console.log(this.clustering_assignment); }
    //     }
    //     console.log("clustering_assignment= ",this.clustering_assignment);
    //     console.log("end CLUSTERING");
        
    //     this.endOfMoDMap(localDBG);
    // }


}


class FileUtils {

    fromTextToJSON(text_mapfile) {
        let lines = text_mapfile.split("\n");
        console.log(lines);
        
        let group_points = [];

        let group_sizes = lines[0].split(",");
        let group_colors = lines[1].split(",");
        let num_of_fields = parseInt(lines[2]);
        let fields = lines[3].split(",");

        let parts = lines[6].split("#");
        let group_names = parts[0].split(", ");
        let map_description = parts[1];
        let total_num_of_sequences = parseInt(parts[2]);
        

        let all_sequences_info = [];
        for(let i=0; i<total_num_of_sequences; i++){
            let min_idx = 7 + (5 + num_of_fields) * i;
            let max_idx = 7 + (5 + num_of_fields) * (i + 1) - 1;
            console.log(min_idx, max_idx);
            
            group_points.push([
                parseFloat(lines[min_idx]), 
                parseFloat(lines[min_idx + 1]),
                parseFloat(lines[min_idx + 2]),
                parseFloat(lines[min_idx + 3]),
                parseFloat(lines[min_idx + 4])
            ])
            console.log(group_points);
            
            let data_info = lines.slice(min_idx + 5, max_idx);
            console.log(data_info);
            for(let j=0; j<data_info.length; j++){
                if(fields[j] == "Index"){
                    all_sequences_info.push({
                        "accession_number": -1, 
                        "sequence_length": -1, 
                        "taxa_info": null, 
                        "fcgr": null,
                        "header": null
                    })
                }else if(fields[j] == "Acc"){
                    all_sequences_info[i]["accession_number"] = data_info[j];
                }else if(fields[j] == "Name"){
                    all_sequences_info[i]["header"] = data_info[j];
                }else if(fields[j] == "Length"){
                    all_sequences_info[i]["sequence_length"] = parseInt(data_info[j]);
                }else if(fields[j] == "Taxa"){
                    all_sequences_info[i]["taxa_info"] = data_info[j];
                }else{
                    throw Error("Unknown field in textToJson convertion!");
                }
            }
        }
        console.log(all_sequences_info);

        let all_accession_nums = [];
        for(let i=0; i<all_sequences_info.length; i++) {
            all_accession_nums.push(all_sequences_info[i]["accession_number"]);
        }
        console.log(all_accession_nums);

        let map_groups = [];
        let group_start_idx = 0;
        for (let i= 0; i<group_sizes.length; i++){
            map_groups.push({
                "name": group_names[i], 
                "color": group_colors[i], 
                "accession_numbers": all_accession_nums.slice(group_start_idx, group_start_idx + group_sizes[i]), 
            })
            group_start_idx = group_start_idx + group_sizes[i];
        }

        let dataset = {
            "description": map_description,
            "kmer": 9,
            "enable_taxa_info": all_sequences_info[0]["taxa_info"] != undefined,
            "groups": map_groups,
            "points": group_points,
            "all_sequences_info": all_sequences_info,
        }
        return dataset;

    }

    fromJSONToText(json_mapfile){
        let map_description = json_mapfile['description'];
        let map_groups = json_mapfile['groups'];
        let map_points = json_mapfile['points'];
        let map_all_sequences_info = json_mapfile['all_sequences_info'];
        let map_colors = [];
        let map_names = [];
        let map_group_sizes = [];
        let map_accession_numbers = [];

        for(let i= 0; i<map_groups.length; i++){
            let cur_group = map_groups[i];
            map_colors.push(cur_group['color']);
            map_names.push(cur_group['name']);
            map_group_sizes.push(cur_group['accession_numbers'].length);
            for(let j=0; j<cur_group['accession_numbers'].length; j++){
                map_accession_numbers.push(cur_group['accession_numbers'][j]);
            }
        }

        // PRODUCING MAPFILE
        let mapfile = '';
        // line 1
        for(let i=0; i<map_group_sizes.length; i++){
            mapfile += map_group_sizes[i] + ",";
        }
        mapfile = mapfile.slice(0,-1);
        mapfile += "\n";
        // line 2 
        for(let i=0; i<map_colors.length; i++){
            mapfile += map_colors[i] + ",";
        }
        mapfile = mapfile.slice(0,-1);
        mapfile += "\n";
        // line 3+4
        mapfile += "5\nIndex,Acc,Name,Length,Taxa\n";
        // line 5
        for(let i=0; i<map_colors.length; i++){
            mapfile += map_colors[i] + ",";
        }
        mapfile = mapfile.slice(0,-1);
        mapfile += "\n";
        // line 6
        for(let i=0; i<map_names.length; i++){
            mapfile += map_names[i] +" ("+map_group_sizes[i]+"),";
        }
        mapfile = mapfile.slice(0,-1);
        mapfile += "\n";
        // line 7
        for(let i=0; i<map_names.length; i++){
            mapfile += map_names[i] + ', ';
        }
        mapfile = mapfile.slice(0,-2);  //because of space above!
        
        // old mapfile
        mapfile += '#' + map_description + '#' + map_accession_numbers.length+"#Approx.Inf.Dist#NA\n";
        // rest of the lines
        for(let i=0; i<map_accession_numbers.length; i++){
            for(let j=0; j<map_points[i].length; j++){
                mapfile += map_points[i][j] +"\n";	
            }
            mapfile += i+"\n";                      // INDEX
            mapfile += map_all_sequences_info[i]['accession_number']+"\n";     // ACC
            mapfile += map_all_sequences_info[i]['header'].split("|").pop().trim()+"\n";     // NAME
            mapfile += map_all_sequences_info[i]['sequence_length']+"\n";     // LENGTH
            mapfile += map_all_sequences_info[i]['taxa']+"\n";     // TAXA
        }
        mapfile = mapfile.slice(0,-1);
        // console.log("c/p=\n",mapfile);
        return mapfile;
    }

    async fetch_file(url_to_fetch, is_text = false){
		let fetched_file = fetch(url_to_fetch, {method: 'GET'})
		.then(response => {
			// console.log(response);
			const reader = response.body.getReader();
			return new ReadableStream({
				start(controller) {
				return pump();
				function pump() {
					return reader.read().then(({ done, value }) => {
					// When no more data needs to be consumed, close the stream
					if (done) {
						controller.close();
						return;
					}
					// Enqueue the next data chunk into our target stream
					controller.enqueue(value);
					return pump();
					});
				}
				}  
			})
		})
        .then(stream => new Response(stream))
        .then(function(response){
            if(is_text){
                return response.text();
            }else{
                return response.json();
            }
        })
		.then(function (response) {
            return response;
        })  

		return fetched_file;
	}

    put(key, value){
        if(typeof(Storage) !== "undefined") {
            try{
                localStorage.setItem(key, value);
                return key;
            }catch(err){
                console.error(err.message);
                throw Error('['+err.message+']. Consider emptying localStorage.');
            }
        } else {
            console.error("ERROR! NO LOCAL STORAGE AVAILABLE!?");
            throw Error("ERROR! NO LOCAL STORAGE AVAILABLE!?");
            // list local storage for debuging
            // for (var a in localStorage) {console.log(a);}
            // console.log(Object.keys(localStorage));
            // for (var a in localStorage){ delete localStorage[a]; }
        }
    }

    get(key){
        return localStorage.getItem(key);
    }



}

class UIUtils {
    
    build_alert(alert_html, alert_type){
        if(alert_type === undefined) {
            alert_type = "info";
        }
        return '<div class="alert alert-' + alert_type + '" role="alert">' + alert_html + '</div>';
    }
}