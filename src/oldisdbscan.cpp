
//calc KNNDists full dataset in memory
void calcKnnDists( Rcpp::NumericMatrix &data,  Rcpp::NumericMatrix &knnDists,  Rcpp::IntegerMatrix &knnIdx, int &k, int &numRows, int &numCols){


    Rcpp::NumericVector currentTmpRow(numCols);
    Rcpp::NumericVector tmpRow(numCols);
    Rcpp::IntegerVector pointMaxIdx(numRows); // da passare per riferimento nella versione a chunk
    Rcpp::NumericVector pointMaxValue(numRows);

    for (int i = 0; i < numRows; i++){
        pointMaxIdx(i) = -1;
        pointMaxValue(i) = 99999;
    }

    double total, diff;
    int i,j,m,n;
    for( i = 0; i < numRows; i++){

        currentTmpRow = data.row(i);

        for(  j = 0; j < numRows; j++){ //cosi calcolo solo metà
            tmpRow = data.row(j);
            total = 0;

            if ( i != j){
                for(  m = 0; m < numCols ;m++){
                    diff = tmpRow(m) - currentTmpRow(m);
                    total += diff * diff;

                }
                pointMaxValue(i) = knnDists(i,0);
                pointMaxIdx(i) = 0;
                for (  n = 0;  n < k; n++){
                    if (knnDists(i,n) == 0){
                        pointMaxIdx(i) = n;
                        pointMaxValue(i) = knnDists(i,n);
                        break;
                    }
                    else if ( pointMaxValue(i) < knnDists(i,n)    ){
                        pointMaxIdx(i) = n;
                        pointMaxValue(i) = knnDists(i,n);
                    }
                }
                if ( total < pointMaxValue(i) || pointMaxValue(i) == 0){
                    pointMaxValue(i) = total;
                    knnDists(i,pointMaxIdx(i)) = total;
                    knnIdx(i,pointMaxIdx(i)) = j;
                }
            }
        }
    }

}

//calculate KnnDist Chunk
void calcKnnDists(Rcpp::NumericVector &currentRow ,Rcpp::NumericMatrix &data,  Rcpp::NumericMatrix &knnDists,  Rcpp::IntegerMatrix &knnIdx, Rcpp::IntegerVector &pointMaxIdx , Rcpp::NumericVector &pointMaxValue  ,int &k, int &numRows, int &numCols, int &offset, int &rowsToRead, int &distFor){

    Rcpp::NumericVector currentTmpRow(numCols);
    Rcpp::NumericVector tmpRow(numCols);

    double total, diff;

    for( int j = 0; j < rowsToRead; j++){ //cosi calcolo solo metà
        tmpRow = data.row(j);
        total = 0;

        if ( distFor != (j+offset)){
            for( int m = 0; m < numCols ;m++){
                diff = tmpRow(m) - currentRow(m);
                total += diff * diff;

            }
            pointMaxValue(distFor) = knnDists(distFor,0);
            pointMaxIdx(distFor) = 0;
            for ( int n = 0;  n < k; n++){
                if (knnDists(distFor,n) == 0){
                    pointMaxIdx(distFor) = n;
                    pointMaxValue(distFor) = knnDists(distFor,n);
                    break;
                }
                else if ( pointMaxValue(distFor) < knnDists(distFor,n)    ){
                    pointMaxIdx(distFor) = n;
                    pointMaxValue(distFor) = knnDists(distFor,n);
                }
            }
            if ( total < pointMaxValue(distFor) || pointMaxValue(distFor) == 0){
                pointMaxValue(distFor) = total;
                knnDists(distFor,pointMaxIdx(distFor)) = total;
                knnIdx(distFor,pointMaxIdx(distFor)) = j+offset;
            }
        }
    }
}

//calcolo ISK
Rcpp::IntegerMatrix findIsk(Rcpp::IntegerMatrix &knnIdx, int &k, int &numRows, Rcpp::IntegerVector &nOfIsk){

    Rcpp::IntegerMatrix iskNeighbours(numRows,k);
    Rcpp::IntegerVector currentRow(k);
    Rcpp::IntegerVector tmpRow(k);


    int neighbour_id;

    for (int i = 0; i < numRows; i++){
        currentRow = knnIdx.row(i);
        //iskCurrentRow = iskNeighbours.row(i);
        Rcpp::IntegerVector iskCurrentRow;
        for(int j = 0; j < k; j++){
            neighbour_id = currentRow(j);
            tmpRow = knnIdx.row(neighbour_id);

            for(int m=0; m < k; m++ ){
                if (i == tmpRow(m)){
                    iskCurrentRow.push_back(neighbour_id);
                }
            }
        }
        nOfIsk(i) = iskCurrentRow.size();
        if (nOfIsk(i) < k){
            for(int m = nOfIsk(i); m < k; m++){
                iskCurrentRow.push_back(NA_INTEGER);
            }
            iskNeighbours.row(i) = iskCurrentRow;
        }
        else{
            iskNeighbours.row(i) = iskCurrentRow;
        }

    }
    return iskNeighbours;
}




//BEGIN STRATIFICATION


void calculateINFLO(Rcpp::NumericMatrix &knnDists, Rcpp::IntegerMatrix &iskMatrix, Rcpp::IntegerVector &nOfIsk, Rcpp::NumericVector &INFLOVector, int &k, int &numRows,  double &max_inflo){

    double sum;
    double current_knn_max = -1;
    double knn_max = -1;
    Rcpp::IntegerVector inf_inflo;
    int tmp;
    for (int i = 0; i < numRows; i++){

        current_knn_max = knnDists(i,k-1);

        if(nOfIsk(i) > 0 ){
            sum = 0.0;
            for (int j = 0; j < nOfIsk(i); j++){
                tmp = iskMatrix(i,j);
                knn_max = knnDists(tmp,k-1);
                sum += 1/knn_max;
            }

            sum = sum / nOfIsk(i);
            sum = sum * current_knn_max;
            INFLOVector(i) = sum;

            if (sum > max_inflo){
                max_inflo = sum;
            }
        }
        else{
            inf_inflo.push_back(i);
            INFLOVector(i) = std::numeric_limits<double>::infinity();
        }
        current_knn_max = -1;
    }
    for (auto &idxPoint : inf_inflo){
        INFLOVector(idxPoint) = max_inflo;
    }
}

//calculateINFLOCorrect
void calculateINFLOCorrect( Rcpp::IntegerVector &nOfRNN,Rcpp::IntegerVector &nOfIsk, Rcpp::NumericVector &INFLOVector, int &numRows){
    for(int i = 0;  i < numRows; i++){
        if(nOfIsk(i) > 0){
            INFLOVector(i) =  INFLOVector(i) / nOfRNN(i);
        }
    }
}

//calculate normalized INFLO
void calculateINFLONorm(Rcpp::NumericVector &INFLOVector, int &numRows,  double &max_inflo,double &max_NNdist ){
    for(int i = 0;  i < numRows; i++){
        INFLOVector(i) = (INFLOVector(i) / max_inflo ) * max_NNdist;
    }
}

//create partition
void createPartition( Rcpp::NumericVector &INFLOPlusNNdist, Rcpp::NumericVector &INFLOVector, Rcpp::IntegerVector &map,  Rcpp::IntegerVector &layerVector ,double &mean_current_INFLOPlusNNdist, int& cut, int &layer, double &INFLO_threshold, bool& stop, int &numRows){
    int new_cut = 0;

    for(int i = 0;  i < numRows; i++ ){
        if(INFLOPlusNNdist(i) < mean_current_INFLOPlusNNdist){
            new_cut = i;
        }
    }

    double interval_mean_INFLO = 0;
    if ( 1+cut >= numRows ){
        for (int i = cut-1; i < new_cut + 1; i++ ){
            interval_mean_INFLO = interval_mean_INFLO + INFLOVector(map(i));
        }
        interval_mean_INFLO = interval_mean_INFLO / (new_cut - cut +1);
    }
    else{
        for (int i = cut; i < new_cut + 1; i++ ){
            interval_mean_INFLO = interval_mean_INFLO + INFLOVector(map(i));
        }
        interval_mean_INFLO = interval_mean_INFLO / (new_cut - cut +1);
    }

    for (int i = cut; i < new_cut +1; i++){
        layerVector(map(i)) = layer;
    }

    if (interval_mean_INFLO < INFLO_threshold){
        cut =  new_cut +1;
        stop = false;
    }
    else{
        cut = new_cut + 1;
        stop = true;
    }
}


//stratifier
void stratifier(Rcpp::NumericVector &INFLOVector, Rcpp::NumericVector &NNdistVector, Rcpp::NumericVector &INFLOPlusNNdist, Rcpp::IntegerVector &map, Rcpp::IntegerVector &newOrder, Rcpp::IntegerVector &layerVector, int &numRows){

    double mean_INFLO = std::accumulate(INFLOVector.begin(), INFLOVector.end(), 0.0)/ numRows;
    double sum = std::accumulate(INFLOVector.begin(), INFLOVector.end(), 0.0);
    double mean = sum/ numRows;

    Rcpp::NumericVector diff(numRows);
    std::transform(INFLOVector.begin(), INFLOVector.end(), diff.begin(), [mean](double x) { return x - mean; });

    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double var_INFLO = (sq_sum / numRows);

    int layer =  0;
    int cut = 0;

    Rcpp::IntegerVector cuts;
    double INFLO_treshold =mean_INFLO + var_INFLO;
    //begin partitioning
    while(size_t(cut) < numRows){
        double mean_current_INFLOPlusNNdist = std::accumulate( INFLOPlusNNdist.begin() + cut , INFLOPlusNNdist.end(), 0.0) / (numRows - cut);
        bool stop;


        createPartition(INFLOPlusNNdist,INFLOVector,map,layerVector,mean_current_INFLOPlusNNdist,cut, layer,INFLO_treshold, stop, numRows);

        layer++;
        cuts.push_back(cut);

        if(stop){
            for(int i = cut; i < numRows; i++){
                layerVector(map(i)) = layer;
            }
            cut = numRows;
            cuts.push_back(cut);
        }
    }
}


Rcpp::IntegerVector calcCluster(Rcpp::IntegerMatrix &iskNeighbours, int &numRows, int &k,  Rcpp::IntegerVector &nOfIsk, Rcpp::LogicalVector &border){

    Rcpp::IntegerVector allClusters(numRows);
    Rcpp::IntegerVector needCheck(k);

    for(int i = 0; i < numRows; i++){
        allClusters(i) = -1;
    }

    bool start = true;
    int cluster = 0;

    int not_assigned = numRows;
    int pos_init = 0;

    while (not_assigned > 0){
        int i = 0;
        bool found = false;

        while ( i < numRows && !found && !start ){
            if(allClusters(i) == -1){
                pos_init = i;
                found = true;
            }
            else{
                i++;
            }
        }

        int temp_not_assigned = not_assigned;
        Rcpp::IntegerVector tmp_list;
        tmp_list.push_back(pos_init);
        allClusters(pos_init) = cluster;
        not_assigned--;

        expandCluster(iskNeighbours, allClusters,  tmp_list, pos_init, cluster, k,not_assigned, nOfIsk, border );

        if (temp_not_assigned - not_assigned >= k){
            cluster++;
        }
        else{
            for (auto &idxPoint : tmp_list){
                if (idxPoint != pos_init){
                    if(nOfIsk(idxPoint) == 0){
                        allClusters(idxPoint) = -1;
                        not_assigned++;
                    }
                }
            }
            allClusters(pos_init) = -2;
        }
        start  = false;
    }
    return allClusters;
}


//cluster post stratification
Rcpp::IntegerVector calcCluster(Rcpp::IntegerMatrix &iskNeighbours, int &numRows, int &k,  Rcpp::IntegerVector &nOfIsk ,Rcpp::LogicalVector &border, Rcpp::IntegerVector &map){
    Rcpp::IntegerVector allClusters(numRows);
    Rcpp::IntegerVector needCheck(k);
    for(int i = 0; i < numRows; i++){
        allClusters(i) = -1;
    }

    bool start = true;
    int cluster = 0;

    int not_assigned = numRows;
    int pos_init = 0;
    while (not_assigned > 0){
        int i = 0;
        bool found = false;

        while ( i < numRows && !found && !start ){
            if(allClusters(map(i)) == -1){
                pos_init = map(i);
                found = true;
            }
            else{
                i++;
            }
        }

        int temp_not_assigned = not_assigned;
        Rcpp::IntegerVector tmp_list;
        tmp_list.push_back(pos_init);
        allClusters(pos_init) = cluster;
        not_assigned--;

        expandCluster(iskNeighbours, allClusters,  tmp_list, pos_init, cluster, k,not_assigned, nOfIsk, border);

        if (temp_not_assigned - not_assigned >= k){
            cluster++;
        }
        else{
            for (auto &idxPoint : tmp_list){
                if (idxPoint != pos_init){
                    if(nOfIsk(idxPoint) == 0){
                        allClusters(idxPoint) = -1;
                        not_assigned++;
                    }
                }
            }
            allClusters(pos_init) = -2;
        }
        start  = false;
    }
    return allClusters;

}


