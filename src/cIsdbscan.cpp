#define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::depends(ClusterR)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

#include "beachmat/numeric_matrix.h"
#include "beachmat/integer_matrix.h"
//#include "functions.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <math.h>

using namespace arma;


//get the number of rows
template<typename T>
int get_nrow(const T& data){

    auto matrix_type=beachmat::find_sexp_type(data);

    if(matrix_type== INTSXP){
        auto final_matrix=beachmat::create_integer_matrix(data);
        //const size_t& nc = final_matrix->get_ncol();
        const size_t& nr = final_matrix->get_nrow();
        int n_row = nr;
        return n_row;

    }else if(matrix_type== REALSXP){
        auto final_matrix=beachmat::create_numeric_matrix(data);
        //const size_t& nc = final_matrix->get_ncol();
        const size_t& nr = final_matrix->get_nrow();
        int n_row = nr;
        return n_row;
    }else{

        return 0;
    }
}




//get the number of columns
template<typename T>
int get_ncol(const T& data){

    auto matrix_type=beachmat::find_sexp_type(data);

    if(matrix_type== INTSXP){
        auto final_matrix=beachmat::create_integer_matrix(data);
        const size_t& nc = final_matrix->get_ncol();
        //const size_t& nr = final_matrix->get_nrow();
        int n_col = nc;
        return n_col;

    }else if(matrix_type== REALSXP){
        auto final_matrix=beachmat::create_numeric_matrix(data);
        const size_t& nc = final_matrix->get_ncol();
        //const size_t& nr = final_matrix->get_nrow();
        int n_col = nc;
        return n_col;
    }else{

        return 0;
    }
}

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
    #pragma omp parallel for private(j,m,n,total,diff,currentTmpRow,tmpRow)
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
    #pragma omp barrier

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




//sort knnDists and knnIdx
void sortKnnDists( Rcpp::NumericMatrix &knnDists, Rcpp::IntegerMatrix &knnIdx, int &k, int &numRows){
    bool swapped;
    double swap_tmp;
    int idx_swap_tmp;
    int j,m;

    for(int i = 0; i < numRows; i++){
        for(j = 0; j < k - 1; j++){
            swapped = false;
            for(m = 0; m < k - j - 1; m++){
                if(knnDists(i,m) > knnDists(i,m+1)){
                    swap_tmp = knnDists(i,m);
                    knnDists(i,m) = knnDists(i,m+1);
                    knnDists(i,m+1) = swap_tmp;
                    idx_swap_tmp = knnIdx(i,m);
                    knnIdx(i,m) = knnIdx(i,m+1);
                    knnIdx(i,m+1) =  idx_swap_tmp;
                    swapped = true;
                }
            }
            if (swapped == false){
                break;
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

void calculateNNKdist(Rcpp::NumericMatrix &knnDists, Rcpp::NumericVector &NNdistVector,  int &k, int &numRows, double &max_NNdist){
    double sum = 0;

    for( int i = 0; i < numRows; i++){
        for(int j = 0; j < k; j++){
            sum += knnDists(i,j);
        }
        NNdistVector(i) = sum;
        if (sum > max_NNdist){
            max_NNdist = sum;
        }
        sum = 0;
    }
}


void calculateRNN(Rcpp::IntegerMatrix &knnIdx, Rcpp::IntegerVector &nOfRNN, int &k, int &numRows){
    int count = 0;
    for(int i = 0; i < numRows; i++){
        for(int j = 0; j < numRows; j++){
            for (int m = 0; m < k; m++){
                if(knnIdx(j,m) == i){
                    count++;
                    break;
                }
            }
        }
        nOfRNN(i) = count;
        count = 0;
    }
}


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

//sort points
void sortPoints(Rcpp::NumericVector &INFLOVector, Rcpp::NumericVector &NNdistVector, Rcpp::IntegerVector &map, Rcpp::IntegerVector &newOrder,Rcpp::NumericVector &INFLOPlusNNdist, int &numRows ){


    for (int i = 0; i < numRows; i++){
        newOrder(i) = i;
        INFLOPlusNNdist(i) = INFLOVector(i) + NNdistVector(i);
    }

    bool swapped;
    double swap_tmp;
    int idx_swap_tmp;

    for (int i = 0; i < numRows - 1; i++){
        swapped = false;
        for (int j = 0; j < numRows - i - 1; j++){
            if(INFLOPlusNNdist(j) > INFLOPlusNNdist(j+1)){
                swap_tmp = INFLOPlusNNdist(j);
                INFLOPlusNNdist(j) = INFLOPlusNNdist(j+1);
                INFLOPlusNNdist(j+1) = swap_tmp;
                idx_swap_tmp = newOrder(j);
                newOrder(j) = newOrder(j+1);
                newOrder(j+1) = idx_swap_tmp;
                swapped = true;
            }
        }
        if (swapped == false){
            break;
        }
    }

    for(int i = 0; i < numRows; i ++){
        map(newOrder(i)) = i;
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






//END STRATIFICATION METHODS



//Expand Cluster
void expandCluster(Rcpp::IntegerMatrix &iskNeighbours, Rcpp::IntegerVector &allClusters,  Rcpp::IntegerVector &tmp_list, int pos_init, int cluster, int k, int &not_assigned, Rcpp::IntegerVector &nOfIsk, Rcpp::LogicalVector &border ){
    Rcpp::IntegerVector tmpIsk;
    if(nOfIsk(pos_init) > (k*2.0/3.0)){
        tmpIsk = iskNeighbours.row(pos_init);
        for (auto &iskIdx : tmpIsk){
            if(iskIdx != NA_INTEGER ){
                if (allClusters(iskIdx) == -1){
                    tmp_list.push_back(iskIdx);
                    allClusters(iskIdx) = cluster;
                    not_assigned--;
                    expandCluster(iskNeighbours,allClusters,tmp_list,iskIdx,cluster,k,not_assigned,nOfIsk, border);
                }
            }

        }
    }
    else{
        border(pos_init) = true;
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

//'
//' cIsdbscan
//'
//' cIsdbscan for both matrix and HDF5Matrix
//'
//'@param data numeric matrix or integer matrix or HDF5Matrix
//'@param k the number of clusters
//'@param batch_size the size of the mini batches
//'@return a list with the following attributes: centroids, WCSS_per_cluster, best_initialization, iters_per_initialization
//'@details
//'This function performs k-means clustering using mini batches.
//'
//'\strong{kmeans++}: kmeans++ initialization. Reference : http://theory.stanford.edu/~sergei/papers/kMeansPP-soda.pdf AND http://stackoverflow.com/questions/5466323/how-exactly-does-k-means-work
//'
//'\strong{random}: random selection of data rows as initial centroids
//'
//'@references
//'https://github.com/mlampros/ClusterR
//'
//' @export
// [[Rcpp::export]]
Rcpp::List cIsdbscan(SEXP data, int k, int batch_size, bool stratif){

    int data_n_rows = get_nrow(data);
    int data_n_cols = get_ncol(data);


    Rcpp::IntegerVector nOfIsk(data_n_rows);
    Rcpp::NumericMatrix submat(data_n_rows, data_n_cols);
    Rcpp::NumericMatrix knnDists(data_n_rows, k);
    Rcpp::IntegerMatrix knnIdx(data_n_rows, k);
    Rcpp::NumericVector tmp(data_n_cols);
    Rcpp::IntegerVector nOfRNN(data_n_rows);
    Rcpp::IntegerMatrix calIsk(data_n_rows,data_n_cols);
    Rcpp::LogicalVector border(data_n_rows);

    //wHEN DATA FITS MEMORY
    if (data_n_rows <= batch_size){
        auto final_matrix=beachmat::create_numeric_matrix(data);

        for ( unsigned int i = 0; i < data_n_rows; i++ ){
            final_matrix->get_row(i, tmp.begin());
            submat.row(i) = tmp;
        }
        calcKnnDists(submat, knnDists, knnIdx, k, data_n_rows, data_n_cols);
        sortKnnDists(knnDists, knnIdx, k, data_n_rows);
        calIsk = findIsk(knnIdx,k,data_n_rows,nOfIsk);
        if (!stratif){
            Rcpp::IntegerVector clusteringRes = calcCluster(calIsk,data_n_rows,k,nOfIsk,border);
            return Rcpp::List::create(Rcpp::Named("clusters") = clusteringRes);
        }
        else{
            double max_NNK = -1;
            double max_inflo = -1;
            double max_NNdist = -1;
            Rcpp::NumericVector NNdistVector(data_n_rows); //da creare solo se stratification == true;
            Rcpp::NumericVector INFLOVector(data_n_rows);
            Rcpp::IntegerVector map(data_n_rows);
            Rcpp::IntegerVector newOrder(data_n_rows);
            Rcpp::NumericVector INFLOPlusNNdist(data_n_rows);
            Rcpp::IntegerVector layerVector(data_n_rows);


            calculateNNKdist(knnDists,NNdistVector,k,data_n_rows, max_NNdist);
            calculateRNN(knnIdx,nOfRNN,k,data_n_rows);
            calculateINFLO(knnDists,knnIdx,nOfIsk,INFLOVector, k, data_n_rows, max_inflo);
            calculateINFLOCorrect(nOfRNN, nOfIsk, INFLOVector, data_n_rows);
            calculateINFLONorm(INFLOVector, data_n_rows, max_inflo, max_NNdist);
            sortPoints(INFLOVector,NNdistVector, map,newOrder, INFLOPlusNNdist, data_n_rows);
            stratifier(INFLOVector,NNdistVector,INFLOPlusNNdist,map,newOrder,layerVector,data_n_rows);

            Rcpp::IntegerVector clusteringRes = calcCluster(calIsk,data_n_rows,k,nOfIsk,border,map);
            return Rcpp::List::create(Rcpp::Named("clusters") = clusteringRes, Rcpp::Named("layer") = layerVector, Rcpp::Named("border") = border);
        }
    }
    else{

        auto final_matrix=beachmat::create_numeric_matrix(data);
        int rowsToRead = std::floor(batch_size/2);
        int howManyTimesTofitData = std::ceil(data_n_rows/rowsToRead);
        Rcpp::NumericVector tmpForj(data_n_cols);
        int offset;
        int counter = 0;

        Rcpp::IntegerVector pointMaxIdx(data_n_rows);
        Rcpp::NumericVector pointMaxValue(data_n_rows);

        for (int i = 0; i < data_n_rows; i++){
            pointMaxIdx(i) = -1;
            pointMaxValue(i) = 99999;
        }

        for (int i = 0; i < data_n_rows; i++){
            final_matrix->get_row(i, tmp.begin());
            offset = rowsToRead * counter;
            for(int j = rowsToRead * counter; j < (rowsToRead * (counter+1)) && j < data_n_rows; j++){
                final_matrix->get_row(j, tmpForj.begin());
                submat.row(j-offset) = tmpForj;
            }
            calcKnnDists(tmp ,submat, knnDists, knnIdx, pointMaxIdx, pointMaxValue, k, data_n_rows, data_n_cols, offset, rowsToRead, i);
            counter++;
        }

        sortKnnDists(knnDists, knnIdx, k, data_n_rows);
        calIsk = findIsk(knnIdx,k,data_n_rows,nOfIsk);
        if (!stratif){
            Rcpp::IntegerVector clusteringRes = calcCluster(calIsk,data_n_rows,k,nOfIsk,border);
            return Rcpp::List::create(Rcpp::Named("clusters") = clusteringRes);
        }
        else{
            double max_NNK = -1;
            double max_inflo = -1;
            double max_NNdist = -1;
            Rcpp::NumericVector NNdistVector(data_n_rows); //da creare solo se stratification == true;
            Rcpp::NumericVector INFLOVector(data_n_rows);
            Rcpp::IntegerVector map(data_n_rows);
            Rcpp::IntegerVector newOrder(data_n_rows);
            Rcpp::NumericVector INFLOPlusNNdist(data_n_rows);
            Rcpp::IntegerVector layerVector(data_n_rows);


            calculateNNKdist(knnDists,NNdistVector,k,data_n_rows, max_NNdist);
            calculateRNN(knnIdx,nOfRNN,k,data_n_rows);
            calculateINFLO(knnDists,knnIdx,nOfIsk,INFLOVector, k, data_n_rows, max_inflo);
            calculateINFLOCorrect(nOfRNN, nOfIsk, INFLOVector, data_n_rows);
            calculateINFLONorm(INFLOVector, data_n_rows, max_inflo, max_NNdist);
            sortPoints(INFLOVector,NNdistVector, map,newOrder, INFLOPlusNNdist, data_n_rows);
            stratifier(INFLOVector,NNdistVector,INFLOPlusNNdist,map,newOrder,layerVector,data_n_rows);

            Rcpp::IntegerVector clusteringRes = calcCluster(calIsk,data_n_rows,k,nOfIsk,border,map);
            return Rcpp::List::create(Rcpp::Named("clusters") = clusteringRes, Rcpp::Named("layer") = layerVector, Rcpp::Named("border") = border);
        }



    }


}


