/*
If you found this code useful, please cite:
"On the computation of integrals over fixed-size rectangles of arbitrary dimension", Pattern Recognition Letters, May 2016.

Authors:
Omar Ocegueda (corresponding author: jomaroceguedag@gmail.com)
Oscar Dalmau
Eleftherios Garyfallidis
Maxime Descoteaux
Mariano Rivera
*/

#include<stdio.h>
#include<time.h>
#include <ctime>
#include<stdlib.h>
#define ABS(x) (((x)<0)?(-(x)):(x))
#define MAX(a, b) (((a)<(b))?(b):(a))
#define MIN(a, b) (((a)<(b))?(a):(b))

/// Create a new 4D volume of shape m[0], m[1], m[2], m[3]
double ****new_4d(int *m){
    double ****buffer = new double***[m[0]];
    for(int i=0; i<m[0]; ++i){
        buffer[i] = new double**[m[1]];
        for(int j=0; j<m[1]; ++j){
            buffer[i][j] = new double*[m[2]];
            for(int k=0;k<m[2]; ++k){
                buffer[i][j][k] = new double[m[3]];
            }
        }
    }
    return buffer;
}

/// Initializes the input 4D volume with random integers in [0, 1000)
void init_random(double ****I, int *m){
    srand(time(NULL));
    for(int i=0; i<m[0]; ++i){
        for(int j=0; j<m[1]; ++j){
            for(int k=0; k<m[2]; ++k){
                for(int l=0; l<m[3]; ++l){
                    I[i][j][k][l] = rand()%1000;
                }
            }
        }
    }
}

/// Computes the sum of V along the rectangle of shape m[...] and last corner (i, j, k, l)
double integrate_rectangle_4d(double ****V, int i, int j, int k, int l, int *m){
    double sum = 0;
    for(int ii=0; ii<m[0]; ++ii){
        if(ii > i){
            continue;
        }
        for(int jj=0; jj<m[1]; ++jj){
            if(jj > j){
                continue;
            }
            for(int kk=0; kk<m[2]; ++kk){
                if(kk>k){
                    continue;
                }
                for(int ll=0; ll<m[3]; ++ll){
                    if(ll<=l){
                        sum += V[i-ii][j-jj][k-kk][l-ll];
                    }
                }
            }
        }
    }
    return sum;
}

/// Computes the sum of the 4D volume V (of shape n[...]), over all rectangles of shape m[...]
/// Using the **direct** algorithm (for each rectangle, directly sum all its voxels)
/// The callback function will be called for each rectangle once its sum has been determined
/// This is to illustrate the 'on-the-fly' implementation where no extra temporary memory is
/// required.  
void compute_sums_4d_direct(
    double ****V,
    int *vsize,
    int *m,
    void (*callback)(int, int, int, int, double)){
    
    for(int i=0; i<vsize[0]; ++i){
        for(int j=0 ;j<vsize[1]; ++j){
            for(int k=0; k<vsize[2]; ++k){
                for(int l=0; l<vsize[3]; ++l){
                    // Integrate over rectangle with last corner (i, j, k, l)
                    double sum = integrate_rectangle_4d(V, i, j, k, l, m);
                    callback(i, j, k, l, sum);
                }
            }
        }
    }
}


/// Computes the sum of the 4D volume V (of shape n[...]), over all rectangles of shape m[...]
/// Using the **fast** algorithm (take advantage of the sum of rectangles with last corner in the
/// previous hyper-slice)
/// The callback function will be called for each rectangle once its sum has been determined
/// This is to illustrate the 'on-the-fly' implementation where no extra temporary memory is
/// required (except the two hyper slices to keep track of previous and current sums).
void compute_sums_4d_fast(
    double ****I,
    int *n,
    int *m,
    void (*callback)(int, int, int, int, double)){
    
    // Temporary buffer (only two hyper-slices of the 4D buffer)
    int temp_size[4] = {2, n[1], n[2], n[3]}; 
    double ****temp = new_4d(temp_size);

    // Keep track of the current slice position within the temporary buffer: either 0 or 1
    int current_slice = 1;
    for(int i=0; i<n[0]; ++i){
        current_slice = 1 - current_slice;
        // i_prev is the position of the previous slice within the temporary
        // buffer. This would be simply i-1, but we are using only 
        // 2 hyper-slices to save memory. The previous slice is the one that
        // is *not* the current one: 1 - current_slice
        int prev_i = 1 - current_slice;
        
        
        for(int j=0; j < n[1]; ++j){
            int prev_j = j - 1;
            for(int k=0; k < n[2]; ++k){
                int prev_k = k - 1;
                for(int l=0; l < n[3]; ++l){
                    int prev_l= l - 1;
                    // Initialize the sum with the last corner
                    double sum = I[i][j][k][l];
                    
                    //////////////////////////////////////////////////////////
                    // Add signed sub-volumes displaced by every q in {0,1}^4
                    //////////////////////////////////////////////////////////
                    if(i > 0){
                        sum += temp[prev_i][j][k][l];  // q=(1,0,0,0)
                        if(j > 0){
                            sum -= temp[prev_i][prev_j][k][l];  // q=(1,1,0,0)
                            if(k > 0){
                                sum += temp[prev_i][prev_j][prev_k][l];  // q=(1,1,1,0)
                                if(l > 0){
                                    sum -= temp[prev_i][prev_j][prev_k][prev_l];  // q=(1,1,1,1)
                                }
                            }
                            if(l > 0){
                                sum += temp[prev_i][prev_j][k][prev_l];  // q=(1,1,0,1)
                            }
                        }

                        if(k > 0){
                            sum -= temp[prev_i][j][prev_k][l];  // q=(1,0,1,0)
                            if(l>0){
                                sum += temp[prev_i][j][prev_k][prev_l];  // q=(1,0,1,1)
                            }
                        }
                        
                        if(l > 0){
                            sum -= temp[prev_i][j][k][prev_l];  // q=(1,0,0,1)
                        }
                    }
                    if(j > 0){
                        sum += temp[current_slice][prev_j][k][l];  // q=(0,1,0,0)
                        if(k > 0){
                            sum -= temp[current_slice][prev_j][prev_k][l];  // q=(0,1,1,0)
                            if(l > 0){
                                sum += temp[current_slice][prev_j][prev_k][prev_l];  // q=(0,1,1,1)
                            }
                        }
                        if(l > 0){
                            sum -= temp[current_slice][prev_j][k][prev_l];  // q=(0,1,0,1)
                        }
                    }
                    if(k > 0){
                        sum += temp[current_slice][j][prev_k][l];  // q=(0,0,1,0)
                        if(l>0){
                            sum -= temp[current_slice][j][prev_k][prev_l];  // q=(0,0,1,1)
                        }
                    }
                    if(l > 0){
                        sum += temp[current_slice][j][k][prev_l];  // q=(0,0,0,1)
                    }
                    
                    ///////////////////////////////////////////////////////////
                    // Add signed corners displaced by every q*m, q in {0,1}^4
                    ///////////////////////////////////////////////////////////
                    if(i >= m[0]){
                        sum -= I[i-m[0]][j][k][l];  // q=(1,0,0,0)
                        if(j >= m[1]){
                            sum += I[i-m[0]][j-m[1]][k][l];  // q=(1,1,0,0)
                            if(k >= m[2]){
                                sum -= I[i-m[0]][j-m[1]][k-m[2]][l];  // q=(1,1,1,0)
                                if(l >= m[3]){
                                    sum += I[i-m[0]][j-m[1]][k-m[2]][l-m[3]];  // q=(1,1,1,1)
                                }
                            }
                            if(l >= m[3]){
                                sum -= I[i-m[0]][j-m[1]][k][l-m[3]];  // q=(1,1,0,1)
                            }
                        }
                        
                        if(k >= m[2]){
                            sum += I[i-m[0]][j][k-m[2]][l];  // q=(1,0,1,0)
                            if(l >= m[3]){
                                sum -= I[i-m[0]][j][k-m[2]][l-m[3]];  // q=(1,0,1,1)
                            }
                        }
                        
                        if(l >= m[3]){
                            sum += I[i-m[0]][j][k][l-m[3]];  // q=(1,0,0,1)
                        }
                    }
                    if(j >= m[1]){
                        sum -= I[i][j-m[1]][k][l];  // q=(0,1,0,0)
                        if(k >= m[2]){
                            sum += I[i][j-m[1]][k-m[2]][l];  // q=(0,1,1,0)
                            if(l >= m[3]){
                                sum -= I[i][j-m[1]][k-m[2]][l-m[3]];  // q=(0,1,1,1)
                            }
                        }
                        if(l >= m[3]){
                            sum += I[i][j-m[1]][k][l-m[3]];  // q=(0,1,0,1)
                        }
                    }
                    if(k >= m[2]){
                        sum -= I[i][j][k-m[2]][l];  // q=(0,0,1,0)
                        if(l >= m[3]){
                            sum += I[i][j][k-m[2]][l-m[3]];  // q=(0,0,1,1)
                        }
                    }
                    if(l >= m[3]){
                        sum -= I[i][j][k][l-m[3]];  // q=(0,0,0,1)
                    }

                    // Process current sum and save it in the temporary buffer
                    callback(i, j, k, l, sum);
                    temp[current_slice][j][k][l] = sum;
                }
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// Compare the two implementations: direct vs. fast
//////////////////////////////////////////////////////////////////////////////////////

// The process for each voxel is simply "save" the sum in a temporary buffer (note that
// this extra buffer may not be necessary in some applications where we only need to
// process the rectangles 'on-the-fly', which may save a significant amount of memory).
double ****A_direct;
double ****A_fast;

void process_rectangle_direct(int i, int j, int k, int l, double sum){
    A_direct[i][j][k][l] = sum;
}

void process_rectangle_fast(int i, int j, int k, int l, double sum){
    A_fast[i][j][k][l] = sum;
}


void test_fast_integrals(void){
    int n[] = {50, 50, 50, 50};  // Volume size
    int m[] = {3, 4, 5, 6};  // Rectangle size
    double ****I = new_4d(n);
    init_random(I, n);
    
    // Run both implementations
    A_direct = new_4d(n);
    A_fast = new_4d(n);
    
    
    printf("Computing using fast implementation...");
    fflush(stdout);
    clock_t begin = clock();
    compute_sums_4d_fast(I, n, m, process_rectangle_fast);
    clock_t end = clock();
    double elapsed = double(end - begin) / CLOCKS_PER_SEC;
    printf ("Elapsed: %.2lf s.\n", elapsed);
    fflush(stdout);

    printf("Computing using direct implementation...");
    fflush(stdout);
    begin = clock();
    compute_sums_4d_direct(I, n, m, process_rectangle_direct);
    end = clock();
    elapsed = double(end - begin) / CLOCKS_PER_SEC;
    printf ("Elapsed: %.2lf s.\n", elapsed);
    fflush(stdout);
    
    // Verify accuracy
    double max_difference = -1;
    double min_I = I[0][0][0][0];
    double max_I = I[0][0][0][0];
    for(int i=0; i<n[0]; ++i){
        for(int j=0; j<n[1]; ++j){
            for(int k=0; k<n[2]; ++k){
                for(int l=0; l<n[3]; ++l){
                    double d = ABS(A_direct[i][j][k][l] - A_fast[i][j][k][l]);
                    max_difference = MAX(max_difference, d);
                    min_I = MIN(min_I, I[i][j][k][l]);
                    max_I = MAX(max_I, I[i][j][k][l]);
                }
            }
        }
    }
    printf("Dynamic range of input:[%lf, %lf]\n", min_I, max_I);
    printf("Maximum difference:%lf\n", max_difference);
}

int main(){
    test_fast_integrals();
    return 0;
}


