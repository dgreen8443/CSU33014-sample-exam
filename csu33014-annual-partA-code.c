//
// CSU33014 Summer 2020 Additional Assignment
// Part A of a two-part assignment
//

// Please examine version each of the following routines with names
// starting partA. Where the routine can be vectorized, please
// complete the corresponding vectorized routine using SSE vector
// intrinsics. Where is cannot be vectorized...

// Note the restrict qualifier in C indicates that "only the pointer
// itself or a value directly derived from it (such as pointer + 1)
// will be used to access the object to which it points".


#include <immintrin.h>
#include <xmmintrin.h>
#include <stdio.h>

#include "csu33014-annual-partA-code.h"

void print_vector(__m128 vector){
    
    float * array = malloc(sizeof(float)*4);
    _mm_storeu_ps(array, vector);
    printf("%f %f %f %f \n",array[0], array[1],array[2],array[3]); 

}
float add_3(__m128 vector){
  float ret = 0;
  for (int i = 0; i < 3; i++){
    ret 
  }
}
void nan_remove(__m128 vector){
    float * array = malloc(sizeof(float)*4);
    _mm_storeu_ps(array, vector);
    for(int i = 0; i < 4; i++){
        if (array[i] != 0){
            array[i] = 1;
        }
    }
    vector = _mm_loadu_ps(&array[0]);
    
}

/****************  routine 0 *******************/

// Here is an example routine that should be vectorized
void partA_routine0(float * restrict a, float * restrict b,
		    float * restrict c) {
  for (int i = 0; i < 1024; i++ ) {
    a[i] = b[i] * c[i];
  }
}

// here is a vectorized solution for the example above
void partA_vectorized0(float * restrict a, float * restrict b,
		    float * restrict c) {
  __m128 a4, b4, c4;
  
  for (int i = 0; i < 1024; i = i+4 ) {
    b4 = _mm_loadu_ps(&b[i]);
    c4 = _mm_loadu_ps(&c[i]);
    a4 = _mm_mul_ps(b4, c4);
    _mm_storeu_ps(&a[i], a4);
  }
}

/***************** routine 1 *********************/

// in the following, size can have any positive value
float partA_routine1(float * restrict a, float * restrict b,
		     int size) {
  float sum = 0.0;
  
  for ( int i = 0; i < size; i++ ) {
    sum = sum + a[i] * b[i];
  }
  return sum;
}

// insert vectorized code for routine1 here
float partA_vectorized1(float * restrict a, float * restrict b,
		     int size) {
  // replace the following code with vectorized code
  float sum = 0.0;
  int over = size % 4;
  for(int i = 0; i < over; i++){
      sum = sum + (a[i] * b[i]);
  }
  __m128 a4, b4, c4, s4;
  float * sum_array = malloc(sizeof(float) * 4);
  s4 = _mm_set1_ps(0.0f);
  for (int i = over ; i < size; i = i + 4)
  {
    a4 = _mm_loadu_ps(&a[i]);
    b4 = _mm_loadu_ps(&b[i]);
    c4 = _mm_mul_ps(a4, b4);
    s4 = _mm_add_ps(c4, s4);

  }
  _mm_storeu_ps(sum_array, s4);
  for (int i = 0; i < 4; i++){
    sum = sum + sum_array[i];
  }
  free(sum_array);
  return sum;
}

/******************* routine 2 ***********************/

// in the following, size can have any positive value
void partA_routine2(float * restrict a, float * restrict b, int size) {
  for ( int i = 0; i < size; i++ ) {
    a[i] = 1 - (1.0/(b[i]+1.0));
  }
}

// in the following, size can have any positive value
void partA_vectorized2(float * restrict a, float * restrict b, int size) {
  // replace the following code with vectorized code
  __m128 a4, b4, c4, one4, tempone4, temp4;
  
  one4 = _mm_set1_ps(1);
  
  int over = size % 4;
  for(int i = 0; i < over; i++){
      a[i] = 1.0 - (1.0/(b[i]+ 1.0));
  }
  for (int i = over; i < size; i = i + 4)
  {
  
    b4 = _mm_loadu_ps(&b[i]);
    temp4 = _mm_add_ps(b4, one4);
    tempone4 = _mm_rcp_ps(temp4);

    a4 = _mm_sub_ps(one4, tempone4);
    _mm_storeu_ps(&a[i], a4);
    
  }
  
}

/******************** routine 3 ************************/

// in the following, size can have any positive value
void partA_routine3(float * restrict a, float * restrict b, int size) {
  for ( int i = 0; i < size; i++ ) {
    if ( a[i] < 0.0 ) {
      a[i] = b[i];
    }
  }
}

// in the following, size can have any positive value
void partA_vectorized3(float * restrict a, float * restrict b, int size) {
  // replace the following code with vectorized code
  
  __m128 a4, b4, mask4, zero4, temp4_a, temp4_b;
  zero4 = _mm_set1_ps(0.0f);
  int over = size % 4;
  for(int i = 0; i < over; i++){
      if (a[i] < 0.0){
          a[i] = b[i];
      }
  }
  for (int i = over; i < size; i = i + 4)
  {
      a4 = _mm_loadu_ps(&a[i]);
      b4 = _mm_loadu_ps(&b[i]);
      mask4 = _mm_cmplt_ps(a4,zero4);
      temp4_a = _mm_andnot_ps(mask4,a4);
      temp4_b = _mm_and_ps(b4, mask4);
      a4 = _mm_or_ps(temp4_a, temp4_b);
      _mm_storeu_ps(&a[i], a4);
  }
}

/********************* routine 4 ***********************/

// hint: one way to vectorize the following code might use
// vector shuffle operations
void partA_routine4(float * restrict a, float * restrict b,
		       float * restrict c) {
  for ( int i = 0; i < 2048; i = i+2  ) {
    a[i] = b[i]*c[i] - b[i+1]*c[i+1];
    a[i+1] = b[i]*c[i+1] + b[i+1]*c[i];
  }
}

void partA_vectorized4(float * restrict a, float * restrict b,
		       float * restrict  c) {
  // replace the following code with vectorized code
  __m128 a4, b4, c4, d4, e4, f4, g4, h4,b14,c14,d14,g14;
  /* for ( int i = 0; i < 2048; i = i+4  ) {
    b4 = _mm_loadu_ps(&b[i]);
    c4 = _mm_loadu_ps(&c[i]);
    d4 = _mm_mul_ps(b4,c4);
    e4 = d4;
    g4 = _mm_shuffle_ps(e4, d4, _MM_SHUFFLE(2, 3, 0, 1)); //[bi+1 * ci+1| bi *ci |bi+3 * ci+ 3| bi+2 * ci+2 ]
    h4 = _mm_sub_ps(d4, g4); //pos 3, 1 are a[i] and a[i+2] respectively
    print_vector(h4);
    b14 = _mm_loadu_ps(&b[i]);
    c14 = _mm_loadu_ps(&c[i]);
    d14 = _mm_mul_ps(b14, c14);

    g14 = _mm_sub_ps(d14, _mm_shuffle_ps(d14, d14, _MM_SHUFFLE(2, 3, 0, 1)));
    print_vector(g14);
    e4 = _mm_shuffle_ps(c4, c4, _MM_SHUFFLE(2, 3, 0, 1)); // [ci+1|ci|ci+3|ci+2]
    f4 = _mm_mul_ps(b4, e4); //[bi * ci+1|bi+1 * ci| bi+2 * ci+ 3| bi+3 *ci+2]
    g4 = _mm_shuffle_ps(f4, f4, _MM_SHUFFLE(2, 3, 0, 1)); // [bi+1* ci| bi *ci+1]

    d4 = _mm_add_ps(f4, g4); // pos 3, 1 are a[i+1] and a[i+3] respectively
    g4 = _mm_shuffle_ps(h4, d4, _MM_SHUFFLE(3,1,2,0)); //[bi*ci - bi+1ci+1|bi+2ci+2 - bi+3ci+3|bici+1 + bi+1ci|bi+2ci+3 + bi+3ci+2]
    a4 = _mm_shuffle_ps(g4, g4, _MM_SHUFFLE(3,1,2,0));
    _mm_storeu_ps(&a[i], a4);
    
  } */
  for (int i = 0; i < 2048; i = i + 4)
  {

    b14 = _mm_loadu_ps(&b[i]);
    c14 = _mm_loadu_ps(&c[i]);
    d14 = _mm_mul_ps(b14, c14);

    g14 = _mm_sub_ps(d14, _mm_shuffle_ps(d14, d14, _MM_SHUFFLE(2, 3, 0, 1)));
        
    f4 = _mm_mul_ps(b14, _mm_shuffle_ps(c14,c14,_MM_SHUFFLE(2,3,0,1)));                              //[bi * ci+1|bi+1 * ci| bi+2 * ci+ 3| bi+3 *ci+2]
    d4 = _mm_add_ps(f4, _mm_shuffle_ps(f4, f4, _MM_SHUFFLE(2, 3, 0, 1)) ); // pos 3, 1 are a[i+1] and a[i+3] respectively
    g4 = _mm_shuffle_ps(g14, d4, _MM_SHUFFLE(3, 1, 2, 0)); //[bi*ci - bi+1ci+1|bi+2ci+2 - bi+3ci+3|bici+1 + bi+1ci|bi+2ci+3 + bi+3ci+2]
    a4 = _mm_shuffle_ps(g4, g4, _MM_SHUFFLE(3, 1, 2, 0));
    _mm_storeu_ps(&a[i], a4);
  }
}

/********************* routine 5 ***********************/

// in the following, size can have any positive value
void partA_routine5(unsigned char * restrict a,
		    unsigned char * restrict b, int size) {
  for ( int i = 0; i < size; i++ ) {
    a[i] = b[i];
  }
}

void partA_vectorized5(unsigned char * restrict a,
		       unsigned char * restrict b, int size) {
  // replace the following code with vectorized code
  int over = size % 4;
  for (int i = 0; i < over ; i++){
    a[i] = b[i];
  }
  __m128 a4, b4;
  for ( int i = over; i < size; i =i+4 ) {
    b4 = _mm_loadu_ps((float *)&b[i]);
    a4 = b4;
    _mm_storeu_ps((float * )&a[i], a4);
  }
}

/********************* routine 6 ***********************/

void partA_routine6(float * restrict a, float * restrict b,
		       float * restrict c) {
  a[0] = 0.0;
  for ( int i = 1; i < 1023; i++ ) {
    float sum = 0.0;
    for ( int j = 0; j < 3; j++ ) {
      sum = sum +  b[i+j-1] * c[j];
    }
    a[i] = sum;
  }
  a[1023] = 0.0;
}

void partA_vectorized6(float * restrict a, float * restrict b,
		       float * restrict c) {
  // replace the following code with vectorized code
  __m128 a4, b4, c4, d4, e4, sum4;
  a[0] = 0.0;
 
   
  c4 = _mm_set1_ps(c[0]);
  d4 = _mm_set1_ps(c[1]);
  e4 = _mm_set1_ps(c[2]);
  for (int i = 4; i < 1023; i = i+4){
    // sum4 = _mm_set1_ps(0.0f);
    b4 = _mm_loadu_ps(&b[i-1]);
    a4 = _mm_add_ps(b4, c4);
    a4 = _mm_add_ps(a4, d4);
    a4 = _mm_add_ps(a4, e4);
    _mm_storeu_ps(&a[i], a4);
     
  }
  a[1023] = 0.0;
}
/*
////////////////////////////////////////////////////// works
void partA_vectorized6(float * restrict a, float * restrict b,
		       float * restrict c) {
  // replace the following code with vectorized code
  a[0] = 0.0;
  __m128 b4, c4, product4;
  // load c[0] into c4
  c4 = _mm_loadu_ps(&c[0]);
  float product[4] = {0,0,0,0};
  for (int i = 1; i < 1023; i++)
  {
    // load b[i + 0 - 1] to b[i + 2 - 1] into a b4
    b4 = _mm_loadu_ps(&b[i-1]);
    // multiply values in b4 all by c
    product4 = _mm_mul_ps(b4, c4);
    // put the sum of the first 3 products in a[i] n the same order 
    // as the sequential code so the answer is the same (we can simply ignore the 
    // fourth value)
    _mm_storeu_ps(&product[0], product4);
    a[i] = (product[0] + product[1]) + product[2];
  }
  a[1023] = 0.0;
}
*/

/*

// in the following, size can have any positive value
float partA_routine1(float * restrict a, float * restrict b,
		     int size) {
  float sum = 0.0;
  
  for ( int i = 0; i < size; i++ ) {
    sum = sum + a[i] * b[i];
  }
  //printf("sum = %f", sum);
  return sum;
}

// insert vectorized code for routine1 here
float partA_vectorized1(float * restrict a, float * restrict b,
		     int size) {
  // replace the following code with vectorized code
  //printf("size = %d\n", size);
  //printf("[%f, %f, %f, %f ... %f]\n", a[0], a[1], a[2], a[3],a[size-1]);
  //printf("[%f, %f, %f, %f ... %f]\n", b[0], b[1], b[2], b[3],b[size-1]);
  float sum[4];
  printf("size = %d\n", size);
  printf("mod 4 = %d\n", size % 4);

  __m128 sum4, a4, b4;
  sum4 = _mm_set1_ps(0.0);
  int i;
  int remainder = size%4;
  for (i = 0; i < size-remainder; i = i+4 ) {
    a4 = _mm_loadu_ps(&a[i]);
    b4 = _mm_loadu_ps(&b[i]);
    sum4 = _mm_add_ps(sum4, _mm_mul_ps(a4, b4));
  }
  _mm_storeu_ps(sum, sum4);
  int remres = 0;
  for (; i<size; i++){  
    printf("%d\n", i);
    remres = remres + a[i] * b[i];
  }
  //printf("[%f,%f,%f,%f]\n",sum[0],sum[1],sum[2],sum[3]);
  float result = sum[0] + sum[1] + sum[2] + sum[3] + remres;
  //printf("result = %f", result);
  return result ;
}


// in the following, size can have any positive value
void partA_routine2(float * restrict a, float * restrict b, int size) {
  for ( int i = 0; i < size; i++ ) {
    a[i] = 1 - (1.0/(b[i]+1.0));
  }
}

// in the following, size can have any positive value
void partA_vectorized2(float * restrict a, float * restrict b, int size) {
  // replace the following code with vectorized code
  __m128 a4, b4, ones;
  
  //printf("size = %d\n", size);
  //printf("mod 4 = %d\n", size % 4);

  ones = _mm_set1_ps(1.0);
  for ( int i = 0; i < size; i= i+4 ) {
    a4 = _mm_loadu_ps(&a[i]);
    b4 = _mm_loadu_ps(&b[i]);
    a4 = _mm_sub_ps(ones, _mm_div_ps(ones, _mm_add_ps(b4, ones)));
    _mm_storeu_ps(&a[i], a4);
  }

}


// in the following, size can have any positive value
void partA_routine3(float * restrict a, float * restrict b, int size) {
  //struct timeval start_time, stop_time;
  //long long compute_time;
  //TODO delete the timing data here before submission, only used to check if my function is faster
  //gettimeofday(&start_time, NULL);
  for ( int i = 0; i < size; i++ ) {
    if ( a[i] < 0.0 ) {
      a[i] = b[i];
    }
  }
  //gettimeofday(&stop_time, NULL);
  //compute_time = (stop_time.tv_sec - start_time.tv_sec) * 1000000L +
  //  (stop_time.tv_usec - start_time.tv_usec);
  //("Serial time to beat: %lld microseconds\n", compute_time);
}

// in the following, size can have any positive value
void partA_vectorized3(float * restrict a, float * restrict b, int size) {
  // replace the following code with vectorized code
  //struct timeval start_time, stop_time;
  //long long compute_time;
  //TODO delete the timing data here before submission, only used to check if my function is faster
  //gettimeofday(&start_time, NULL);
  __m128 a4, b4, zeros;
  //printf("size = %d\n", size);
  //printf("mod 4 = %d\n", size % 4);
  zeros = _mm_set1_ps(0.0);
  float cmpA[4];

  for ( int i = 0; i < size; i=i+4 ) {
    a4 = _mm_loadu_ps(&a[i]);
    b4 = _mm_loadu_ps(&b[i]);
    __m128 cmp = _mm_cmplt_ps(a4, zeros);
    _mm_storeu_ps(&cmpA[0], cmp);
    if (cmpA[0] != 0) a[i] = b[i];
    if (cmpA[1] != 0) a[i+1] = b[i+1];
    if (cmpA[2] != 0) a[i+2] = b[i+2];
    if (cmpA[3] != 0) a[i+3] = b[i+3];
  }
  //gettimeofday(&stop_time, NULL);
  //compute_time = (stop_time.tv_sec - start_time.tv_sec) * 1000000L +
  //  (stop_time.tv_usec - start_time.tv_usec);
  //printf("Vectorized Time: %lld microseconds\n", compute_time);
}


// hint: one way to vectorize the following code might use
// vector shuffle operations
void partA_routine4(float * restrict a, float * restrict b,
		       float * restrict c) {
  for ( int i = 0; i < 2048; i = i+2  ) {
    a[i] = b[i]*c[i] - b[i+1]*c[i+1];
    a[i+1] = b[i]*c[i+1] + b[i+1]*c[i];
    
  }
}

void partA_vectorized4(float * restrict a, float * restrict b,
		       float * restrict  c) {
  // replace the following code with vectorized code
  __m128 a4, b4, c4;


  for ( int i = 0; i < 2048; i = i+4  ) {
    b4 = _mm_loadu_ps(&b[i]);
    c4 = _mm_loadu_ps(&c[i]);

    //extract index from vector (code to extract and index from vector from @doug65536 https://stackoverflow.com/questions/5526658/intel-sse-why-does-mm-extract-ps-return-int-instead-of-float)
    //_MM_SHUFFLE pushes the idex to the least significant bit and _mm_)cvtss_f32 extrast the least significant 32 bits
    // extarct (as followed: //_mm_extract_ps(b4, 0);) doesnt work as it extracts to a general register which will interpret as an int not as a float
    float b0 = _mm_cvtss_f32(_mm_shuffle_ps(b4, b4, _MM_SHUFFLE(0, 0, 0, 0)));     //first element in vector b[i]
    float b1 = _mm_cvtss_f32(_mm_shuffle_ps(b4, b4, _MM_SHUFFLE(0, 0, 0, 1)));     //second element in vector b[i+1]
    float b2 = _mm_cvtss_f32(_mm_shuffle_ps(b4, b4, _MM_SHUFFLE(0, 0, 0, 2)));     //third element in vector b[i] but b[i] in "next iteration"
    float b3 = _mm_cvtss_f32(_mm_shuffle_ps(b4, b4, _MM_SHUFFLE(0, 0, 0, 3)));     //fourth element in vector b[i+1] but b[i+1] in "next iteration"
    float c0 = _mm_cvtss_f32(_mm_shuffle_ps(c4, c4, _MM_SHUFFLE(0, 0, 0, 0)));    
    float c1 = _mm_cvtss_f32(_mm_shuffle_ps(c4, c4, _MM_SHUFFLE(0, 0, 0, 1)));  
    float c2 = _mm_cvtss_f32(_mm_shuffle_ps(c4, c4, _MM_SHUFFLE(0, 0, 0, 2)));  
    float c3 = _mm_cvtss_f32(_mm_shuffle_ps(c4, c4, _MM_SHUFFLE(0, 0, 0, 3)));  

    //do the arithmetic on the vector
    float a0 = b0*c0 - b1*c1;       
    float a1 = b0*c1 + b1*c0;                        
    float a2 = b2*c2 - b3*c3; 
    float a3 = b2*c3 + b3*c2; 

    //set to a4 and add to array a
    a4 = _mm_set_ps(a3,a2,a1,a0);       //I have no idea why this has to be backwards, but it does (maybe just for this code but honestly i dont know and im too lazy to find out)
    _mm_storeu_ps(&a[i], a4);
  }
}

// in the following, size can have any positive value
void partA_routine5(unsigned char * restrict a,
		    unsigned char * restrict b, int size) {
  for ( int i = 0; i < size; i++ ) {
    a[i] = b[i];
  }
}

void partA_vectorized5(unsigned char * restrict a,
		       unsigned char * restrict b, int size) {
  // replace the following code with vectorized code
  
  //printf("size = %d\n", size);
  //printf("mod 4 = %d\n", size % 4);

  //as char cannot be used in vetors, use integer (and corrosponding ASCII table char)
  __m128i b4;  
  
  for ( int i = 0; i < size; i=i+4 ) {
    b4 = _mm_loadu_si128((const __m128i*)&b[i]);    //cast char to its appropriate ascii integer, i.e. A=65, B=66 etc...
    _mm_storeu_si128(&a[i], b4);
  }
}

*/