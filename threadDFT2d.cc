// Threaded two-dimensional Discrete FFT transform
// AISWARIA SUKUMARAN NAIR
// ECE8893 Project 2

#include <iostream>
#include <string>
#include <math.h>
#include<pthread.h>
#include <sys/time.h>
#include<stdint.h>
#include "Complex.h"
#include "InputImage.h"
using namespace std;

// You will likely need global variables indicating how
// many threads there are, and a Complex* that points to the
// 2d image being transformed.


// |Global Variables|
//----------------------
//int sum;   //shared by the threads
//void *runner(void* param);   //thread

Complex* ImageData;
Complex* WeightMatrixValue;
//Complex* t_column, t_row;
//Complex* intermediate;

//Image details
unsigned ImageWidth;
unsigned ImageHeight;
//unsigned ThreadCount;
unsigned ImageSize;
//unsigned ImageDetails;
unsigned N;
// Global variable that threads use to count elements
unsigned elementCount = 0;
bool flag=0;
// Total number of threads
unsigned nThreads = 16;
unsigned runningThreads=0;		//currently excuting threads


// The mutex and condition variables allow the main thread to
// know when all helper threads are completed.

//pthread_mutex_t startCountMutex;
pthread_mutex_t exit_mutex;
pthread_mutex_t img_mutex;
pthread_mutex_t executing_mutex;
pthread_mutex_t cout_mutex;
//pthread_mutex_t elementCountMutex;
pthread_cond_t exit_cond;

//Barrier declaration
pthread_barrier_t barrier;

//int startCount;

//void WeightValuesMatrix(Complex * Array, unsigned size);
void MyBarrier_Init();
void MyBarrier();

void reorderArray(Complex * Array2D, unsigned size, unsigned width);
void Transform1D(Complex* h, unsigned N, unsigned width);
void TransformWrapper(Complex* h, unsigned width, unsigned height, unsigned startIndex, unsigned numberOfRows, bool flag);
void WeightValuesMatrix(Complex * Array2D, unsigned size, bool flag);
void InverseTransform1D(Complex* h, unsigned N, unsigned width);
void Transform2D(const char* inputFN, bool flag) ;
void* Transform2DThreadFwd(void* v);
void* Transform2DThreadReverse(void* v);
//void InverseTransform2D(const char* inputFN) ;




int main(int argc, char** argv)
{

  string fn("Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  
 // string fn2("after2d-correct.txt"); // default file name
//  if (argc > 1) fn2 = string(argv[1]);  // if name specified on cmd line
  // MPI initialization here
  
  	MyBarrier_Init();

	//Lock the exit mutex
	pthread_mutex_lock(&exit_mutex);
	
  	Transform2D(fn.c_str(), 1); // Perform the transform
 	 //Transform2D(fn2.c_str(),0); // Perform the transform
  
  return 0;

}  
	

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform

unsigned ReverseBits(unsigned v)
{ //  Provided to student
  unsigned n = N; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value
   
  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

// GRAD Students implement the following 2 functions.
// Undergrads can use the built-in barriers in pthreads.

// Call MyBarrier_Init once in main

void MyBarrier_Init()// you will likely need some parameters)
{	
	//Pthreads
	if(pthread_barrier_init(&barrier, 0, nThreads))
	{
		cout << "Could not initialize pthreads barrier.\n";
		fflush(stdout);
		//exit(-1);
	}
//Initialize mutexes & conditions
	pthread_mutex_init(&img_mutex, 0);
	pthread_mutex_init(&exit_mutex, 0);
	pthread_mutex_init(&cout_mutex, 0);
	pthread_mutex_init(&executing_mutex, 0);
	pthread_cond_init(&exit_cond, 0);	
	
}


// Each thread calls MyBarrier after completing the row-wise 



void MyBarrier() // Again likely need parameters
{
	
	//Pthreads
	int32_t rc = pthread_barrier_wait(&barrier);
	
    if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
    {
		pthread_mutex_lock(&cout_mutex);
        cout << "Could not wait on barrier\n";
		fflush(stdout);
		pthread_mutex_unlock(&cout_mutex);
        
}
 }
 
void Transform1D(Complex* h, unsigned N, unsigned width)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)
  // "N/NumberOfPoints"--How many transforms
  
  //Re-order array elemts
  //Calculate Transform
  unsigned temp, wMatrixIndex;
  	
	for(unsigned NumberOfPoints = 1; NumberOfPoints < N; NumberOfPoints *= 2)
	{
		for(unsigned index = 0; index < NumberOfPoints; ++index)
		{
			wMatrixIndex = index * N /(NumberOfPoints*2);
			
			for(unsigned j = index; j < N; j+=NumberOfPoints*2 )
			{
				temp = j + NumberOfPoints;
				
				
				if(width!=0)
				{
				pthread_mutex_lock(&img_mutex);
				Complex intermediate = WeightMatrixValue[wMatrixIndex] * h[temp*width];
				h[temp*width] = h[j*width] - intermediate;
				h[j*width] = h[j*width] + intermediate;	
				pthread_mutex_unlock(&img_mutex);
			}
			else
			{
			 	pthread_mutex_lock(&img_mutex);
				Complex intermediate = WeightMatrixValue[wMatrixIndex] * h[temp];
				h[temp] = h[j] - intermediate;
				h[j] = h[j] + intermediate;	
				}
				pthread_mutex_unlock(&img_mutex);
			}
		}
	}
	
cout<<"Tranform1D Called"<<endl;	
	
} 


void InverseTransform1D(Complex* h, unsigned N, unsigned width)
{
  // Implement the efficient Danielson-Lanczos DFT here.
  // "h" is an input/output parameter
  // "N" is the size of the array (assume even power of 2)
  // "N/NumberOfPoints"--How many transforms
  
  //Re-order array elemts
  //Calculate Transform
  unsigned temp, wMatrixIndex;
 //	N=ImageWidth;
 
	for(unsigned NumberOfPoints = 1; NumberOfPoints < N; NumberOfPoints *= 2)
	{
		for(unsigned index = 0; index < NumberOfPoints; ++index)
		{
			wMatrixIndex = index * N /(NumberOfPoints*2);
			
			for(unsigned j = index; j < N; j+=NumberOfPoints*2 )
			{
				temp = j + NumberOfPoints;
				
				
				if(width!=0)
				{
				pthread_mutex_lock(&img_mutex);
				Complex intermediate = WeightMatrixValue[wMatrixIndex] * h[temp*width];
				h[temp*width] = h[j*width] - intermediate;
				h[j*width] = h[j*width] + intermediate;	
				pthread_mutex_unlock(&img_mutex);
			}
			else
			{
			 	pthread_mutex_lock(&img_mutex);
				Complex intermediate = WeightMatrixValue[wMatrixIndex] * h[temp];
				h[temp] = h[j] - intermediate;
				h[j] = h[j] + intermediate;	
				}
				pthread_mutex_unlock(&img_mutex);
			}
		}
	}
	
	
	for(unsigned iter = 0; iter < N; iter++)
  {
    h[iter].real = h[iter].real / N;
	h[iter].imag = h[iter].imag / N;
  }  
  cout<<"Tranform1D Inverse Called"<<endl;
	
} 


void TransformWrapper(Complex* h, unsigned width, unsigned height, unsigned startIndex, unsigned numberOfRows,bool flag)
{
	if(height!=0)
	{
		if(flag)
	{
	for(unsigned i = startIndex; i < (startIndex + numberOfRows) ; i++)
	{
		reorderArray(&h[i], height, width);
		Transform1D(&h[i], height, width);
	}
	}
	else
	{
	for(unsigned i = startIndex; i < (startIndex + numberOfRows) ; i++)
	{
		reorderArray(&h[i], height, width);
		InverseTransform1D(&h[i], height, width);
	}
	}
	
	}
	else
	{
		if(flag)
	{
		
		for(unsigned i = startIndex; i < startIndex + numberOfRows; i++)
	{
		reorderArray(&h[i*width],width,0);
		Transform1D(&h[i*width], width,0);
	}
	}

else
	{
		
		for(unsigned i = startIndex; i < startIndex + numberOfRows; i++)
	{
		reorderArray(&h[i*width],width,0);
		InverseTransform1D(&h[i*width], width,0);
	}
	}

	}
	cout<<"Tranform Wrapper Called"<<endl;
	
}
 // This is the thread starting point.  "v" is the thread number
  // Calculate 1d DFT for assigned rows
  // wait for all to complete
  // Calculate 1d DFT for assigned columns
  // Decrement active count and signal main if all complete
  
  void Transform2D(const char* inputFN, bool flag) 
{
	InputImage image(inputFN);  // Create the helper object for reading the image

	//Set global parameters / get data
	ImageHeight = image.GetHeight();
	ImageWidth = image.GetWidth();
	ImageSize = ImageHeight * ImageWidth;
	ImageData = image.GetImageData();
	N = ImageWidth;
	//Precalculate weighting values to save time later
	WeightMatrixValue = new Complex[ImageHeight/2];
	
//	if(flag)
	WeightValuesMatrix(WeightMatrixValue, ImageWidth, flag);
	//else
	//WeightValuesMatrix(WeightMatrixValue, ImageWidth, 0);
	
	// Create 16 threads
	for(uint8_t i = 0; i < nThreads; i++)
	{
		pthread_mutex_lock(&executing_mutex);
		runningThreads++;
		pthread_mutex_unlock(&executing_mutex);
		
		pthread_t individualThread;
		
		if(flag)
		pthread_create(&individualThread, 0, Transform2DThreadFwd, (void*)i);
		else
		pthread_create(&individualThread, 0, Transform2DThreadReverse, (void*)i);
	}

	// Wait for all threads complete
	pthread_cond_wait(&exit_cond, &exit_mutex);
	
	// Write the transformed data
	if(flag)
	image.SaveImageData("Tower-DFT2D.txt", ImageData, 1024, 1024);  
	else
	image.SaveImageData("Tower-IDFT2D.txt", ImageData, 1024, 1024);
	
	//Free up the used memory to avoid memory leaks
	cout<<"Transform 2D Called"<<endl;
	delete ImageData; 
	delete WeightMatrixValue;
}

/*void InverseTransform2D(const char* inputFN) 
{
	InputImage image(inputFN);  // Create the helper object for reading the image

	//Set global parameters / get data
	ImageHeight = image.GetHeight();
	ImageWidth = image.GetWidth();
	ImageSize = ImageHeight * ImageWidth;
	ImageData = image.GetImageData();
	N = ImageWidth;
	//Precalculate weighting values to save time later
	WeightMatrixValue = new Complex[ImageHeight / 2];
	
	WeightValuesMatrix(WeightMatrixValue, ImageWidth, 0);
	
	// Create 16 threads
	for(unsigned i = 0; i < nThreads; i++)
	{
		pthread_mutex_lock(&executing_mutex);
		runningThreads++;
		pthread_mutex_unlock(&executing_mutex);
		
		pthread_t individualThread;
		
		pthread_create(&individualThread, 0, Transform2DThread, (void*)i);
	}

	// Wait for all threads complete
	pthread_cond_wait(&exit_cond, &exit_mutex);
	
	// Write the transformed data
	
	image.SaveImageData("Tower-IDFT2D.txt", ImageData, 1024, 1024);
	
	//Free up the used memory to avoid memory leaks
	
	delete ImageData; 
	delete WeightMatrixValue;
}
  
*/  

void* Transform2DThreadFwd(void* v)
{ 
  
    uint64_t thread = (uint64_t) v;
//	uint16_t rank = (uint16_t) t;
	uint16_t rank = (uint16_t) thread;
	unsigned numberOfRows = ImageHeight / nThreads;
	unsigned startIndex = numberOfRows * rank;
	
	
    TransformWrapper(ImageData, ImageWidth, 0, startIndex,numberOfRows,1);
	//Enter the barrier and wait
	MyBarrier();
    TransformWrapper(ImageData,  ImageWidth, ImageHeight, startIndex, numberOfRows,1);
    
    
//Decrement thread count
	pthread_mutex_lock(&executing_mutex);
	
//Indicate to main fucntion if last thread
	if(--runningThreads == 0)
	{
		
		pthread_mutex_lock(&exit_mutex);
		pthread_cond_signal(&exit_cond);
		pthread_mutex_unlock(&exit_mutex);
	}
	
	pthread_mutex_unlock(&executing_mutex);
cout<<"TranformThreadFwd  Called"<<endl;
	return 0;

}

//----------------------------Inverse------------------------

void* Transform2DThreadReverse(void* v)
{ 
  
    uint64_t thread = (uint64_t) v;

	uint16_t rank = (uint16_t) thread;
	unsigned numberOfRows = ImageHeight / nThreads;
	unsigned startIndex = numberOfRows * rank;
	
	
    
	TransformWrapper(ImageData, ImageWidth, 0, startIndex,  numberOfRows,0);
	//Enter the barrier and wait
	MyBarrier();
    TransformWrapper(ImageData,  ImageWidth, ImageHeight, startIndex, numberOfRows,0);
    
//Dcerememnt count of threads
	pthread_mutex_lock(&executing_mutex);
	
	//Check to see if we're the last thread to exit
	if(--runningThreads == 0)
	{
		//Indicator to main function
		pthread_mutex_lock(&exit_mutex);
		pthread_cond_signal(&exit_cond);
		pthread_mutex_unlock(&exit_mutex);
	}
	
	pthread_mutex_unlock(&executing_mutex);
	
cout<<"Transform Thread Reverse Called"<<endl;
	return 0;

}

//-----------------------------------------------------------------------------------------


void WeightValuesMatrix(Complex * Array2D, unsigned size, bool flag)
{
	//Precalculate weight values
	
	if(flag)
	{
		
	for(unsigned i = 0; i < size; i++)
	{
		double realvalue, imagvalue;
		
		realvalue = cos(2 * M_PI * i / size);
		imagvalue = -1 * sin(2 * M_PI * i / size);
		
		Array2D[i].real = realvalue;
		Array2D[i].imag = imagvalue;
	}
	}
	else
	{

		for(unsigned i = 0; i < size; i++)
	{
			double realvalue, imagvalue;
		
		realvalue = cos(2 * M_PI * i / size);
		imagvalue = sin(2 * M_PI * i / size);
		
		Array2D[i].real = realvalue/N;
		Array2D[i].imag = imagvalue/N;
	}
	
	}
	
cout<<"Weight Values Matrix Called"<<endl;	
	
}

void reorderArray(Complex * Array2D, unsigned size, unsigned width)
{
	//Re-ordering the array to use the algorithm
	unsigned j;
	for(unsigned i = 0; i < size; i++)
	{
//	Get the bits, reverse them and then form the thread!
		j = ReverseBits(i); 
		if(j > i)
		{
		
			if (width!=0)
		{
			pthread_mutex_lock(&img_mutex);
   			Complex thread = Array2D[i * width];
			Array2D[i * width] = Array2D[j * width];
			Array2D[j * width] = thread;
			pthread_mutex_unlock(&img_mutex);
		}
		else
		{
			pthread_mutex_lock(&img_mutex);
			Complex thread = Array2D[i];
			Array2D[i] = Array2D[j];
			Array2D[j] = thread;
			pthread_mutex_unlock(&img_mutex);
			}
			
		
		}
	}
	cout<<"Reorder array Called"<<endl;
}
