#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "time.h"


/**** Timing Macros *****/
struct timespec tick_clockData;
struct timespec tock_clockData;

#define TICK()	{ clock_gettime(CLOCK_MONOTONIC, &tick_clockData);}
#define TOCK(str) { clock_gettime(CLOCK_MONOTONIC, &tock_clockData);\
                	printf("%s\t%f milliseconds\n",str,\
                	(((double) tock_clockData.tv_sec + tock_clockData.tv_nsec / 1000000000.0) - \
                	((double) tick_clockData.tv_sec + tick_clockData.tv_nsec / 1000000000.0))\
                	* 1000);}
/**** OpenCL necessary Defines ****/
#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#include <CL/cl2.hpp>

/**** OpenCL API variables ****/
cl_int err;
cl_uint check_status = 0;

cl_platform_id platform_id;
cl_platform_id platforms[16];
cl_uint platform_count;
cl_uint platform_found = 0;
char cl_platform_vendor[1001];

cl_uint num_devices;
cl_uint device_found = 0;
cl_device_id devices[16];
char cl_device_name[1001];
cl_device_id device_id;

cl_context context;
cl_command_queue q;
cl_program program;

cl_kernel gibbsSamplerKernel; // Handler for hardware kernel
// Alternative for 2 compute units
// cl_kernel gibbsSamplerKernel_1;
// cl_kernel gibbsSamplerKernel_2;

cl_mem pt_in[1];
cl_mem pt_out[1];
cl_int status;

/*** Variables used as kernel arguments ***/
cl_mem input_buffer;
cl_mem output_buffer;

// Input Array
char *motifs_1d;
// Output Array
char *best_motifs_output;


/***********************************************************
 * Function:  load_file_to_memory
 * ------------------------------
 * Loads the hardware binary file (*.xclbin) into memory.
 * *********************************************************/
cl_uint load_file_to_memory(const char *filename, char **result)
{
	cl_uint size = 0;
	FILE *f = fopen(filename, "rb");
	if (f == NULL) {
    	*result = NULL;
    	return -1; // -1 means file opening fail
	}
	fseek(f, 0, SEEK_END);
	size = ftell(f);
	fseek(f, 0, SEEK_SET);
	*result = (char *)malloc(size+1);
	if (size != fread(*result, sizeof(char), size, f)) {
    	free(*result);
    	return -2; // -2 means file reading fail
	}
	fclose(f);
	(*result)[size] = 0;
	return size;
}


int lenHelper(unsigned x)
{
	if (x >= 1000000000) return 10;
	if (x >= 100000000)  return 9;
	if (x >= 10000000)   return 8;
	if (x >= 1000000)	return 7;
	if (x >= 100000) 	return 6;
	if (x >= 10000)  	return 5;
	if (x >= 1000)   	return 4;
	if (x >= 100)    	return 3;
	if (x >= 10)     	return 2;
	return 1;
}


/***********************************************************
 * Function:  main
 * ---------------------------------------------------------
 * Main function.
 * *********************************************************/
int main(int argc, char *argv[]){
  // Input and output array size
  size_t n0,len;
  cl_uint iplat, n_i0;
  char *hw_binary_path,*kernelbinary;
  char buffer[2048];
  int argcounter;

  FILE *input;
  int k, t, n; // k is the size of k-mers ; n is the times we need to iterate ; t is the number of motifs;
  char *first_motif_buffer, *check;
  //char *motifs_1d;
  int bytes_read;
  int motif_length, number_of_motifs;
  int i;

  if ( argc < 2)
  {
	printf("1 Argument needed : <*.xclbin path>");
	return EXIT_FAILURE;
  }
  hw_binary_path = argv[1];

	/**********************************************
   *
   *    		 Xilinx OpenCL Initialization
   *
 	* We must follow specific steps to get the necessary
 	* information and handlers, in order to be able
 	* to use the available accelerator device (FPGA).
 	* After every step, we always check for any errors
 	* that might have occured. In case of error, the
 	* program aborts and exits immediately.
   * *********************************************/



	/**************************************************
  * Step 1:
	* Get available OpenCL platforms and devices.
	* In our case, is a Xilinx FPGA device.
	* If the underlying platform has other accelerators
	* available, we could use them too (e.g. GPU, CPU).
  **************************************************/
    err = clGetPlatformIDs(16, platforms, &platform_count);
  if (err != CL_SUCCESS)
  {
  	printf("Error: Failed to find an OpenCL platform!\n");
  	return EXIT_FAILURE;
  }

  printf("INFO: Found %d platforms\n", platform_count);

    for (iplat=0; iplat<platform_count; iplat++)
  {
   	 err = clGetPlatformInfo(platforms[iplat], CL_PLATFORM_VENDOR, 1000, (void *)cl_platform_vendor,NULL);
	if (err != CL_SUCCESS)
	{
  	printf("Error: clGetPlatformInfo(CL_PLATFORM_VENDOR) failed!\n");
  	return EXIT_FAILURE;
	}
	/*** We are interested for Xilinx devices ***/
	if (strcmp(cl_platform_vendor, "Xilinx") == 0)
	{
  	printf("INFO: Selected platform %d from %s\n", iplat, cl_platform_vendor);
  	platform_id = platforms[iplat];
  	platform_found = 1; // There is only one available.
	}
  }

  if (!platform_found)
  {
	printf("ERROR: Platform Xilinx not found. Exit.\n");
	return EXIT_FAILURE;
  }

    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ACCELERATOR, 16, devices, &num_devices);
  printf("INFO: Found %d devices\n", num_devices);
  if (err != CL_SUCCESS)
  {
	printf("ERROR: Failed to create a device group!\n");
	return EXIT_FAILURE;
  }

    device_id = devices[0]; // we have only one device

    // ---------------------------------------------------------------
    // Step 2 : Create Context
    // ---------------------------------------------------------------
    context = clCreateContext(0,1,&device_id,NULL,NULL,&err);
    if (!context)
  {
	printf("Error: Failed to create a compute context!\n");
	return EXIT_FAILURE;
  }

    // ---------------------------------------------------------------
    // Step 3 : Create Command Queue
    // ---------------------------------------------------------------
    q = clCreateCommandQueue(context, device_id, CL_QUEUE_PROFILING_ENABLE, &err);
	// ---------------------------------------------------------------
    // Step 3 for 2 compute units: Create Command Queue
    // ---------------------------------------------------------------
	// q = clCreateCommandQueue(context, device_id, CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err);
    if (!q)
  {
	printf("Error: Failed to create a command q! Error code: %i\n",err);
	return EXIT_FAILURE;
  }

    // ---------------------------------------------------------------
    // Step 4 : Load Hardware Binary File (*.xclbin) from disk
    // ---------------------------------------------------------------
  n_i0 = load_file_to_memory(hw_binary_path, (char **) &kernelbinary);
  if (n_i0 < 0)
  {
	printf("failed to load kernel from xclbin: %s\n", hw_binary_path);
	exit(EXIT_FAILURE);
  }
  n0 = n_i0;

	// ---------------------------------------------------------------
    // Step 5 : Create program using the loaded hardware binary file
    // ---------------------------------------------------------------
  program = clCreateProgramWithBinary(context, 1, &device_id, &n0, (const unsigned char **) &kernelbinary, &status, &err);
    free(kernelbinary);

  if ((!program) || (err!=CL_SUCCESS))
  {
	printf("Error: Failed to create compute program from binary %d!\n", err);
	exit(EXIT_FAILURE);
  }

    // -------------------------------------------------------------
    //  Step 6 for 1 Compute Unit: Create Kernels - the actual handler of the kernel that
	//       	we will be using. We first create a program, and then
	//       	obtain the kernel handler from the program.
    // -------------------------------------------------------------
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
  {
	printf("Error: Failed to build program executable!\n");
	clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
	printf("%s\n", buffer);
	exit(EXIT_FAILURE);
  }
    gibbsSamplerKernel = clCreateKernel(program, "gibbsSamplerKernel", &err);
  if (!gibbsSamplerKernel || err != CL_SUCCESS)
  {
	printf("Error: Failed to create compute gibbsSamplerKernel!\n");
   	 exit(EXIT_FAILURE);
  }
  //// -------------------------------------------------------------
  //  Step 6 for 2 Compute Units : Create Kernels - one handler for each
  //       	instance of the kernel that we will be using. The program
  //       	is the same as the step with 1 compute unit.
  // -------------------------------------------------------------
  // gibbsSamplerKernel_1 = clCreateKernel(program, "gibbsSamplerKernel:{gibbsSamplerKernelInstance_1}", &err);
  // if (!gibbsSamplerKernel_1 || err != CL_SUCCESS) {
  // 	printf("Error: Failed to create compute krnl_gibbsSamplerKernel_1!\n");
  //     exit(EXIT_FAILURE);
  // }
  // err = 0;
  // gibbsSamplerKernel_2 = clCreateKernel(program, "gibbsSamplerKernel:{gibbsSamplerKernelInstance_2}", &err);
  // if (! gibbsSamplerKernel_2 || err != CL_SUCCESS) {
  // 	printf("Error: Failed to create compute gibbsSamplerKernel_2!\n");
  //     exit(EXIT_FAILURE);
  // }


  // -------------------------------------------------------------------------
  // Step 7 : Create buffers.
  // We do not need to allocate separate memory space (malloc), because
  // on a MPSoC system (e.g. ZCU102 board), we map the memory space that is
  // allocated at clCreateBuffer, to a usable memory space for our
  // host application. We also do not need to use free for any reason.
  // See Xilinx UG1393 for detailed information.
  // ------------------------------------------------------------------------

  /****************************
   * Data initialization
   * **************************/
  TICK();
  // Load input data to memory
  //read_input();

  input = fopen("input_brca1.txt", "r+");
  if (input == NULL)
  {
	printf("ERROR: Input not found!\n");
	exit(EXIT_FAILURE);
  }

  fseek(input, 0 , SEEK_SET);
  bytes_read = 0;

  fscanf(input, "%d %d %d", &k, &t, &n);
  printf("k : %d, t : %d, N : %d\n",k, t, n);

  bytes_read = lenHelper(k) + 1 + lenHelper(t) + 1 + lenHelper(n) + 1;
  printf("Bytes read: %d\n", bytes_read);
  fseek(input, bytes_read , SEEK_SET);

  first_motif_buffer = (char*)malloc(5000*sizeof(char));
  check = fgets(first_motif_buffer, 5000*sizeof(char), input);
  if(check == NULL)
	exit(EXIT_FAILURE);

  motif_length = strlen(first_motif_buffer);
  first_motif_buffer[motif_length - 1] = '\0';
  motif_length = strlen(first_motif_buffer);

  printf("Motif length: %d\n", motif_length);
  //printf("first motif: %s\n", first_motif_buffer);

  motifs_1d = (char*)malloc(t*(motif_length+1) * sizeof(char));

  number_of_motifs = 1;
  /*motifs = (char**)malloc(number_of_motifs * sizeof(char*));
  (*motifs) = (char*)malloc((motif_length + 1) * sizeof(char));
  strcpy((*motifs), first_motif_buffer);*/
  strcpy(motifs_1d, first_motif_buffer);
  //printf("%s\n", motifs_1d);
  free(first_motif_buffer);

  while (check != NULL)
  {
	number_of_motifs++;

	check = fgets((motifs_1d + (number_of_motifs - 1)*(motif_length+1)), (motif_length + 1)*sizeof(char), input);
	if((check == NULL) || (*check == '\n'))
	{
  	number_of_motifs--;
  	break;
	}

	(*(motifs_1d + (number_of_motifs - 1)*(motif_length+1) + motif_length)) = '\0';
	//printf("%s\n", (motifs_1d + (number_of_motifs - 1)*(motif_length+1)));
	fseek(input, 1 , SEEK_CUR);
  }

  TOCK("load_time:");
  /*** Input image array buffer ***/
  input_buffer = clCreateBuffer(context,  CL_MEM_READ_ONLY,  number_of_motifs*(motif_length+1), NULL, &err);
  if (err != CL_SUCCESS)
  {
   printf("Return code for clCreateBuffer - input_buffer: %d",err);
  }
    motifs_1d = (char *)clEnqueueMapBuffer(q,input_buffer,CL_TRUE,CL_MAP_WRITE,0,number_of_motifs*(motif_length+1),0,NULL,NULL,&err);

  /*** Output image array buffer ***/
  output_buffer = clCreateBuffer(context,  CL_MEM_READ_WRITE,  t*(k+1)*sizeof(char), NULL, &err);
  if (err != CL_SUCCESS)
  {
   printf("Return code for clCreateBuffer - output_buffer: %d",err);
  }
    best_motifs_output = (char *)clEnqueueMapBuffer(q,output_buffer,CL_TRUE,CL_MAP_READ,0,t*(k+1)*sizeof(char),0,NULL,NULL,&err);

  /*****
   * Kernel execution.
   * Note: We always check for errors.
   * **************************/

  /*****
   * In order to use multiple compute units,
   * consider what arguments earch kernel needs,
   * and how you must split the input and output
   * arrays.
   * ******************************************/
  TICK();

  // Set HW Kernel arguments
  //r = FILTER_RADIUS;
  //size_x = SIZE_X;
  //size_y = SIZE_Y;
  argcounter = 0;
  err = 0;
  err |= clSetKernelArg(gibbsSamplerKernel,argcounter++, sizeof(cl_mem), &input_buffer);
  err |= clSetKernelArg(gibbsSamplerKernel,argcounter++, sizeof(cl_mem), &output_buffer);
  err |= clSetKernelArg(gibbsSamplerKernel,argcounter++, sizeof(int), &k);
  err |= clSetKernelArg(gibbsSamplerKernel,argcounter++, sizeof(int), &t);
  err |= clSetKernelArg(gibbsSamplerKernel,argcounter++, sizeof(int), &n);
  err |= clSetKernelArg(gibbsSamplerKernel,argcounter++, sizeof(int), &motif_length);
  if (err != CL_SUCCESS)
  {
   	 printf("Error: Failed to set gibbsSamplerKernel arguments! %d\n", err);
	 }

  // Enqueue input memory objects migration - Host -> Device
  pt_in[0] = input_buffer;
  pt_out[0] = output_buffer;

  err = clEnqueueMigrateMemObjects(q,(cl_uint)1, pt_in, 0 ,0,NULL, NULL);
  if (err)
  {
	printf("Error: Failed to migrate memobjects to device! %d\n", err);
	return EXIT_FAILURE;
  }

  // Start kernel execution
    err = clEnqueueTask(q, gibbsSamplerKernel, 0, NULL, NULL);
  if (err)
  {
	printf("Error: Failed to execute kernel! %d\n", err);
	return EXIT_FAILURE;
  }

  // Enqueue output memory objects migration - Device -> Host
    err = clEnqueueMigrateMemObjects(q,(cl_uint)1, pt_out, CL_MIGRATE_MEM_OBJECT_HOST,0,NULL, NULL);
    if (err != CL_SUCCESS)
  {
	printf("Error: Failed to migrate membojects from device: %d!\n", err);
	return EXIT_FAILURE;
  }

  // Wait for execution to finish
  clFinish(q);
  TOCK("gibbs_sampler time:");

  printf("Best Motifs 1d: \n");
  for(i = 0 ; i < t ; i++)
  {
	printf("%s\n", best_motifs_output+i*(k+1));
  }
  // Compare output results with golden
  //TICK();
  //compare();
  //TOCK("compare_time:");


  /*****
   * Clean up code.
   * We always free any memory allocated, and release
   * all OpenCL objects.
   * ***********************************************/

  // Release memory objects and other OpenCL objects
  err = clEnqueueUnmapMemObject(q,input_buffer,motifs_1d,0,NULL,NULL);
    if(err != CL_SUCCESS)
  {
   	 printf("Error: Failed to unmap device memory input!\n");
    }

  err = clEnqueueUnmapMemObject(q,output_buffer,best_motifs_output,0,NULL,NULL);
    if(err != CL_SUCCESS)
  {
   	 printf("Error: Failed to unmap device memory output!\n");
    }

  /*err = clEnqueueUnmapMemObject(q,gaussian_buffer,gaussian,0,NULL,NULL);
    if(err != CL_SUCCESS)
  {
   	 printf("Error: Failed to unmap device memory gaussian!\n");
    }*/
  clFinish(q);


  clReleaseMemObject(input_buffer);
    clReleaseMemObject(output_buffer);
    //clReleaseMemObject(gaussian_buffer);

    clReleaseProgram(program);
  clReleaseKernel(gibbsSamplerKernel);
  clReleaseCommandQueue(q);
  clReleaseContext(context);

  return 0;

}
