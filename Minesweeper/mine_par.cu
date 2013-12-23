#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>

#include <map>
#include <vector>
#include <iostream>
#include <cmath>

#include <stdlib.h>     /* srand, rand */
#include <time.h> 

#include "CycleTimer.h"

#define IDX(x,y) ((x) + (width * (y)))
#define GRID(x,y) grid[((x) + (width * (y)))]
#define BOMB(x,y) (GRID((x),(y)) == 10 || GRID((x),(y)) == (-10))
#define HIDDEN(x,y) (GRID((x),(y)) < 0)
#define IX(i) ((i)%width)
#define IY(i) ((i)/width)


void printGrid(int width, int height, int *grid) {
  int i,j;
  std::cout<<"   ";
  for(i=0; i<width; i++)
    std::cout<<(i%10)<<"   ";
  std::cout<<"\n";

  for(j=0; j<height; j++) {
    std::cout<<(j%10)<<" ";
    for(i=0; i<width; i++) {
      if(grid[IDX(i,j)]==0){
	std::cout<<"    ";
        continue;
      }
      std::cout<<grid[IDX(i,j)]<<" ";
      if(grid[IDX(i,j)]>=0)
	std::cout<<" ";
      if(abs(grid[IDX(i,j)])<10)
	std::cout<<" ";
    }
    std::cout<<"\n";
  }
  std::cout<<"   ";
  for(i=0; i<width; i++)
    std::cout<<(i%10)<<"   ";
  std::cout<<"\n";
}

void printGridCompact(int width, int height, int *grid) {
  int i,j;
  std::cout<<"  ";
  for(i=0; i<width; i++)
    std::cout<<(i%10)<<" ";
  std::cout<<"\n";

  for(j=0; j<height; j++) {
    std::cout<<(j%10)<<" ";
    for(i=0; i<width; i++) {
      if(grid[IDX(i,j)]==0){
	std::cout<<"  ";
        continue;
      }
      std::cout<<grid[IDX(i,j)]<<" ";
    }
    std::cout<<"\n";
  }
  std::cout<<"  ";
  for(i=0; i<width; i++)
    std::cout<<(i%10)<<"   ";
  std::cout<<"\n";
}


/*
__global__ void
saxpy_kernel(int N, float alpha, float* x, float* y, float* result) {
  
  // compute overall index from position of thread in current block,
  // and given the block we are in
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  
  if (index < N)
    result[index] = alpha * x[index] + y[index];
}
*/

/**
 * If on a border will set to index number, 0 otherwise. 
 * If the cell is on the grid border but not on a real border will be negative
 */
__global__ void identify_boundaries(int width, int height, int *grid, int *borders) {
  //int index = blockIdx.x * blockDim.x + threadIdx.x;
  //int index = threadIdx.x;
  //int x = IX(index), y = IY(index);
  int x = threadIdx.x, y = blockIdx.x;
  int index = x + (y*width);
  int i,j;

  borders[index] = 0;

  if(!HIDDEN(x,y)) 
    return;

  if(x==0 || x==(width-1) || y==0 || y==(height-1)) {
    borders[index] = -(index + 1);
    //return;
  }

  for(i=max(0, x-1); i<=min(width-1, x+1); i++) {
    for(j=max(0, y-1); j<=min(height-1, y+1) ; j++) {
      if(!HIDDEN(i,j)) {
	borders[index] = index + 1;
	return;
      }
    }
  }

}

/**
 * If a tile is on a boundary it will attempty to make itself a smaller boundary
 *   by looking at it's neighbors. If it finds a smaller neighbor then set done to true.
 * If the tile is on the grid border but not the solving border then the value will be negative.
 */
__global__ void consolidate_boundaries(int width, int height, 
				       int *grid, int *borders, bool *done) {
  //int index = threadIdx.x;
  //int x = IX(index), y = IY(index);
  int x = threadIdx.x, y = blockIdx.x;
  int index = x + (y*width);

  int i,j;
  int newMin = borders[index];
  int other;

  bool border = newMin < 0;
  if(border)
    newMin = -newMin;

  for(j=max(0, y-2); j<=min(height-1, y+2); j++) {
    for(i=max(0, x-2); i<=min(width-1, x+2); i++) {
      other = abs(borders[IDX(i,j)]);
      if(other != 0 && other<newMin)
	newMin = other;
    }
  }

  if(border)
    newMin = -newMin;

  if(newMin != borders[index]) {
    borders[index] = newMin;
    *done = false;
  }
}

/**
 * Every revealed cell next to a border cell will become negative of that border id.
 * All other cells will become 0 (except other border cells). 
 * When called *border should have positive values on borders, negative on grid
 *   borders, and 0 otherwise.
 */
__global__ void update_borders(int width, int height, int *border, int *grid, bool *flags) {
  //int index = threadIdx.x;
  //int x = IX(index), y = IY(index);
  int x = threadIdx.x, y = blockIdx.x;
  int index = x + (y*width);

  int i, j, xt, yt;

  // if a cell is flagged then automatically becomes 0
  if(flags[index]) {
    border[index] = 0;
    return;
  }

  // if on a border do nothing
  if(border[index] > 0) {
    return;
  }

  // if a cell is negative then it was on a grid border and not an actual border
  // if hidden then don't do anything with the cell
  if(border[index]<0 || HIDDEN(x,y)) {
    border[index] = 0;
    return;
  }

  // each known cell tries to identify with a border
  for(j=-1; j<=1; j++) {
    yt = y+j;
    if(yt<0 || yt>=height)
      continue;
    for(i=-1; i<=1; i++) {
      xt = x+i;
      if(xt<0 || xt>=width)
	continue;
      if(border[IDX(xt,yt)]>0 && !flags[IDX(xt,yt)]) {
	border[index] = -border[IDX(xt,yt)];
	return;
      }
    }
  }  
}

/**
 * sets all border cells to 0
 */
__global__ void grid_borders_zero(int width, int height, int *arr) {
  //int index = threadIdx.x;
  int index = threadIdx.x + (width*blockIdx.x);
  arr[index]=0;
}

/**
 * Sets all cells that are 3 to 0
 */
__global__ void finalize_knowns(int width, int height, int *knowns) {
  //int index = threadIdx.x;
  int index = threadIdx.x + (width*blockIdx.x);
  if(knowns[index]==3)
    knowns[index]=0;
}


/**
 * Sets all cells that are uncovered or flagged to 0
 */
__global__ void clean_knowns(int width, int height, int *grid, bool *flags, int *knowns) {
  //int index = threadIdx.x;
  int index = threadIdx.x + (width*blockIdx.x);
  if(flags[index] || grid[index]>=0)
    knowns[index]=0;
}

/**
 * Fill in the temp grid as follows: positive number of unknown bombs if uncovered. 
 * Subtract known (flagged) bombs from self (if positive number). 
 * 0 otherwise
 */
__global__ void make_temp_grid(int width, int height, int *grid, bool *flags, 
			       int *old_knowns, int *temp_grid) {
  //int index = threadIdx.x;
  //int x = IX(index), y = IY(index);
  int x = threadIdx.x, y = blockIdx.x;
  int index = x + (y*width);

  int i, j, xt, yt;

  if(flags[index] || grid[index]<0) {
    temp_grid[index] = 0;
    return;
  }

  int grid_num = grid[index];

  for(j=-1; j<=1; j++) {
    yt = y+j;
    if(yt<0 || yt>=height)
      continue;
    for(i=-1; i<=1; i++) {
      xt = x+i;
      if(xt<0 || xt>=width)
	continue;
      if(flags[IDX(xt,yt)] || old_knowns[IDX(xt,yt)]==1)
	grid_num--;
    }
  }

  temp_grid[index] = grid_num;
}

/**
 * Uses values in blockId and threadId to know where bombs are placed
 *   in this call. The combined 4 ints can be seen as one long bit vector
 *   which store where bombs are, index will return true if there is a 1 (bomb)
 *   at that index
 */
__device__ bool is_bomb(int index, int thread_max, int x_block_max, int y_block_max) {
  // only need threadIdx, if index<10 look at the first 1024 values which are stored
  //   in threadIdx.x
  if(index < thread_max) {
    return ((threadIdx.x>>index) & 1) == 1;
  }

  // there's probably a cleaner method to do this, but its pretty short and not necessarily slow
  index -= thread_max;
  if(index < x_block_max) {
    return ((blockIdx.x>>index) & 1) == 1;
  }

  index -= x_block_max;
  if(index < y_block_max) {
    return ((blockIdx.y>>index) & 1) == 1;
  }

  index -= y_block_max;
  //if(index < max_block_dim) {
    return ((blockIdx.z>>index) & 1) == 1;
    //}
    //return true;
}

/**
 * Each thread takes one possible solution, as defined by threadIdx and blockIdx. 
 * If the solution is not valid do nothing. 
 * All values in knowns should be 0 at start of the call. 
 * If the value can be a bomb the first bit will be 1, else 0. 
 * If the value can be empty the second bit will be 1, else 0. 
 */
__global__ void find_solvable_tiles(int width, int height, int *grid, int *knowns, 
				    int *border, int border_size, int *temp_grid, 
				    int *known_border, int known_border_size, 
				    int thread_max, int x_block_max, int y_block_max) {
  // in the for loop keeps track of which known cell we are checking against
  int c, i, j, x, y, cell, bombs;

  for(c = 0; c<known_border_size; c++) { 
    cell = known_border[c];
    bombs = 0;
    x = IX(cell);
    y = IY(cell);

    for(i=max(x-1,0); i<=min(x+1,width-1); i++) {
      for(j=max(y-1,0); j<=min(y+1,height-1); j++) {
	if(temp_grid[IDX(i,j)]<0)
	  if(is_bomb((-1-temp_grid[IDX(i,j)]), thread_max, x_block_max, y_block_max))
	    bombs++;
      }
    }

    if(temp_grid[cell]!=bombs)
      return;
  }

  // if the program gets here then this is a valid solution
  // commit into knowns (first bit->1 if bomb, second bit->1 if not)
  for(c=0; c<border_size; c++) {
    cell = border[c];
    //atomicAdd(&knowns[cell], 1);
    if(is_bomb(c, thread_max, x_block_max, y_block_max)) {
      atomicOr(&knowns[cell], 1);
    } else {
      atomicOr(&knowns[cell], 2);
    }
  }
}

/**
 * Scans ever element in the border if it has solution
 * Callable right after find_solvable_tiles with same data
 */
__global__ void border_has_solution(int width, int *knowns, int *border, 
				    int border_id, bool *solution) {
  //int index = threadIdx.x;
  int index = threadIdx.x + (width * blockIdx.x);
  if(border[index] != border_id)
    return;
  int c = knowns[index];
  if(c==1 || c==2)
    *solution = true;
}

__global__ void preliminary_knowns(int width, int height, int *grid, bool *flags, 
				   int *border, int *knowns) {
  //int index = threadIdx.x;
  int x = threadIdx.x, y = blockIdx.x;
  int index = x + (y*width);

  if(index >= (width*height) || border[index]>=0)
    return;
  //int x = IX(index), y = IY(index);
  int bomb_count = 0, unknown_count = 0;

  for(int j=max(0,y-1); j<=min(height-1, y+1); j++) 
    for(int i=max(0,x-1); i<=min(width-1, x+1); i++) {
      int c = IDX(i,j);
      if(knowns[c]==1 || flags[c])
	bomb_count++;
      else if(knowns[c]==0 && grid[c]<0) 
	unknown_count++;
    }

  if(unknown_count == 0)
    return;

  int all_bombs = 0; // will be 1 if all should be bombs, -1 if clear, 0 ow
  if((grid[index]-bomb_count) == unknown_count)
    all_bombs = 1;
  else if(bomb_count == grid[index])
    all_bombs = -1;
  
  // if know neighbors then fill them
  if(all_bombs!=0) {
    for(int j=max(0,y-1); j<=min(height-1, y+1); j++) 
      for(int i=max(0,x-1); i<=min(width-1, x+1); i++) {
	int c = IDX(i,j);
	if(flags[c] || grid[c]>=0 || border[c]<=0 || knowns[c]!=0)
	  continue;
	if(all_bombs<0)//all clear
	  atomicOr(&knowns[c], 2);
	else //all bombs
	  atomicOr(&knowns[c], 1);
      }
  }
}

/**
 * Recursive call to display cell
 */
void hint_display_cell(int x, int y, int width, int height, int *grid) {
  if(GRID(x,y)>=0)
    return;

  grid[IDX(x,y)] = -(GRID(x,y)+1);

  if(GRID(x,y)==0) {
    for(int i = std::max(x-1,0); i<=x+1 && i<width; i++) {
      for(int j = std::max(y-1,0); j<=y+1 && j<height; j++) {
        if(HIDDEN(i,j))
          hint_display_cell(i,j, width, height, grid);
      }
    }
  }
}

/**
 * Knowns is the working array that this function will overwrite with new data. 
 * Assumes that all flags are correct. 
 */
bool known_cells(int width, int height, int *grid, bool *flags, int *knowns, 
		 int* border, bool random_sol, bool display_time, bool display_knowns,
		 bool display_debug_grids) {
  double totalStartTime = CycleTimer::currentSeconds();

  cudaDeviceProp *prop = (cudaDeviceProp *)malloc(sizeof(cudaDeviceProp));
  cudaGetDeviceProperties(prop, 0);
  //std::cout<<"Threads: "<<prop->maxThreadsPerBlock<<" Blocks: "<<prop->maxGridSize[0]<<
  //" "<<prop->maxGridSize[1]<<" "<<prop->maxGridSize[2]<<"\n";

  bool changed_grid = false;

  int size = width*height;

  //int *border = new int[size];
  int *device_border;

  int *device_grid;
  bool *device_flags;
  int *device_knowns;
  int *device_temp_grid;

  bool *device_border_done;
  bool *device_solution;

  int bool_array_size = sizeof(bool)*size;
  int int_array_size = sizeof(int)*size;

  cudaMalloc(&device_grid, int_array_size);
  cudaMalloc(&device_flags, bool_array_size);
  cudaMalloc(&device_knowns, int_array_size);
  cudaMalloc(&device_temp_grid, int_array_size);
  cudaMalloc(&device_border, int_array_size);

  cudaMalloc(&device_border_done, sizeof(bool));
  cudaMalloc(&device_solution, sizeof(bool));

  cudaMemcpy(device_grid, grid, int_array_size, cudaMemcpyHostToDevice);
  cudaMemcpy(device_knowns, knowns, int_array_size, cudaMemcpyHostToDevice);
  cudaMemcpy(device_flags, flags, bool_array_size, cudaMemcpyHostToDevice);

  // each cell identifies itself as border or grid boundary
  identify_boundaries<<<height, width>>>(width, height, device_grid, device_border);
  cudaThreadSynchronize();
  
  bool border_done = false;

  // make all numbers on the same border the same (lowest) number
  while(!border_done) {
    border_done = true;
    cudaMemcpy(device_border_done, &border_done, sizeof(bool), cudaMemcpyHostToDevice);
    consolidate_boundaries<<<height, width>>>(width, height, device_grid, 
					device_border, device_border_done);
    cudaThreadSynchronize();
    cudaMemcpy(&border_done, device_border_done, 
	       sizeof(bool), cudaMemcpyDeviceToHost);
  }

  //cudaThreadSynchronize();
  update_borders<<<height, width>>>(width, height, device_border, device_grid, device_flags);
  cudaThreadSynchronize();
  // borders are now identified by numbers
  cudaMemcpy(border, device_border, int_array_size, cudaMemcpyDeviceToHost);

  if(display_debug_grids) {
    std::cout<<"\nBorder:\n";
    printGrid(width, height, border);
  }

  // start filling out temp grid
  // uncovered cells in device_temp_grid will now have the number of unflagged mines touching them
  clean_knowns<<<height, width>>>(width, height, device_grid, device_flags, device_knowns);
  cudaThreadSynchronize();


  double prelimStartTime = CycleTimer::currentSeconds();
  preliminary_knowns<<<height, width>>>(width, height, device_grid, device_flags, 
				  device_border, device_knowns);
  cudaThreadSynchronize();
  preliminary_knowns<<<height, width>>>(width, height, device_grid, device_flags, 
					device_border, device_knowns);
  cudaThreadSynchronize();
  double prelimEndTime = CycleTimer::currentSeconds();

  if(display_debug_grids) {
    cudaMemcpy(knowns, device_knowns, int_array_size, cudaMemcpyDeviceToHost);
    std::cout<<"After prelim:\n";
    printGrid(width, height, knowns);
  }

  make_temp_grid<<<height, width>>>(width, height, device_grid, device_flags, 
			      device_knowns, device_temp_grid);
  cudaThreadSynchronize();
  
  cudaMemcpy(knowns, device_knowns, int_array_size, cudaMemcpyDeviceToHost);

  int *temp_grid = new int[size];
  // get the temp grid from CPU
  cudaMemcpy(temp_grid, device_temp_grid, int_array_size, cudaMemcpyDeviceToHost);

  std::map<int, std::vector<int> > borders;
  std::map<int, std::vector<int> > known_borders;

  // fill in the lists for the borders
  // One list for the known and another for the unknown borders
  for(int i=0; i<size; i++) {
    // if border is 0 or there is a previous solution, don't be part of the border
    if(border[i]==0 || knowns[i]!=0)
      continue;
    // fill in vector for the known boundary
    if(border[i]<0) {
      if(known_borders.count(-border[i])==0)
	known_borders[-border[i]] = std::vector<int>(1,i);
      else
	known_borders[-border[i]].push_back(i);
      continue;
    }

    if(borders.count(border[i])==0)
      borders[border[i]] = std::vector<int>(1,i);
    else
      borders[border[i]].push_back(i);
    // temp grid will have the -index-1 when the cell has an unknown, will index into border array
    temp_grid[i] = -(borders[border[i]].size());
  }

  // changed temp_grid above, now put back into device_temp and delete
  cudaMemcpy(device_temp_grid, temp_grid, int_array_size, cudaMemcpyHostToDevice);
  if(display_debug_grids) {
    std::cout<<"\nTemp Grid:\n";
    printGrid(width, height, temp_grid);
  }
  delete(temp_grid);

  // Calculate the max number of cells each dimension can be responsible for.
  // Do this by calculating log_2 of the max thread and difference max dimensions
  //   based on GPU information
  int thread_max = 0, x_block_max = 0, y_block_max = 0, z_block_max = 0;
  int temp = prop->maxThreadsPerBlock;
  while (temp >>= 1) ++thread_max;
  temp = prop->maxGridSize[0];
  while (temp >>= 1) ++x_block_max;
  temp = prop->maxGridSize[1];
  while (temp >>= 1) ++y_block_max;
  temp = prop->maxGridSize[2];
  while (temp >>= 1) ++z_block_max;
  z_block_max--;
  x_block_max--;
  y_block_max--;
  if(display_time) 
    std::cout<<"Max Threads: "<<thread_max<<" Max Blocks: "<<x_block_max<<" "<<
      y_block_max<<" "<<z_block_max<<"\n";
  //std::cout<<prop->maxThreadsPerBlock<<" "<<prop->maxGridSize[0]<<" "<<
  //prop->maxGridSize[1]<<" "<<prop->maxGridSize[2]<<"\n";

  double startTime = CycleTimer::currentSeconds();
  double calcTime = 0;

  // for each boundary find possible solutions, actuall work here
  for (std::map<int,std::vector<int> >::iterator it=borders.begin(); 
       it!=borders.end(); ++it) {
    std::vector<int> border_vec = it->second;
    std::vector<int> known_border_vec = known_borders[it->first];

    /*
    std::cout<<"Border "<< it->first <<" size: "<<border_vec.size()<<
      " known size: "<<known_border_vec.size()<<"\nBorder: ";
    for(int i=0; i<border_vec.size(); i++)
      std::cout<<border_vec[i]<<" ";
    std::cout<<"\nKnown: ";
    for(int i=0; i<known_border_vec.size(); i++)
      std::cout<<known_border_vec[i]<<" ";
    std::cout<<"\n";
    */

    //std::cout << it->first << " => " << it->second << '\n';
    int *border_elements;
    cudaMalloc(&border_elements, sizeof(int)*border_vec.size());
    cudaMemcpy(border_elements, &border_vec[0], sizeof(int)*border_vec.size(), 
	       cudaMemcpyHostToDevice);

    int *known_border_elements;
    cudaMalloc(&known_border_elements, sizeof(int)*known_border_vec.size());
    cudaMemcpy(known_border_elements, &known_border_vec[0], 
	       sizeof(int)*known_border_vec.size(), cudaMemcpyHostToDevice);

    // these are the number of grid cells the different dimensions will represent
    // actual numbers that will be sent to the gpu will be 2^x
    int block_threads = thread_max;
    int block_x = 0;
    int block_y = 0;
    int block_z = 0;

    int border_left = border_vec.size();
    bool skip = false;

    if(border_left <= thread_max) {
      block_threads = border_left;
    } else {
      border_left -= thread_max;
      
      // set block_x
      if(border_left <= x_block_max) {
	block_x = border_left;
      } else {
	block_x = x_block_max;
	border_left -= x_block_max;

	// if still more left set block_y
	if(border_left <= y_block_max) {
	  block_y = border_left;
	} else {
	  block_y = y_block_max;
	  border_left -= y_block_max;
	  
	  // if still more left set block_z
	  if(border_left <= z_block_max) {
	    block_z = border_left;
	  } else {
	    //if you get to this point the border is quite large, just
	    // give up
	    skip = true;
	  }
	}
      }
    }
    
    if(!skip) {
      dim3 grid_block(1<<block_x, 1<<block_y, 1<<block_z);
      
      if(display_time)
		std::cout<<"Border size: "<<border_vec.size()<<"\nThreads: "<<
	  pow(2.0,block_threads)<< " Grid: "<<grid_block.x<<" "<<
	  grid_block.y<<" "<<grid_block.z<<"\n";

      double calcStartTime = CycleTimer::currentSeconds();
      find_solvable_tiles<<< grid_block,
		1<<block_threads >>>(width, height, device_grid, device_knowns, 
				  border_elements, border_vec.size(), device_temp_grid,
				  known_border_elements, known_border_vec.size(), 
				  thread_max, x_block_max, y_block_max);
      cudaThreadSynchronize();
      double calcEndTime = CycleTimer::currentSeconds();
      if(display_time)
	std::cout<<"Calc Time: "<<(calcEndTime-calcStartTime)<<"\n\n";
      calcTime += (calcEndTime-calcStartTime);

      // check if there is a solution for this border, if not unveil a random tile
      //random_sol = false;
      if(random_sol) {
	bool solution=false;
	cudaMemcpy(device_solution, &solution, sizeof(bool), cudaMemcpyHostToDevice);
	border_has_solution<<<height, width>>>(width, device_knowns, device_border, it->first, 
					       device_solution);
	cudaThreadSynchronize();
	cudaMemcpy(&solution, device_solution, sizeof(bool), cudaMemcpyDeviceToHost);
	if(!solution) {
	  //std::cout<<"NO SOLUTION for border "<< it->first <<"\n";
	  changed_grid = true;
	}
	srand(time(NULL));
	int c = 0;//don't want to repeat forever
	while(!solution && c<border_vec.size()) {
	  int t = border_vec[rand()%border_vec.size()];
	  int x = IX(t), y = IY(t);
	  if(!BOMB(x,y)) {
	    std::cout<<"Revealed grid "<<x<<" "<<y<<"\n";
	    hint_display_cell(x, y, width, height, grid);
	    solution = true;
	  }
	  c++;
	}
      }
    }

    // clean up
    cudaFree(border_elements);
    cudaFree(known_border_elements);
  }

  if(display_time) {
    std::cout<<"\nTotal Calc Time: "<<calcTime<<"\n";
    double endTime = CycleTimer::currentSeconds();
    std::cout<<"Prelim run time: "<<(prelimEndTime-prelimStartTime)<<
      "\nMain work loop: "<<(endTime-startTime)<<"\n";
  }
  //cudaThreadSynchronize();

  //cudaMemcpy(knowns, device_knowns, int_array_size, cudaMemcpyDeviceToHost);
  //std::cout<<"\nBefore Clean:\n";
  //printGrid(width,height,knowns);

  finalize_knowns<<<height, width>>>(width, height, device_knowns);
  clean_knowns<<<height, width>>>(width, height, device_grid, device_flags, device_knowns);
  cudaThreadSynchronize();

  cudaMemcpy(knowns, device_knowns, int_array_size, cudaMemcpyDeviceToHost);
  
  cudaFree(device_grid);
  cudaFree(device_flags);
  cudaFree(device_knowns);
  cudaFree(device_temp_grid);
  cudaFree(device_border_done);
  cudaFree(device_border);
  cudaFree(device_solution);

  free(prop);
  delete(border);

  if(display_time) {
    double totalEndTime = CycleTimer::currentSeconds();
    std::cout<<"Total execution time: "<<(totalEndTime-totalStartTime)<<"\n";
  }

  if(display_knowns || display_debug_grids) {
    std::cout<<"\nFinal knowns:\n";
    printGridCompact(width,height,knowns);
    //std::cout<<"Done\n";
  }

  return changed_grid;
}



void
saxpyCuda(int N, float alpha, float* xarray, float* yarray, float* resultarray) {
  /*
    int totalBytes = sizeof(float) * 3 * N;

    // compute number of blocks and threads per block
    const int threadsPerBlock = 512;
    const int blocks = (N + threadsPerBlock - 1) / threadsPerBlock;

    float* device_x;
    float* device_y;
    float* device_result;

    //
    // TODO allocate device memory buffers on the GPU using cudaMalloc
    //
    size_t size = sizeof(float)*N;
    cudaMalloc(&device_x, size);
    cudaMalloc(&device_y, size);
    cudaMalloc(&device_result, size);

    //
    // TODO copy input arrays to the GPU using cudaMemcpy
    //
    cudaMemcpy(device_x, xarray, size, cudaMemcpyHostToDevice);
    cudaMemcpy(device_y, yarray, size, cudaMemcpyHostToDevice);

    // start timing after allocation of device memory
    double startTime = CycleTimer::currentSeconds();

    // run kernel
    saxpy_kernel<<<blocks, threadsPerBlock>>>(N, alpha, device_x, device_y, device_result);
    cudaThreadSynchronize();

    // end timing after result has been copied back into host memory

    double endTime = CycleTimer::currentSeconds();

    //
    // TODO copy result from GPU using cudaMemcpy
    //
    cudaMemcpy(resultarray, device_result, size, cudaMemcpyDeviceToHost);

    cudaError_t errCode = cudaPeekAtLastError();
    if (errCode != cudaSuccess) {
        fprintf(stderr, "WARNING: A CUDA error occured: code=%d, %s\n", errCode, cudaGetErrorString(errCode));
    }

    double overallDuration = endTime - startTime;
    printf("Overall: %.3f ms\t\t[%.3f GB/s]\n", 1000.f * overallDuration, toBW(totalBytes, overallDuration));

    // TODO free memory buffers on the GPU
    cudaFree(device_x);
    cudaFree(device_y);
    cudaFree(device_result);
  */
}

void
printCudaInfo() {

    // for fun, just print out some stats on the machine

    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);

    printf("Found %d CUDA devices\n", deviceCount);

    for (int i=0; i<deviceCount; i++) {
        cudaDeviceProp deviceProps;
        cudaGetDeviceProperties(&deviceProps, i);
        printf("Device %d: %s\n", i, deviceProps.name);
        printf("   SMs:        %d\n", deviceProps.multiProcessorCount);
        printf("   Global mem: %.0f MB\n",
               static_cast<float>(deviceProps.totalGlobalMem) / (1024 * 1024));
        printf("   CUDA Cap:   %d.%d\n", deviceProps.major, deviceProps.minor);
    }
}
