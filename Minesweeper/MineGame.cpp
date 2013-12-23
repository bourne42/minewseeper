
#include <iostream>
#include <stdlib.h> // for rand
#include <time.h>
#include <algorithm> // std::max

#include "MineGame.h"

#define IDX(x,y) ((x) + (width * (y)))
#define GRID(x,y) grid[((x) + (width * (y)))]
#define FLAG(x,y) flags[((x) + (width * (y)))]
#define GRIDN(x,y) (-GRID(x,y) - 1)
#define BOMB(x,y) (GRID((x),(y)) == 10 || GRID((x),(y)) == (-10))
#define HIDDEN(x,y) (GRID((x),(y)) < 0)

MineGame::MineGame(int width, int height, int num_mines, bool game_to_file, 
		   bool display_time, bool display_knowns, bool display_debug_grids) : 
  width(width), height(height), num_mines(num_mines), display_time(display_time), 
  display_knowns(display_knowns), display_debug_grids(display_debug_grids)
{
  grid = new int[width*height];
  flags = new bool[width*height];
  knownCells = new int[width*height];

  //border = new int[width*height];

  //has_game_to_write = false;

  //blah = new int[width*height]();
  //asdf = new int[width*height]();
  restart();
  //t = new int[width*height];
  //std::cout<<t<<"\n";
}

MineGame::~MineGame(void)
{
}
/*
void MineGame::complete_write() {
  out_file.close();
  write_game = false;
  }*/

void MineGame::restart() {
  /*if(has_game_to_write && write_game)
    complete_write();
  if(write_game) {
    std::time_t t = std::time();
    std::cout<<"New file: "<<asctime(t)<<"\n";
    out_file.open(asctime(t));
    }*/

  for(int i=0; i<width*height; i++) {
    flags[i] = false;
    knownCells[i] = 0;
  }
  setup_grid();
  finished = 0;
  incorrect = 0;
  first_move = true;
}

void MineGame::setup_grid()
{
  int c = 0, x, y, i;//, j;
  srand(time(NULL)); // initialize random seed

  //set all to 0
  for(y = 0; y<height; y++) {
    for(x = 0; x<width; x++) {
      grid[IDX(x,y)] = 0;
    }
  }

  // try to place bombs in random places
  for(i=0; i<num_mines*2 && c<num_mines; i++){
    x = rand() % width;
    y = rand() % height;
    if(!BOMB(x,y)) {
      GRID(x,y)=-10;
      c++;
    }
  }

  // There are more effeciant ways to do this (avoid some recomputation).
  // I will leave as is because more important things to improce and this
  //   is only done once.
  /*for(y = 0; y<height; y++) {
    for(x = 0; x<width; x++) {
      if(BOMB(x,y))
	continue;
      for(i = std::max(x-1,0); i<=x+1 && i<width; i++) {
	for(j = std::max(y-1,0); j<=y+1 && j<height; j++) {
	  if(BOMB(i,j))
	    GRID(x,y)++;
	}
      }
      GRID(x,y) = -(GRID(x,y)+1);
    }
    }*/
  calculate_neighbors();

  /*
  if(write_game) {
    out_file<<width<<" "<<height<<"\n";
    
    for(y = 0; y<height; y++) {
      for(x = 0; x<width; x++) {
	std::cout<<x;
	if(x!=(width-1))
	  std::cout<<" ";
      }
      std::cout<<"\n";
    }
  }
  */
}

void MineGame::calculate_neighbors() {
  int x, y, i, j;
  for(y = 0; y<height; y++) {
    for(x = 0; x<width; x++) {
      if(BOMB(x,y))
	continue;
      grid[IDX(x,y)] = 0;
      for(i = std::max(x-1,0); i<=x+1 && i<width; i++) {
	for(j = std::max(y-1,0); j<=y+1 && j<height; j++) {
	  if(BOMB(i,j))
	    GRID(x,y)++;
	}
      }
      grid[IDX(x,y)] = -(GRID(x,y)+1);
    }
  }
}

bool MineGame::check_game_won() {
  if(finished!=0)
    return finished == 1;

  int x,y;
  
  for(y = 0; y<height; y++) {
    for(x = 0; x<width; x++) {
      if(!HIDDEN(x,y) && BOMB(x,y)){
	finished = -1;
	return false;
      }
      if(HIDDEN(x,y) && !BOMB(x,y))
	return false;
    }
  }

  finished = 1;
  return true;
}

void MineGame::press_cell(int x, int y, bool flag) {
  if(finished!=0)
    return;

  /*has_game_to_write = true;
  if(write_game) {
    out_file<<x<<" "<<y;
    if(flag)
      out_file<<" f";
    out_file<<"\n";
    }*/

  // don't do anything if the cell has been unveiled
  // if it has been unveiled and all its bombs are accounted for then clear all neighbors
  if(!HIDDEN(x,y)) {
    bool done = true;
    // check if flagged all bombs
    if(GRID(x,y)>0) {
      int bomb_count = 0;
      for(int i=std::max(0,x-1); i<=std::min(width-1, x+1); i++) 
	for(int j=std::max(0,y-1); j<=std::min(height-1, y+1); j++) {
	  if(FLAG(i,j))
	    bomb_count++;
	}
      // all bombs accounted for, clear rest
      if(GRID(x,y)==bomb_count) {
	done = false;
	for(int i=std::max(0,x-1); i<=std::min(width-1, x+1); i++) 
	  for(int j=std::max(0,y-1); j<=std::min(height-1, y+1); j++) {
	    if(!FLAG(i,j))
	      display_cell(i,j);
	  }
      }
    }

    if(done) {
      std::cout<<"Square is already uncovered\n";
      return;
    }
  }

  if(first_move) {
    first_move = false;
    if(BOMB(x,y)) { // move bomb
      int c = 0;
      while(c<=20) {
	int xt = rand() % width;
	int yt = rand() % height;
	if(!BOMB(xt,yt)) {
	  GRID(xt,yt)=-10;
	  c = 21;
	}
      }
      grid[IDX(x,y)]=0;
      calculate_neighbors();
      //std::cout<<"Moved bomb, your welcome\n";
    }
  }
  
  if(flag) {
    if(!BOMB(x,y) && FLAG(x,y))
      incorrect--;
    else if(!BOMB(x,y) && !FLAG(x,y)) {
      std::cout<<"WRONG FLAG\n";
      incorrect++;
    }
    flags[IDX(x,y)] = !FLAG(x,y);
  }
  else
    display_cell(x,y);
  
  if(check_game_won()) {
    if(finished == 1)
      std::cout<<"Won, restarting...\n";
    else
      std::cout<<"Lost, restarting...\n";
    restart();
    return;
  }

  // don't recompute on a flag because will have no new data
  // not necessary because small time saved, but will avoid unnecessary work
  if(incorrect==0){
    bool asdf = true;
    int c = 0;
    while(asdf && c<20) {
    int *border = new int[width*height];
      /*
      int c = 0;
      while(known_cells(width, height, grid, flags, knownCells, border, true) && c<2){
	std::cout<<"Try again... ";
	c++;
      }
      */
    if(!known_cells(width, height, grid, flags, knownCells, border, true, 
		    display_time, display_knowns, display_debug_grids))
      asdf = false;
    else {
      std::cout<<"Grid try "<<c<<"\n";
      print_grid();
    }
    //if(c!=0)
    //std::cout<<"Tried for solution\n";
    c++;
    //break;
    //std::cout<<"\nTried "<<c<<" times for a solution\n";
    }
  }
}

void MineGame::display_cell(int x, int y) {
  int i,j;

  if(GRID(x,y)>=0)
    return;

  if(BOMB(x,y)) {
    grid[IDX(x,y)] = 10;
    lose_game();
    return;
  }

  grid[IDX(x,y)] = -(GRID(x,y)+1);
  
  if(GRID(x,y)==0) {
    for(i = std::max(x-1,0); i<=x+1 && i<width; i++) {
      for(j = std::max(y-1,0); j<=y+1 && j<height; j++) {
	if(HIDDEN(i,j))
	  display_cell(i,j);
      }
    }
  }
}

void MineGame::lose_game() {
  finished = -1;
  //std::cout<<"Lost, restarting...\n";
  //restart();
}

int MineGame::game_finished() {
  return finished;
}

int MineGame::getHeight() {
  return height;
}

int MineGame::getWidth() {
  return width;
}

bool MineGame::flag(int x, int y) {
  return FLAG(x,y);
}

int MineGame::grid_val(int x, int y) {
  return GRID(x,y);
}

int MineGame::known_val(int x, int y) {
  return knownCells[IDX(x,y)];
}


void MineGame::print_grid() {
  int i,j;
  std::cout<<"  ";
  for(i=0; i<width; i++)
    std::cout<<(i%10)<<" ";
  std::cout<<"\n";

  for(j=0; j<height; j++) {
    std::cout<<(j%10)<<" ";
    for(i=0; i<width; i++) {
      //continue;
      if (FLAG(i,j))
	std::cout<<"X ";
      else if (HIDDEN(i,j))
	std::cout<<"- ";
      else if (GRID(i,j)==0)
	std::cout<<"  ";
      else if (BOMB(i,j))
	std::cout<<"B ";
      else
	std::cout<<GRID(i,j)<<" ";
    }
    std::cout<<" "<<j<<"\n";
  }
  std::cout<<"  ";
  for(i=0; i<width; i++)
    std::cout<<(i%10)<<" ";
  std::cout<<"\n";

  /*
  // print out debug
  std::cout<<"\nDEBUG:\n  ";
  for(i=0; i<width; i++)
    std::cout<<i<<"   ";
  std::cout<<"\n";

  for(j=0; j<height; j++) {
    std::cout<<j<<" ";
    for(i=0; i<width; i++) {
      if(knownCells[IDX(i,j)]==0){
	std::cout<<"    ";
	continue;
      }
      std::cout<<knownCells[IDX(i,j)]<<" ";
      if(knownCells[IDX(i,j)]>=0)
	std::cout<<" ";
      if(abs(knownCells[IDX(i,j)])<10)
	std::cout<<" ";
    }
    std::cout<<"\n";
  }
  //std::cout<<"  ";
  //for(i=0; i<width; i++)
  //std::cout<<i<<" ";
  std::cout<<"\n";
  */
}
