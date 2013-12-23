
#include <iostream>
#include <fstream>
//#include <vector>

void saxpyCuda(int N, float alpha, float* x, float* y, float* result);
void printCudaInfo();

bool known_cells(int width, int height, int *grid, bool *flags, int *knowns, int *border, 
		 bool random_sol, bool display_time, bool display_knowns,
                 bool display_debug_grids);

class MineGame
{
 public:
  MineGame(int width, int height, int num_mines, bool write_game, bool display_time, 
	   bool display_knowns, bool display_debug_grids);
  ~MineGame();
  void run_game();
  void press_cell(int x, int y, bool flag);
  void print_grid();

  /**
   * Restarts the game by placing new bombs
   */
  void restart();

  /** 
   * Returns finished
   * 0 if not finished, 1 if won, -1 if lost
   */
  int game_finished();

  int getHeight();
  int getWidth();
  bool flag(int x, int y);
  int grid_val(int x, int y);
  int known_val(int x, int y);

  //void complete_write();

 private:
  //bool write_game;
  //bool has_game_to_write;
  //ofstream out_file;
  int width, height;
  bool first_move;

  /**
   * 0 if not finished, -1 if lost, 1 if won
   */
  int finished;
  int num_mines;

  bool display_time, display_knowns, display_debug_grids;

  // stores the number of incorrect flags
  int incorrect;

  /**
   * 1d array that will index into 2d grid
   * 0-8 indicate no bomb and have the number of adjacent bombs, shown
   * (-1)-(-9) indicate no bomb and is still hidden
   * -10 is a bomb, 10 is a displayed bomb
   */
  int *grid;

  /**
   * different array to store which cells have flags
   * true if flagged, false ow
   */
  bool *flags;

  //int *t;

  // this will be filled in by a GPU call, keep stored here so don't re-malloc everytime
  int *knownCells;

  // GPU will fill in this to identify borders
  //int *border;

  /**
   * places mines into the grid
   * assumes grid has been initialized in constructor
   */
  void setup_grid();

  /**
   *checks if all tiles that aren't bombs are unchecked, will modify fihished
   * returns true if won
   * NOTE: right now this is serial, mass check. It would be straightforward 
   *   to make this fairly cheap and not just call it on every move, maybe later
   */
  bool check_game_won();

  /**
   * This will show the cell at the given coordinates. 
   * If it is a bomb will trigger death. 
   * If 0 then will recursively call self to unveil larger area. 
   */
  void display_cell(int x, int y);

  /**
   * Set finished to -1, may do more in future?
   */
  void lose_game();

  /**
   * Fills in the numbers for given bomb placement
   * Called in restart() and if need to move a bomb
   */
  void calculate_neighbors();

};
