#include <stdlib.h>
#include <stdio.h>
//#include <getopt.h>
#include <string>
#include <iostream>
#include <sstream>
#include <string.h>

#include <fstream> // file io

//#include <ncurses.h>
//#include <locale.h>

#include "MineGame.h"

void play_text(MineGame game) {
  int x = 0, y = 0;

  std::string str, token;
  size_t p0 = 0, p1 = std::string::npos;
  int t = 0;
  bool flag;

  while(game.game_finished() == 0) {
    std::cout<<"\n";
    game.print_grid();
    /*
    std::cout<<"X: ";
    std::cin>>x;
    if(x<0)
      break;
    std::cout<<"Y: ";
    std::cin>>y;
    if(y<0)
      break;
    std::cout<<"\n";
    */
    std::cout<<"X Y (f): ";
    getline (std::cin, str);
    flag = false;

    p0 = 0;
    p1 = std::string::npos;
    t = 0;

    while(p0 != std::string::npos) {
      p1 = str.find_first_of(" ", p0);
      if(p1 != p0)
	{
	  token = str.substr(p0, p1 - p0);
	  if(t == 0)
	    x = atoi(token.c_str());
	  else if(t == 1) {
	    if(x>=0)
	      y = atoi(token.c_str());
	    else { // print to file
	      std::cout<<"Out file: "<<token<<"\n";
	      std::ofstream out_file;
	      out_file.open(token.c_str());

	      for(int j=0; j<game.getHeight(); j++) {
		for(int i=0; i<game.getWidth(); i++) {
		  if(game.known_val(i,j)==1)
		    out_file<<"-3";
		  else if(game.known_val(i,j)==2)
		    out_file<<"-4";
		  else if(game.flag(i,j))
		    out_file<<"-2";
		  else if(game.grid_val(i,j)<0)
		    out_file<<"-1";
		  else
		    out_file<<game.grid_val(i,j);

		  if(i!=(game.getWidth()-1))
		    out_file<<",";
		}
		if(j!=(game.getHeight()-1))
		  out_file<<"\n";
	      }

	      out_file.close();
	      //continue;
	    }
	  }
	  else
	    flag = true;
	  t++;
	}
      p0 = str.find_first_not_of(" ", p1);
    }

    //if(x<0) {
    //game.complete_write();
    //}
    if(x>=0)
      game.press_cell(x, y, flag);
  }
  std::cout<<"Game done, ";
  if(game.game_finished()==1)
    std::cout<<"You Won\n";
  else if(game.game_finished()==-1)
    std::cout<<"You Lost\n";
  else
    std::cout<<"it just stopped\n";
}

void play_gui(MineGame game) {
  /*
  WINDOW *win, *wrap;
 
  MEVENT evt;

  int tile_width = 20;
  int screen_height = game.getHeight()*tile_width;
  int screen_width = game.getWidth()*tile_width;

  initscr();

  win = newwin(screen_height + 2, 2 * screen_width + 2, 0, 0);
  wrap = newwin(screen_height + 3, 2 * screen_width + 2, 1, 0);
  
  keypad(wrap, 1);
  mousemask(BUTTON1_CLICKED | BUTTON2_CLICKED | BUTTON3_CLICKED, 0);

  while(true) {
    int ch;

    // repaint
    box(win, 0, 0);
    for(int j=0; j<game.getHeight(); j++) {
      for(int i=0; i<game.getWidth(); i++) {
	mvwprintw(win, i + 1, 2 * j + 1, " %c", 'X');
      }
    }
    wrefresh(wrap);
    wrefresh(win);
  }

  mousemask(0, 0);
  keypad(wrap, 0);
  endwin();
  */
}

int main(int argc, char** argv)
{
  int width = 10, height = 10, bombs = 15;
  bool display_time = false, display_knowns = false, display_debug_grids = false;

  for (int i = 1; i < argc; i++) { 
    //if (i + 1 != argc) 
    if (strcmp(argv[i], "-w") == 0) {//argv[i] == "-w") {
      if(i != (argc-1))
	width = atoi(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-h") == 0) {
      if(i != (argc-1))
	height = atoi(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-b") == 0) {
      if(i != (argc-1))
	bombs = atoi(argv[i + 1]);
      i++;
    } else if (strcmp(argv[i], "-t") == 0) {
      display_time = true;
    } else if (strcmp(argv[i], "-k") == 0) {
      display_knowns = true;
    } else if (strcmp(argv[i], "-d") == 0) {
      display_debug_grids = true;
    } else {
      std::cout << argv[i] << " Not valid.\n";
    }
    //std::cout << argv[i] << " ";
  }

  //MineGame game(10, 10, 15, true);
  MineGame game(width, height, bombs, true, 
		display_time, display_knowns, display_debug_grids);
  //play_gui(game);
  play_text(game);
  //delete(&game);
  return 0;
}
