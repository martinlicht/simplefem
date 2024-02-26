#include <pthread.h> // Multithreading
#include <stdio.h>
#include <stdlib.h>  // for atexit()
#include <termios.h> // For changing terminal mode
#include <unistd.h>  // For changing terminal mode

#include "../../basic.hpp"


static termios original; // A struct to save the original state of terminal
static int ESCPressed = 0;      // For thread communication

void disableRAWMode();
void enableRAWMode();
void *asciRead( void* );
void *print( void* );

int main() {
  // Start Multithreading
  pthread_t id_print, id_read;

  pthread_create(&id_print, nullptr, print, nullptr);
  pthread_create(&id_read, nullptr, asciRead, nullptr);

  pthread_join(id_print, nullptr);
  pthread_join(id_read, nullptr);

  return 0;
}

/// Reads keyboard input
void *asciRead( void* ) {
  enableRAWMode(); // local function: Enable Raw Mode
  char ch;
  while ((ch = getchar()) != 27)
    ; // ASCI code for ESC is 27
  ESCPressed = 1;
  //printf("ESC Pressed!\n");
  LOG << "ESC Pressed!\n";
  return nullptr;
}

/// Doing Stuff while listening to keyboard
void *print( void* ) {
  while (!ESCPressed) { // When ESC is not pressed
    sleep(1);
    // printf("I am Printing!\n");
    LOG << "I am Printing!\n";
  }
  // printf("Printing Thread Finished!\n");
  LOG << "Printing Thread Finished!\n";
  display_mallinfo(false);
  return nullptr;
}

/// This function enables RAW mode for terminal.
void enableRAWMode() {
  struct termios raw;
  tcgetattr(STDIN_FILENO, &raw); // Save the state of the terminal to struct raw
                                 // STDIN_FILENO is from <stdlib.h>
                                 // tcgetattr() from <termios.h>
  tcgetattr(STDIN_FILENO, &original);
  atexit(&disableRAWMode); // Revert to canonical mode when exiting the program
                           // atext() from <stdlib.h>
  raw.c_lflag &=
      ~(ECHO | ICANON); // Turn off canonical mode
                        // Turn off ECHO mode so that keyboard is not
                        // printing to terminal
                        // ICANON and ECHO is bitflag. ~ is binary NOT operator

  tcsetattr(STDIN_FILENO, TCSAFLUSH,
            &raw); // Set the terminal to be in raw mode
                   // tcsetattr() from <termios.h>
}

void disableRAWMode() {
  tcsetattr(STDIN_FILENO, TCSAFLUSH,
            &original); // Set terminal to original state
}

