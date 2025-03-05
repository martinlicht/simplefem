#include <cstdio>   // Standard I/O operations (printf, etc.)
#include <cstdlib>  // Standard library (atexit, exit)
#include <ctime>    // Timing functions (clock, CLOCKS_PER_SEC)

#if defined(__linux__) or defined(__unix__)
#include <fcntl.h>   // File control options (fcntl, F_GETFL, F_SETFL, O_NONBLOCK)
#include <pthread.h> // POSIX multithreading (pthread_create, pthread_join)
#include <termios.h> // Terminal control (tcgetattr, tcsetattr)
#include <unistd.h>  // UNIX standard functions (read, sleep, STDIN_FILENO)
#endif 

#include "../../basic.hpp"

#if defined(__linux__) or defined(__unix__)

// Flag used for thread communication 
static bool process_reading_is_finished = false; 

// Stores the original_terminal_settings terminal settings before switching to raw mode.
static termios original_terminal_settings;

// Function prototypes
void terminal_mode_into_raw();
void terminal_mode_restore();
void* process_keyboard_reader( void* /*unused*/ );
void* process_printer( void* /*unused*/ );

/**
 * @brief Changes the terminal mode to non-canonical (raw) mode.
 * 
 * This function modifies terminal settings to disable line buffering and echoing.
 * It also ensures that settings can be restored when the program exits.
 */
void terminal_mode_into_raw()
{
    struct termios raw; 

    // Retrieve and store the current terminal settings for later restoring
    tcgetattr( STDIN_FILENO, &original_terminal_settings );
    
    // Retrieve and store the current terminal settings for modification
    tcgetattr( STDIN_FILENO, &raw );

    // Ensure that terminal settings are restored upon program exit
    atexit( &terminal_mode_restore );

    // Disable echo (ECHO) and canonical mode (ICANON)
    raw.c_lflag &= ~(ECHO | ICANON); 

    // Configure non-blocking read behavior
    raw.c_cc[VMIN]  = 1;  // Minimum number of characters to read
    raw.c_cc[VTIME] = 0;  // Timeout for read operations (0 = no timeout)

    // Apply new terminal settings
    tcsetattr( STDIN_FILENO, TCSAFLUSH, &raw );

    // Set non-blocking mode for stdin
    int fcntlFlags = fcntl(STDIN_FILENO, F_GETFL, 0);
    fcntl( STDIN_FILENO, F_SETFL, fcntlFlags | O_NONBLOCK );
}

/**
 * @brief Restores the original_terminal_settings terminal mode.
 * 
 * This function is automatically called when the program exits, ensuring that the terminal returns to its default settings.
 */
void terminal_mode_restore()
{
    tcsetattr( STDIN_FILENO, TCSAFLUSH, &original_terminal_settings );
}

/**
 * @brief Reads keyboard input in a separate thread.
 * 
 * This function listens for user keystrokes and prints them to the console.
 * If the ESC key (ASCII 27) is pressed, it terminates early.
 * 
 * @param unused (void*): Unused parameter to match pthread signature.
 * @return Always returns nullptr.
 */
void* process_keyboard_reader( void* /*unused*/ ) 
{
    // Switch terminal to raw mode for real-time input handling
    terminal_mode_into_raw(); 

    // Run for approximately 10 seconds
    const auto t_start = clock();
    while( clock() - t_start < 10 * CLOCKS_PER_SEC ) 
    {
        char c;
        int bytesRead = read( STDIN_FILENO, &c, 1 );

        if( bytesRead > 0 ) 
        {
            if( c == 27 ) { // ESC key detected
                printf("You pressed ESC. Quit.\n");
                break;
            } else {
                printf("You pressed %c\n", c);
            }
        }
    }

    // Signal that input reading is complete
    process_reading_is_finished = true;

    LOG << "Keystroke thread finished!\n";
  
    return nullptr;
}

/**
 * @brief Simulates background work while listening for keyboard input.
 * 
 * This function runs in a separate thread and prints time progress updates
 * every second while the keyboard input thread is active.
 * 
 * @param unused (void*): Unused parameter to match pthread signature.
 * @return Always returns nullptr.
 */
void* process_printer( void* /*unused*/ ) 
{
    const auto t_start = clock();
    
    while( not process_reading_is_finished and (clock() - t_start) < 10 * CLOCKS_PER_SEC ) 
    { 
        sleep(1);
        LOG << "Time: " << (clock() - t_start) << nl;
    }

    LOG << "Printer thread finished!\n";

    // Display memory allocation statistics (assumed function from basic.hpp)
    display_mallinfo(false);
  
    return nullptr;
}

#endif 

int main() 
{
#if defined(__linux__) or defined(__unix__)

    // Thread identifiers
    pthread_t thread_id_printer, thread_id_keyboardreader;

    // Create two threads: one for printing updates, one for reading keyboard input
    pthread_create( &thread_id_printer,        nullptr, process_printer,         nullptr );
    pthread_create( &thread_id_keyboardreader, nullptr, process_keyboard_reader, nullptr );

    // Wait for both threads to finish before exiting
    pthread_join( thread_id_printer,        nullptr );
    pthread_join( thread_id_keyboardreader, nullptr );

#endif 

    return 0;
}
