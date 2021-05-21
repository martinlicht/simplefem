
#if false 

// TODO introduce catch2 framework via header only inclusion 

#else

// Declare variables and functions for the tests without defining them

extern const char* test_topic;

extern int test_function( int argc, char* argv[] );

int main( int argc, char* argv[] )
{
	std::cout << "Running unit test....\n";
	std::cout << "\n";
	std::cout << "\t" << test_topic << "\n";
	std::cout << "\n";
	std::cout.flush();
	int ret = test_function( argc, argv );
	std::cout.flush();
	std::cout << "\n";
	std::cout << "Finished unit test!\n";
	std::cout.flush();
	return ret;
}

