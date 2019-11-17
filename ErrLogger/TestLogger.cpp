
/*
Test this with something like:
g++ -std=c++11 -x c++ -pthread -DLOGGING_LEVEL_ALL -DTEST_LOGGING logging.hpp -o logging_test
./logging_test
*/

#define LOGGING_LEVEL_ALL
#include "Logger.hh"

#include <thread>
#include <future>
#include <vector>
#include <functional>

size_t work() {
  std::ostringstream s; s << "hi my name is: " << std::this_thread::get_id();
  
  for(size_t i  = 0; i < 2; ++i) {
    //std::async is pretty uninteresting unless you make things yield
    logging::ERROR(s.str()); std::this_thread::sleep_for(std::chrono::milliseconds(10));
    logging::WARN(s.str()); std::this_thread::sleep_for(std::chrono::milliseconds(10));
    logging::INFO(s.str()); std::this_thread::sleep_for(std::chrono::milliseconds(10));
    logging::DEBUG(s.str()); std::this_thread::sleep_for(std::chrono::milliseconds(10));
    logging::TRACE(s.str()); std::this_thread::sleep_for(std::chrono::milliseconds(10));
    logging::log(logging::timestamp() + " \x1b[35;1m[CUSTOM]\x1b[0m " + s.str() + '\n'); 
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
  return 10;
}

int main(void) {
  //configure logging, if you dont it defaults to standard out logging with colors
  //logging::configure({ {"type", "file"}, {"file_name", "test.log"}, {"reopen_interval", "1"} });
  
  //start up some threads
  std::vector<std::future<size_t> > results;
  for(size_t i = 0; i < 4; ++i) {
    results.emplace_back(std::async(std::launch::async, work));
  }
  
  //dont really care about the results but we can pretend
  int exit_code = 0;
  for(auto& result : results) {
    try {
      size_t count = result.get();
    }
    catch(std::exception& e) {
      std::cout << e.what();
      exit_code++;
    }
  }
  return exit_code;
} // main
